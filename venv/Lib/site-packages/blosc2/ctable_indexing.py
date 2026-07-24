#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Indexing support mixed into :class:`blosc2.CTable`."""

from __future__ import annotations

import ast
import contextlib
import os
from typing import TYPE_CHECKING, Any

import numpy as np

import blosc2
from blosc2 import compute_chunks_blocks
from blosc2.schema import (
    DictionarySpec,
    ListSpec,
    NDArraySpec,
    ObjectSpec,
    StructSpec,
    Utf8Spec,
    VLBytesSpec,
    VLStringSpec,
)

if TYPE_CHECKING:
    from blosc2.ctable import CTable


class _FakeVlMeta:
    """Minimal vlmeta stand-in that accepts writes without touching a real SChunk."""

    def __init__(self):
        self._data: dict = {}

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def get(self, key, default=None):
        return self._data.get(key, default)


class _FakeSchunk:
    """Minimal SChunk stand-in whose vlmeta stores in memory."""

    def __init__(self):
        self.vlmeta = _FakeVlMeta()


def _dict_rank_hash(dictionary) -> str:
    """Stable hash of a dictionary's entries (code position + value).

    Used to detect when a rank-based FULL index has gone stale (the alphabetical
    ranks it encodes no longer match the live dictionary).  Must be stable across
    processes — the hash is persisted in the index descriptor and recomputed on a
    fresh open; ``hash()`` is PYTHONHASHSEED-salted and would spuriously mismatch,
    making the persistent index always fall back.
    """
    import hashlib

    h = hashlib.sha1(usedforsecurity=False)
    for i, value in enumerate(dictionary):
        h.update(repr((i, value)).encode("utf-8"))
    return h.hexdigest()


class _DictRankWrapper:
    """Wrap a dictionary column's codes NDArray, translating codes to alphabetical ranks on read.

    Mirrors the NDArray interface enough for the index builder (dtype, shape, ndim,
    chunks, blocks, __getitem__).  The rank for code *c* is its position in an
    alphabetically-sorted dictionary, so sorting by rank == sorting by decoded string.
    Null codes map to a sentinel rank ``null_rank = len(dictionary)`` (largest, so
    nulls sort last).
    """

    def __init__(
        self,
        codes,
        code_to_rank: np.ndarray,
        null_rank: np.int32,
        null_code: int,
        nullable: bool,
        n_live: int,
    ):
        self._codes = codes  # int32 NDArray
        self._code_to_rank = code_to_rank  # int32 array mapping code -> rank
        self._null_rank = null_rank
        self._null_code = null_code
        self._nullable = nullable
        self.dtype = np.dtype(np.int32)
        # The codes array carries capacity padding beyond the live rows; expose only
        # the live range so the index sidecars match n_rows (no padding → the
        # zero-permutation window read engages instead of falling back).
        self.shape = (n_live,)
        self.ndim = 1
        chunk0 = codes.chunks[0] if codes.chunks else n_live
        block0 = codes.blocks[0] if codes.blocks else n_live
        self.chunks = (min(chunk0, n_live),)
        self.blocks = (min(block0, n_live),)

    def __getitem__(self, key):
        codes_slice = np.asarray(self._codes[key], dtype=np.int32)
        ranks = self._code_to_rank[codes_slice]
        if self._nullable:
            ranks[codes_slice == self._null_code] = self._null_rank
        return ranks


class _CTableBuildProxy:
    """Minimal shim that lets the ``indexing`` module build sidecars for a
    CTable column without touching the column's own ``schunk.vlmeta``.

    Attributes mirror those required by the internal build functions:
    ``urlpath``, ``schunk``, ``shape``, ``ndim``, ``dtype``, ``chunks``,
    ``blocks``, and item access via ``__getitem__``.
    """

    def __init__(self, col_array: blosc2.NDArray, anchor_urlpath: str | None) -> None:
        self._col_array = col_array
        self.urlpath = anchor_urlpath  # controls sidecar placement
        self.schunk = _FakeSchunk()
        self.shape = col_array.shape
        self.ndim = col_array.ndim
        self.dtype = col_array.dtype
        self.chunks = col_array.chunks
        self.blocks = col_array.blocks

    def __getitem__(self, key):
        return self._col_array[key]


class _CTableIndexingMixin:
    # Cost-model constants for cross-column index refinement.
    # Calibrated from profiling with sparse-gather optimisations.
    #   _GATHER_COST_MS_PER_1K_ITEMS_PER_OP  ≈ ms to sparse-gather 1000 items from one operand column
    #   _SCAN_COST_MS_PER_1M_ROWS             ≈ ms to miniexpr-scan 1 million rows
    # If refinement cost exceeds scan cost, fall back to a full scan.
    _GATHER_COST_MS_PER_1K_ITEMS_PER_OP: float = 3.5
    _SCAN_COST_MS_PER_1M_ROWS: float = 4.3

    # Cost-model constants for the segment-summary (SUMMARY index) gate.
    # Surviving blocks are evaluated through miniexpr with a candidate-block
    # bitmap (non-candidate blocks are skipped inside the prefilter), so the
    # per-block overhead is small; what matters is that, once enough blocks
    # survive, evaluating them (plus iterating every chunk) costs more than a
    # plain sequential scan, and the planner falls back.
    #   _SEGMENT_SCAN_MS_PER_1M_ROWS_PER_OP ≈ ms to sequentially scan 1M rows of one column
    #   _SEGMENT_OVERHEAD_MS_PER_OP         ≈ fixed ms per surviving block per operand column
    _SEGMENT_SCAN_MS_PER_1M_ROWS_PER_OP: float = 1.6
    _SEGMENT_OVERHEAD_MS_PER_OP: float = 0.2

    @property
    def _root_table(self) -> CTable:
        """Return the root (non-view) table; *self* if not a view."""
        t = self
        while t.base is not None:
            t = t.base
        return t

    def _invalidate_index_catalog_cache(self) -> None:
        root = self._root_table
        root._cached_index_catalog = None
        root._cached_index_catalog_revision = None

    def _get_index_catalog(self) -> dict:
        root = self._root_table
        revision = root._storage.index_catalog_revision()
        catalog = getattr(root, "_cached_index_catalog", None)
        if catalog is None or getattr(root, "_cached_index_catalog_revision", None) != revision:
            catalog = root._storage.load_index_catalog()
            root._cached_index_catalog = catalog
            root._cached_index_catalog_revision = revision
        return catalog

    def _mark_all_indexes_stale(self) -> None:
        """Bump value_epoch and mark every catalog entry stale on the root table."""
        root = self._root_table
        root._storage.bump_value_epoch()
        catalog = root._get_index_catalog()
        if not catalog:
            return
        changed = False
        for desc in catalog.values():
            if not desc.get("stale", False):
                desc["stale"] = True
                changed = True
        if changed:
            root._storage.save_index_catalog(catalog)
            root._invalidate_index_catalog_cache()

    @staticmethod
    def _validate_index_descriptor(col_name: str, descriptor: dict) -> None:
        """Raise ValueError when an index catalog entry is malformed."""
        if not isinstance(descriptor, dict):
            raise ValueError(f"Malformed index metadata for column {col_name!r}: descriptor must be a dict.")
        token = descriptor.get("token")
        if not isinstance(token, str) or not token:
            raise ValueError(f"Malformed index metadata for column {col_name!r}: missing token.")
        kind = descriptor.get("kind")
        if kind not in {"summary", "bucket", "partial", "full", "opsi"}:
            raise ValueError(f"Malformed index metadata for column {col_name!r}: invalid kind {kind!r}.")
        if kind == "bucket" and not isinstance(descriptor.get("bucket"), dict):
            raise ValueError(f"Malformed index metadata for column {col_name!r}: missing bucket payload.")
        if kind == "partial" and not isinstance(descriptor.get("partial"), dict):
            raise ValueError(f"Malformed index metadata for column {col_name!r}: missing partial payload.")
        if kind == "full" and not isinstance(descriptor.get("full"), dict):
            raise ValueError(f"Malformed index metadata for column {col_name!r}: missing full payload.")

    def _drop_index_descriptor(self, col_name: str, descriptor: dict) -> None:
        """Delete sidecars/cache for a catalog descriptor without touching the column mapping."""
        from pathlib import Path

        from blosc2.indexing import (
            _IN_MEMORY_INDEXES,
            _PERSISTENT_INDEXES,
            _array_key,
            _clear_cached_data,
            _drop_descriptor_sidecars,
            _is_persistent_array,
        )

        token = descriptor["token"]
        col_arr = None
        with contextlib.suppress(Exception):
            col_arr = self._index_target_array(col_name, descriptor)

        if col_arr is not None:
            _clear_cached_data(col_arr, token)

        if col_arr is not None and _is_persistent_array(col_arr):
            arr_key = _array_key(col_arr)
            store = _PERSISTENT_INDEXES.get(arr_key)
            if store is not None:
                store["indexes"].pop(token, None)
        elif col_arr is not None:
            store = _IN_MEMORY_INDEXES.get(id(col_arr))
            if store is not None:
                store["indexes"].pop(token, None)

        _drop_descriptor_sidecars(descriptor)
        self._root_table._expr_index_arrays.pop(token, None)

        expr_values_path = descriptor.get("expr_values_path")
        if expr_values_path is not None:
            with contextlib.suppress(OSError):
                os.remove(expr_values_path)

        anchor = self._storage.index_anchor_path(col_name)
        if anchor is not None:
            proxy_key = ("persistent", str(Path(anchor).resolve()))
            _PERSISTENT_INDEXES.pop(proxy_key, None)
            with contextlib.suppress(OSError):
                os.rmdir(os.path.dirname(anchor))

    def _index_create_kwargs_from_descriptor(self, descriptor: dict) -> dict[str, Any]:
        """Return create_index kwargs that rebuild an existing descriptor."""
        build = "ooc" if bool(descriptor.get("ooc", False)) else "memory"
        kwargs = {
            "kind": descriptor["kind"],
            "optlevel": int(descriptor.get("optlevel", 5)),
            "name": descriptor.get("name") or None,
            "build": build,
            "cparams": descriptor.get("cparams"),
        }
        if descriptor.get("kind") == "full":
            kwargs["method"] = descriptor.get("full", {}).get("build_method", "global-sort")
        if descriptor.get("kind") == "opsi":
            kwargs["opsi_max_cycles"] = descriptor.get("opsi", {}).get("max_cycles")
        target = descriptor.get("target") or {}
        if target.get("source") == "expression":
            kwargs["expression"] = target.get("expression")
        return kwargs

    def _normalize_table_expression_target(
        self, expression: str, operands: dict | None = None
    ) -> tuple[dict, np.dtype]:
        """Normalize a same-table expression target and infer its dtype."""
        if operands is None:
            operands = self._cols
        try:
            ast.parse(expression, mode="eval")
        except SyntaxError as exc:
            raise ValueError("expression is not valid Python syntax") from exc

        owned_ids = {id(arr): name for name, arr in self._root_table._cols.items()}
        dependencies: list[str] = []
        valid = True

        class _Canonicalizer(ast.NodeTransformer):
            def visit_Name(self_inner, node: ast.Name) -> ast.AST:
                nonlocal valid
                operand = operands.get(node.id)
                if operand is None or not isinstance(operand, blosc2.NDArray):
                    return node
                cname = owned_ids.get(id(operand))
                if cname is None:
                    valid = False
                    return node
                dependencies.append(cname)
                return ast.copy_location(ast.Name(id=cname, ctx=node.ctx), node)

        normalized = _Canonicalizer().visit(
            ast.fix_missing_locations(ast.parse(expression, mode="eval")).body
        )
        if not valid or not dependencies:
            raise ValueError("expression indexes require operands from stored columns of the same table")
        dependencies = list(dict.fromkeys(dependencies))
        expression_key = ast.unparse(normalized)
        lazy = blosc2.lazyexpr(expression_key, {dep: self._root_table._cols[dep] for dep in dependencies})
        sample_stop = min(
            len(self._root_table._valid_rows), max(1, int(self._root_table._valid_rows.blocks[0]))
        )
        sample = lazy[:sample_stop]
        if isinstance(sample, blosc2.NDArray):
            sample = sample[:]
        sample = np.asarray(sample)
        dtype = np.dtype(sample.dtype)
        if sample.ndim != 1:
            raise ValueError("expression indexes require expressions returning a 1-D scalar stream")
        target = {
            "source": "expression",
            "expression": expression,
            "expression_key": expression_key,
            "dependencies": dependencies,
        }
        return target, dtype

    def _expression_index_values_path(self, token: str) -> str | None:
        anchor = self._storage.index_anchor_path(token)
        if anchor is None:
            return None
        return os.path.join(os.path.dirname(anchor), "values.b2nd")

    def _build_expression_values_array(self, target: dict, dtype: np.dtype, cparams=None) -> blosc2.NDArray:
        """Build a physical 1-D values array for a table expression target."""
        from blosc2.indexing import _target_token

        root = self._root_table
        capacity = len(root._valid_rows)
        chunks, blocks = compute_chunks_blocks((capacity,), dtype=dtype)
        urlpath = root._expression_index_values_path(_target_token(target))
        if urlpath is not None:
            os.makedirs(os.path.dirname(urlpath), exist_ok=True)
            arr = blosc2.zeros(
                (capacity,), dtype=dtype, urlpath=urlpath, mode="w", chunks=chunks, blocks=blocks
            )
        else:
            arr = blosc2.zeros((capacity,), dtype=dtype, chunks=chunks, blocks=blocks)
        lazy = blosc2.lazyexpr(
            target["expression_key"], {dep: root._cols[dep] for dep in target["dependencies"]}
        )
        step = int(root._valid_rows.chunks[0]) if root._valid_rows.chunks else 65536
        for start in range(0, capacity, step):
            stop = min(start + step, capacity)
            values = lazy[start:stop]
            if isinstance(values, blosc2.NDArray):
                values = values[:]
            arr[start:stop] = np.asarray(values, dtype=dtype)
        root._expr_index_arrays[_target_token(target)] = arr
        return arr

    def _index_target_array(self, lookup_key: str, descriptor: dict) -> blosc2.NDArray:
        """Return the physical array backing a column or expression index."""
        target = descriptor.get("target") or {}
        if target.get("source") != "expression":
            return self._root_table._cols[lookup_key]
        token = descriptor["token"]
        root = self._root_table
        arr = root._expr_index_arrays.get(token)
        if arr is not None:
            return arr
        path = descriptor.get("expr_values_path")
        if path is None:
            raise KeyError(f"No backing array found for expression index {token!r}.")
        arr = blosc2.open(path, mode="r" if root._read_only else "a")
        root._expr_index_arrays[token] = arr
        return arr

    def _resolve_index_catalog_entry(
        self, col_name: str | None = None, *, expression: str | None = None, name: str | None = None
    ) -> tuple[str, dict]:
        """Resolve an index catalog entry by column, expression, or label."""
        catalog = self._root_table._get_index_catalog()
        if col_name is not None and expression is not None:
            raise ValueError("col_name and expression are mutually exclusive")
        if col_name is not None:
            col_name = self._logical_to_physical_name(col_name)
            if col_name not in catalog:
                raise KeyError(f"No index found for column {col_name!r}.")
            return col_name, catalog[col_name]
        if expression is not None:
            from blosc2.indexing import _target_token

            target, _ = self._normalize_table_expression_target(expression)
            token = _target_token(target)
            if token not in catalog:
                raise KeyError(f"No index found for expression {expression!r}.")
            return token, catalog[token]
        if name is not None:
            matches = [(key, desc) for key, desc in catalog.items() if desc.get("name") == name]
            if not matches:
                raise KeyError(f"No index found with name {name!r}.")
            if len(matches) > 1:
                raise ValueError(f"Multiple indexes found with name {name!r}; specify a target explicitly.")
            return matches[0]
        raise TypeError("must specify col_name, expression, or name")

    def _build_index_persistent(
        self,
        col_name: str,
        col_arr: blosc2.NDArray,
        *,
        kind: str,
        optlevel: int,
        name_hint: str | None,
        build: str,
        tmpdir: str | None,
        cparams_obj,
        method: str | None = None,
        opsi_max_cycles: int | None = None,
        summary_levels: tuple[str, ...] | None = None,
        precomputed_summaries: dict | None = None,
    ) -> dict:
        """Build index sidecar files for a persistent-table column; return the descriptor."""
        import tempfile
        from pathlib import Path

        from blosc2.indexing import (
            _PERSISTENT_INDEXES,
            _array_key,
            _build_bucket_descriptor,
            _build_bucket_descriptor_ooc,
            _build_descriptor,
            _build_full_descriptor,
            _build_full_descriptor_ooc,
            _build_levels_descriptor,
            _build_levels_descriptor_ooc,
            _build_opsi_descriptor,
            _build_partial_descriptor,
            _build_partial_descriptor_ooc,
            _copy_descriptor,
            _field_target_descriptor,
            _resolve_full_index_tmpdir,
            _resolve_ooc_mode,
            _target_token,
            _values_for_target,
        )

        anchor = self._storage.index_anchor_path(col_name)
        os.makedirs(os.path.dirname(anchor), exist_ok=True)
        proxy = _CTableBuildProxy(col_arr, anchor)
        proxy_key = _array_key(proxy)
        _PERSISTENT_INDEXES.pop(proxy_key, None)  # clear any stale cache entry

        target = _field_target_descriptor(None)
        token = _target_token(target)
        persistent = True
        dtype = col_arr.dtype
        use_ooc = _resolve_ooc_mode(kind, build)
        if opsi_max_cycles is None:
            opsi_max_cycles = max(1, optlevel if optlevel < 8 else optlevel * 2)

        if use_ooc:
            resolved_tmpdir = _resolve_full_index_tmpdir(proxy, tmpdir)
            levels = _build_levels_descriptor_ooc(
                proxy,
                target,
                token,
                kind,
                dtype,
                persistent,
                cparams_obj,
                summary_levels=summary_levels,
                precomputed_summaries=precomputed_summaries,
            )
            bucket = (
                _build_bucket_descriptor_ooc(
                    proxy, target, token, kind, dtype, optlevel, persistent, cparams_obj
                )
                if kind == "bucket"
                else None
            )
            partial = (
                _build_partial_descriptor_ooc(
                    proxy, target, token, kind, dtype, optlevel, persistent, cparams_obj
                )
                if kind == "partial"
                else None
            )
            full = None
            opsi = None
            if kind == "full":
                with tempfile.TemporaryDirectory(prefix="blosc2-index-ooc-", dir=resolved_tmpdir) as td:
                    full = _build_full_descriptor_ooc(
                        proxy, target, token, kind, dtype, persistent, Path(td), cparams_obj, optlevel
                    )
                    full["build_method"] = "global-sort"
            if kind == "opsi":
                opsi = _build_opsi_descriptor(
                    proxy, target, token, kind, dtype, persistent, cparams_obj, opsi_max_cycles, optlevel
                )
            descriptor = _build_descriptor(
                proxy,
                target,
                token,
                kind,
                optlevel,
                persistent,
                True,
                name_hint,
                dtype,
                levels,
                bucket,
                partial,
                full,
                cparams_obj,
                opsi,
            )
        else:
            values = _values_for_target(proxy, target)
            levels = _build_levels_descriptor(
                proxy,
                target,
                token,
                kind,
                dtype,
                values,
                persistent,
                cparams_obj,
                summary_levels=summary_levels,
            )
            bucket = (
                _build_bucket_descriptor(proxy, token, kind, values, optlevel, persistent, cparams_obj)
                if kind == "bucket"
                else None
            )
            partial = (
                _build_partial_descriptor(proxy, token, kind, values, optlevel, persistent, cparams_obj)
                if kind == "partial"
                else None
            )
            full = None
            opsi = None
            if kind == "full":
                full = _build_full_descriptor(proxy, token, kind, values, persistent, cparams_obj, optlevel)
                full["build_method"] = "global-sort"
            if kind == "opsi":
                opsi = _build_opsi_descriptor(
                    proxy, target, token, kind, dtype, persistent, cparams_obj, opsi_max_cycles, optlevel
                )
            descriptor = _build_descriptor(
                proxy,
                target,
                token,
                kind,
                optlevel,
                persistent,
                False,
                name_hint,
                dtype,
                levels,
                bucket,
                partial,
                full,
                cparams_obj,
                opsi,
            )

        result = _copy_descriptor(descriptor)
        _PERSISTENT_INDEXES.pop(proxy_key, None)  # evict proxy to avoid memory leak
        return result

    def create_index(  # noqa: C901
        self,
        col_name: str | None = None,
        *,
        field: str | None = None,
        expression: str | None = None,
        operands: dict | None = None,
        kind: blosc2.IndexKind = blosc2.IndexKind.BUCKET,
        optlevel: int = 5,
        name: str | None = None,
        build: str = "auto",
        tmpdir: str | None = None,
        granularity: str | None = None,
        **kwargs,
    ) -> blosc2.Index:
        """Build and register an index for a stored column or table expression.

        For tables with **nested (dotted) column names**, pass the dotted leaf
        name directly::

            t.create_index("trip.begin.lon")
            t.where("trip.begin.lon > -87.7").nrows   # index is used automatically

        .. rubric:: Choosing an index kind

        ``BUCKET`` (the default) is the cheapest to build and store.
        It accelerates single‑column ``where`` queries and ``sort_by``
        reuse with approximate ordering derived from value
        quantization.  Sufficient for most workloads.

        ``FULL`` builds a globally sorted index that returns exact
        row positions for any range predicate.  It enables the
        **cross‑column refinement** planner path: when a multi‑column
        conjunction such as ``(tips > 100) & (km > 0) & (sec > 0)``
        indexes only the most selective column, the planner obtains
        compact exact positions from ``FULL`` and evaluates the
        remaining predicates on just those rows.  ``FULL`` is also
        ideal for ``sort_by`` reuse because it carries a complete
        sort order.

        ``PARTIAL`` builds a chunk‑local sorted payload with segment
        navigation.  It is cheaper to build than ``FULL`` (roughly
        half the raw storage) while still providing exact positions
        for cross‑column refinement.  Its exact positions are most
        compact for equality or narrow range queries; wide ranges
        may scan proportionally more candidate segments.

        ``OPSI`` is a specialised tier for approximate ordering;
        prefer ``FULL`` when a globally sorted ordered index is
        needed to accelerate ``sort_by``.

        ``SUMMARY`` stores only per‑segment min/max and is the
        lightest kind; it may still skip segments for broad range
        queries but cannot accelerate ``sort_by``.

        .. rubric:: SUMMARY granularity

        For ``kind=SUMMARY``, ``granularity`` controls the segment size of the
        min/max summaries: ``"block"`` (the default) summarizes at Blosc2 block
        granularity; ``"chunk"`` uses the coarser chunk granularity (fewer
        summary entries, cheaper to build, less selective pruning).
        ``granularity`` is only valid for ``kind=SUMMARY``.
        """
        if self.base is not None:
            raise ValueError("Cannot create an index on a view.")
        if col_name is not None and field is not None:
            raise ValueError("col_name and field are mutually exclusive")
        if expression is not None and (col_name is not None or field is not None):
            raise ValueError("column targets and expression are mutually exclusive")
        if operands is not None and expression is None:
            raise ValueError("operands can only be provided together with expression")
        col_name = field if field is not None else col_name
        if col_name is not None:
            col_name = self._logical_to_physical_name(col_name)

        from blosc2.indexing import (
            _IN_MEMORY_INDEXES,
            _copy_descriptor,
            _normalize_build_mode,
            _normalize_full_build_method,
            _normalize_index_cparams,
            _normalize_index_kind,
            _resolve_summary_levels,
            _target_token,
        )
        from blosc2.indexing import create_index as _ix_create_index

        cparams_obj = _normalize_index_cparams(kwargs.pop("cparams", None))
        method = kwargs.pop("method", None)
        opsi_max_cycles = kwargs.pop("opsi_max_cycles", None)
        precomputed_summaries = kwargs.pop("precomputed_summaries", None)
        if opsi_max_cycles is not None:
            opsi_max_cycles = max(1, int(opsi_max_cycles))
        if kwargs:
            raise TypeError(f"unexpected keyword argument(s): {', '.join(sorted(kwargs))}")

        kind_str = _normalize_index_kind(kind)
        build_str = _normalize_build_mode(build)
        method_str = _normalize_full_build_method(method) if kind_str == "full" else None
        if method is not None and kind_str != "full":
            raise ValueError("method is only supported for kind=IndexKind.FULL")
        if granularity is not None and kind_str != "summary":
            raise ValueError("granularity is only supported for kind=IndexKind.SUMMARY")
        summary_levels = _resolve_summary_levels(granularity)
        catalog = self._get_index_catalog()

        if expression is not None:
            target, dtype = self._normalize_table_expression_target(expression, operands)
            token = _target_token(target)
            if token in catalog:
                raise ValueError(
                    f"Index already exists for expression {expression!r}. "
                    "Call rebuild_index() to replace it or drop_index() first."
                )
            expr_arr = self._build_expression_values_array(target, dtype, cparams=cparams_obj)
            _ix_create_index(
                expr_arr,
                kind=blosc2.IndexKind(kind_str),
                optlevel=optlevel,
                name=name,
                build=build,
                tmpdir=tmpdir,
                cparams=cparams_obj,
                method=method_str,
                opsi_max_cycles=opsi_max_cycles,
                summary_levels=summary_levels,
            )
            store = _IN_MEMORY_INDEXES.get(id(expr_arr))
            if store is None:
                from blosc2.indexing import _load_store

                store = _load_store(expr_arr)
            descriptor = _copy_descriptor(store["indexes"]["__self__"])
            descriptor["target"] = target
            descriptor["token"] = token
            descriptor["dtype"] = str(np.dtype(dtype))
            descriptor["expr_values_path"] = getattr(expr_arr, "urlpath", None)
            value_epoch, _ = self._storage.get_epoch_counters()
            descriptor["built_value_epoch"] = value_epoch
            catalog[token] = descriptor
            self._storage.save_index_catalog(catalog)
            self._invalidate_index_catalog_cache()
            return blosc2.Index._from_table(self, token, descriptor)

        if col_name is None:
            raise TypeError("must specify col_name/field or expression")
        if col_name in self._computed_cols:
            raise ValueError(
                f"Cannot create an index on computed column {col_name!r}: "
                "computed columns have no physical storage."
            )
        if col_name not in self._cols:
            raise KeyError(f"No column named {col_name!r}. Available: {self.col_names}")
        self._ensure_generated_column_not_stale(col_name)
        if col_name in catalog:
            raise ValueError(
                f"Index already exists for column {col_name!r}. "
                "Call rebuild_index() to replace it or drop_index() first."
            )

        col_arr = self._cols[col_name]
        if isinstance(self._schema.columns_by_name[col_name].spec, NDArraySpec):
            spec = self._schema.columns_by_name[col_name].spec
            raise ValueError(
                f"Cannot create an index on ndarray column {col_name!r} with per-row shape {spec.item_shape}. "
                "Materialize a scalar generated column first, e.g. embedding_norm or embedding_max."
            )
        if isinstance(self._schema.columns_by_name[col_name].spec, ListSpec):
            raise ValueError(f"Cannot create an index on list column {col_name!r} in V1.")
        if isinstance(self._schema.columns_by_name[col_name].spec, Utf8Spec):
            raise NotImplementedError(
                f"Cannot create an index on variable-length utf8 column {col_name!r}: "
                "indexing for utf8 columns is not supported yet. "
                "Use a fixed-width string(max_length=N) column if you need an index."
            )
        if isinstance(
            self._schema.columns_by_name[col_name].spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec)
        ):
            raise NotImplementedError(
                f"Cannot create an index on variable-length scalar column {col_name!r}: "
                "indexing for vlstring/vlbytes/struct/object columns is not supported yet."
            )
        # Dictionary columns: index by alphabetical rank instead of insertion-order codes.
        is_dictionary = isinstance(self._schema.columns_by_name[col_name].spec, DictionarySpec)
        dict_rank_meta = None
        if is_dictionary:
            dict_col = col_arr
            n_live = self._n_rows if self._n_rows is not None else len(self._valid_rows)
            dictionary = list(dict_col.dictionary)
            n_entries = len(dictionary)
            order = np.argsort(dictionary, kind="stable")
            code_to_rank = np.empty(n_entries, dtype=np.int32)
            code_to_rank[order] = np.arange(n_entries, dtype=np.int32)
            null_code = dict_col.spec.null_code
            null_rank = np.int32(n_entries)
            # Hash for staleness detection.
            dict_hash = _dict_rank_hash(dictionary)
            dict_rank_meta = {"null_rank": int(null_rank), "dict_hash": dict_hash, "dict_len": n_entries}
            col_arr = _DictRankWrapper(
                dict_col.codes, code_to_rank, null_rank, null_code, dict_col.spec.nullable, n_live
            )
        is_persistent = self._storage.index_anchor_path(col_name) is not None

        if is_persistent:
            descriptor = self._build_index_persistent(
                col_name,
                col_arr,
                kind=kind_str,
                optlevel=optlevel,
                name_hint=name,
                build=build_str,
                tmpdir=tmpdir,
                cparams_obj=cparams_obj,
                method=method_str,
                opsi_max_cycles=opsi_max_cycles,
                summary_levels=summary_levels,
                precomputed_summaries=precomputed_summaries if kind_str == "summary" else None,
            )
        else:
            # In-memory path: materialise ranks as a proper NDArray (small tables only).
            if is_dictionary:
                codes = np.asarray(dict_col.codes[:n_live], dtype=np.int32)
                ranks_arr = code_to_rank[codes]
                if dict_col.spec.nullable:
                    ranks_arr[codes == null_code] = null_rank
                col_arr = blosc2.asarray(ranks_arr)
            _ix_create_index(
                col_arr,
                field=None,
                kind=blosc2.IndexKind(kind_str),
                optlevel=optlevel,
                name=name,
                build=build,
                tmpdir=tmpdir,
                cparams=cparams_obj,
                method=method_str,
                opsi_max_cycles=opsi_max_cycles,
                summary_levels=summary_levels,
                precomputed_summaries=precomputed_summaries if kind_str == "summary" else None,
            )
            store = _IN_MEMORY_INDEXES[id(col_arr)]
            descriptor = _copy_descriptor(store["indexes"]["__self__"])
        if dict_rank_meta is not None:
            full = descriptor.setdefault("full", {})
            if full is None:
                full = descriptor["full"] = {}
            full["dict_rank"] = dict_rank_meta

        value_epoch, _ = self._storage.get_epoch_counters()
        descriptor["built_value_epoch"] = value_epoch

        if is_persistent:
            # Use column name as token so sibling columns in compact stores get
            # distinct keys in the shared _PERSISTENT_INDEXES store (all columns
            # share the same urlpath yet must not collide on "__self__").  File
            # paths were built with token="__self__" (omitted from the filename)
            # and remain unchanged; only the in-memory key is affected.
            # In-memory CTables already have unique _IN_MEMORY_INDEXES entries
            # per column (keyed by id(array)), so no token change is needed.
            descriptor["token"] = col_name

        catalog = self._get_index_catalog()
        catalog[col_name] = descriptor
        self._storage.save_index_catalog(catalog)
        self._invalidate_index_catalog_cache()
        return blosc2.Index._from_table(self, col_name, descriptor)

    def drop_index(
        self, col_name: str | None = None, *, expression: str | None = None, name: str | None = None
    ) -> None:
        """Remove an index and delete any sidecar files."""
        if self.base is not None:
            raise ValueError("Cannot drop an index from a view.")

        lookup_key, descriptor = self._resolve_index_catalog_entry(
            col_name, expression=expression, name=name
        )
        catalog = self._get_index_catalog()
        catalog.pop(lookup_key, None)
        self._validate_index_descriptor(lookup_key, descriptor)
        self._drop_index_descriptor(lookup_key, descriptor)
        self._storage.save_index_catalog(catalog)
        self._invalidate_index_catalog_cache()

    def rebuild_index(
        self, col_name: str | None = None, *, expression: str | None = None, name: str | None = None
    ) -> blosc2.Index:
        """Drop and recreate an index with the same parameters."""
        if self.base is not None:
            raise ValueError("Cannot rebuild an index on a view.")

        lookup_key, old_desc = self._resolve_index_catalog_entry(col_name, expression=expression, name=name)
        self._validate_index_descriptor(lookup_key, old_desc)
        create_kwargs = self._index_create_kwargs_from_descriptor(old_desc)

        self.drop_index(col_name, expression=expression, name=name)
        if "expression" in create_kwargs:
            return self.create_index(expression=create_kwargs.pop("expression"), **create_kwargs)
        return self.create_index(lookup_key, **create_kwargs)

    def compact_index(
        self, col_name: str | None = None, *, expression: str | None = None, name: str | None = None
    ) -> blosc2.Index:
        """Compact an index, merging any incremental append runs."""
        if self.base is not None:
            raise ValueError("Cannot compact an index on a view.")

        from blosc2.indexing import (
            _IN_MEMORY_INDEXES,
            _PERSISTENT_INDEXES,
            _array_key,
            _copy_descriptor,
            _default_index_store,
            _is_persistent_array,
        )
        from blosc2.indexing import compact_index as _ix_compact_index

        lookup_key, descriptor = self._resolve_index_catalog_entry(
            col_name, expression=expression, name=name
        )
        col_arr = self._index_target_array(lookup_key, descriptor)
        catalog = self._get_index_catalog()

        if _is_persistent_array(col_arr):
            anchor = self._storage.index_anchor_path(lookup_key)
            proxy = _CTableBuildProxy(col_arr, anchor)
            proxy_key = _array_key(proxy)
            store = _default_index_store()
            store["indexes"][descriptor["token"]] = descriptor
            _PERSISTENT_INDEXES[proxy_key] = store
            try:
                _ix_compact_index(proxy)
                updated_store = _PERSISTENT_INDEXES.get(proxy_key) or store
                updated_desc = _copy_descriptor(updated_store["indexes"][descriptor["token"]])
            finally:
                _PERSISTENT_INDEXES.pop(proxy_key, None)
            updated_desc["built_value_epoch"] = descriptor.get("built_value_epoch", 0)
            catalog[lookup_key] = updated_desc
            self._storage.save_index_catalog(catalog)
            self._invalidate_index_catalog_cache()
            return blosc2.Index._from_table(self, lookup_key, updated_desc)
        else:
            _ix_compact_index(col_arr)
            store = _IN_MEMORY_INDEXES.get(id(col_arr))
            if store:
                token = descriptor["token"]
                updated_desc = _copy_descriptor(store["indexes"].get(token, descriptor))
                updated_desc["built_value_epoch"] = descriptor.get("built_value_epoch", 0)
                catalog[lookup_key] = updated_desc
                self._storage.save_index_catalog(catalog)
                self._invalidate_index_catalog_cache()
                return blosc2.Index._from_table(self, lookup_key, updated_desc)
            return blosc2.Index._from_table(self, lookup_key, descriptor)

    def index(
        self, col_name: str | None = None, *, expression: str | None = None, name: str | None = None
    ) -> blosc2.Index:
        """Return the index handle for a stored-column or expression target."""
        lookup_key, descriptor = self._resolve_index_catalog_entry(
            col_name, expression=expression, name=name
        )
        return blosc2.Index._from_table(self, lookup_key, descriptor)

    @property
    def indexes(self) -> list[blosc2.Index]:
        """Return a list of :class:`blosc2.Index` handles for all active indexes."""
        catalog = self._root_table._get_index_catalog()
        return [blosc2.Index._from_table(self, col_name, desc) for col_name, desc in catalog.items()]

    def _rewrite_expression_query_for_index(
        self, expression: str, operands: dict, target: dict
    ) -> str | None:
        """Rewrite matching table-expression subtrees to ``_where_x`` for planning."""
        try:
            tree = ast.parse(expression, mode="eval")
        except SyntaxError:
            return None

        class _Rewriter(ast.NodeTransformer):
            def __init__(self, outer):
                self.outer = outer
                self.changed = False

            def generic_visit(self, node):
                normalized = None
                with contextlib.suppress(Exception):
                    normalized, _ = self.outer._normalize_table_expression_target(
                        ast.unparse(node), operands
                    )
                if normalized is not None and normalized.get("expression_key") == target.get(
                    "expression_key"
                ):
                    self.changed = True
                    return ast.copy_location(ast.Name(id="_where_x", ctx=ast.Load()), node)
                return super().generic_visit(node)

        rewriter = _Rewriter(self)
        new_body = rewriter.visit(tree.body)
        if not rewriter.changed:
            return None
        return ast.unparse(new_body)

    def _try_expression_index_where(self, expr_result: blosc2.LazyExpr, catalog: dict) -> np.ndarray | None:
        """Attempt to resolve *expr_result* via a direct table expression index."""
        from blosc2.indexing import evaluate_bucket_query, evaluate_segment_query, plan_query

        expression = expr_result.expression
        operands = dict(expr_result.operands)
        for lookup_key, descriptor in catalog.items():
            target = descriptor.get("target") or {}
            if target.get("source") != "expression" or descriptor.get("stale", False):
                continue
            rewritten = self._rewrite_expression_query_for_index(expression, operands, target)
            if rewritten is None:
                continue
            expr_arr = self._index_target_array(lookup_key, descriptor)
            where_dict = {"_where_x": expr_arr}
            merged_operands = {"_where_x": expr_arr}
            plan = plan_query(rewritten, merged_operands, where_dict)
            if not plan.usable:
                continue
            if plan.exact_positions is not None:
                return np.asarray(plan.exact_positions, dtype=np.int64)
            if plan.bucket_masks is not None:
                _, positions = evaluate_bucket_query(
                    rewritten, merged_operands, {}, where_dict, plan, return_positions=True
                )
                return np.asarray(positions, dtype=np.int64)
            if plan.candidate_units is not None and plan.segment_len is not None:
                _, positions = evaluate_segment_query(
                    rewritten, merged_operands, {}, where_dict, plan, return_positions=True
                )
                return np.asarray(positions, dtype=np.int64)
        return None

    @staticmethod
    def _evaluate_expression_at(expr_result, candidates, *, prefetched: dict | None = None):
        """Evaluate *expr_result* on the operand rows at *candidates*.

        Returns a boolean ``numpy.ndarray`` the same length as *candidates*,
        or ``None`` if evaluation fails.

        Parameters
        ----------
        prefetched:
            Optional dict mapping operand variable names to already-gathered
            NumPy arrays.  When provided, those operands are reused instead of
            re-read from storage.
        """
        try:
            operands = {}
            for var_name, arr in expr_result.operands.items():
                if prefetched is not None and var_name in prefetched:
                    sliced = prefetched[var_name]
                else:
                    sliced = arr[candidates]
                    if hasattr(sliced, "__array__"):
                        sliced = np.asarray(sliced)
                operands[var_name] = sliced
            return blosc2.evaluate(expr_result.expression, operands)
        except Exception:
            return None

    @staticmethod
    def _find_indexed_columns(root_cols, catalog, operands):
        """Return live indexed columns referenced by *operands* in expression order.

        Avoid iterating over ``root_cols.items()`` here: for lazy persistent tables
        that would open every column just to find the indexed operands.
        """
        indexed = []
        seen = set()
        indexed_arrays = {}
        # Any column referenced by *operands* was already opened to build the
        # predicate expression, so it is present in the column cache.  A still
        # un-materialized indexed column therefore cannot be an operand: skip it
        # instead of opening it.  Without this, a predicate on a few columns
        # would eagerly open *every* indexed column (e.g. all 11 SUMMARY indexes
        # on a wide persistent table) — a pure cold-start cost.  Falling back to
        # the old behaviour for non-dict mappings keeps correctness; at worst a
        # genuinely-unloaded operand would miss the index and full-scan.
        cache_is_dict = isinstance(root_cols, dict)
        for col_name, descriptor in catalog.items():
            if col_name not in root_cols:
                continue
            if cache_is_dict and not dict.__contains__(root_cols, col_name):
                continue
            indexed_arrays[col_name] = (root_cols[col_name], descriptor)

        for operand in operands.values():
            if not isinstance(operand, blosc2.NDArray):
                continue
            for col_name, (col_arr, descriptor) in indexed_arrays.items():
                if col_name in seen or col_arr is not operand:
                    continue
                _CTableIndexingMixin._validate_index_descriptor(col_name, descriptor)
                if descriptor.get("stale", False):
                    continue
                indexed.append((col_name, col_arr, descriptor))
                seen.add(col_name)
        return indexed

    @staticmethod
    def _project_units_to_block_bitmap(primary_col_arr, candidate_units, segment_len) -> np.ndarray | None:
        """Project surviving flat segments onto a per-block candidate bitmap in
        *miniexpr* block coordinates (1-D only).

        Returns a ``uint8`` array with one byte per global miniexpr block
        (``nchunk * blocks_per_chunk + nblock``): 1 = the block overlaps a
        surviving segment and must be evaluated, 0 = it can be skipped.  A block
        is marked conservatively (any overlap with a survivor), which is always
        safe because the prefilter re-evaluates the exact predicate per row.
        """
        chunks = getattr(primary_col_arr, "chunks", None)
        blocks = getattr(primary_col_arr, "blocks", None)
        shape = getattr(primary_col_arr, "shape", None)
        if not chunks or not blocks or not shape or len(shape) != 1:
            return None
        chunk_len = int(chunks[0])
        block_len = int(blocks[0])
        nrows = int(shape[0])
        seg_len = int(segment_len)
        if chunk_len <= 0 or block_len <= 0 or nrows <= 0 or seg_len <= 0:
            return None
        cu = np.asarray(candidate_units, dtype=bool)
        nseg = len(cu)
        blocks_per_chunk = -(-chunk_len // block_len)  # ceil
        nchunks = -(-nrows // chunk_len)
        nblocks = nchunks * blocks_per_chunk
        g = np.arange(nblocks, dtype=np.int64)
        nchunk = g // blocks_per_chunk
        nblock = g % blocks_per_chunk
        lo = nchunk * chunk_len + nblock * block_len
        chunk_end = np.minimum((nchunk + 1) * chunk_len, nrows)
        hi = np.minimum(lo + block_len, chunk_end)
        has_rows = lo < chunk_end  # padding-only blocks have no real rows
        # Segments overlapping [lo, hi); test via a cumulative-sum range query.
        s_lo = np.clip(lo // seg_len, 0, nseg)
        s_hi = np.clip(np.where(has_rows, (hi - 1) // seg_len, lo // seg_len) + 1, 0, nseg)
        csum = np.concatenate(([0], np.cumsum(cu.astype(np.int64))))
        any_survivor = (csum[s_hi] - csum[s_lo]) > 0
        return (has_rows & any_survivor).astype(np.uint8)

    @staticmethod
    def _summary_candidate_block_bitmap(primary_col_arr, plan) -> np.ndarray | None:
        """Single-column candidate-block bitmap from a SUMMARY plan."""
        return _CTableIndexingMixin._project_units_to_block_bitmap(
            primary_col_arr, plan.candidate_units, plan.segment_len
        )

    @staticmethod
    def _flatten_conjunction(expression: str):
        """Return the list of ``ast.Compare`` leaves of a pure conjunction
        (``&`` / ``and`` of comparisons), or ``None`` if the expression is not a
        pure conjunction of comparisons (e.g. it contains ``|``, ``~``, or a
        bare boolean operand)."""
        import ast

        try:
            root = ast.parse(expression, mode="eval").body
        except (SyntaxError, ValueError):
            return None
        out = []
        stack = [root]
        while stack:
            node = stack.pop()
            if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
                stack.extend(node.values)
            elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitAnd):
                stack.extend((node.left, node.right))
            elif isinstance(node, ast.Compare):
                out.append(node)
            else:
                return None  # not a pure conjunction of comparisons
        return out

    def _combined_summary_block_bitmap(self, expr_result, indexed_columns, primary_col_arr):
        """Combine per-column SUMMARY block masks for a conjunction.

        For ``(a > x) & (b > y) & (c < z)`` where several of ``a``/``b``/``c``
        carry SUMMARY indexes, each column's summaries yield a candidate-block
        mask (in the shared miniexpr-block grid); ANDing them prunes a block
        unless *every* indexed predicate could match in it.  Returns the combined
        ``uint8`` bitmap, or ``None`` when the expression is not a pure
        conjunction or no indexed column contributes a usable range predicate.
        """
        from blosc2.indexing import (
            _candidate_units_from_summary_handle,
            _clear_cached_data,
            _finest_level,
            _open_level_summary_handle,
            _target_from_compare,
        )

        compares = self._flatten_conjunction(expr_result.expression)
        if compares is None:
            return None
        operands = expr_result.operands
        desc_by_id = {id(arr): desc for _name, arr, desc in indexed_columns}
        per_col: dict[int, np.ndarray] = {}
        for node in compares:
            target = _target_from_compare(node, operands)
            if target is None:
                continue
            base, _target_info, op, value = target
            desc = desc_by_id.get(id(base))
            if desc is None or desc.get("kind") != "summary":
                continue  # non-indexed predicate: contributes no pruning (sound for AND)
            level = _finest_level(desc)
            try:
                # In compact stores every column shares one urlpath, so the
                # sidecar-handle cache key collides across columns; clear it so
                # each column opens *its own* summary sidecar, not a sibling's.
                _clear_cached_data(base, desc["token"])
                summaries = _open_level_summary_handle(base, desc, level)
                cu = _candidate_units_from_summary_handle(summaries, op, value, np.dtype(desc["dtype"]))
            except (RuntimeError, ValueError, TypeError, KeyError):
                continue
            seg_len = desc["levels"][level]["segment_len"]
            bm = self._project_units_to_block_bitmap(primary_col_arr, cu, seg_len)
            if bm is None:
                continue
            # Multiple predicates on the same column AND together (e.g. a>1 & a<5).
            per_col[id(base)] = bm if id(base) not in per_col else (per_col[id(base)] & bm)
        if not per_col:
            return None
        combined = None
        for bm in per_col.values():
            combined = bm if combined is None else (combined & bm)
        return combined

    def _block_pruned_positions(self, expr_result, bitmap, primary_col_arr):
        """Evaluate the predicate through miniexpr with *bitmap* and return the
        matching physical positions.

        Chunks with no candidate block are skipped inside fast_eval
        (``_prune_chunks``), and positions are extracted by **streaming** the
        candidate chunks one at a time (``decompress_chunk`` → ``flatnonzero`` →
        offset).  This never materializes the full boolean mask as one big numpy
        array — only one chunk (~chunk_len bytes) is uncompressed at a time —
        and is faster than a bulk ``mask[:]`` + ``flatnonzero`` even when every
        chunk is a candidate.  Returns the positions array, or ``None`` if the
        bitmap path is unavailable.
        """
        chunks = getattr(primary_col_arr, "chunks", None)
        blocks = getattr(primary_col_arr, "blocks", None)
        shape = getattr(primary_col_arr, "shape", None)
        chunk_mask = None
        if chunks and blocks and shape and len(shape) == 1:
            chunk_len = int(chunks[0])
            block_len = int(blocks[0])
            bpc = -(-chunk_len // block_len)
            if block_len > 0 and len(bitmap) % bpc == 0:
                chunk_mask = bitmap.reshape(len(bitmap) // bpc, bpc).any(axis=1)
        try:
            if chunk_mask is None:
                # Geometry unknown (e.g. non-1-D): fall back to a bulk mask read.
                mask = expr_result.compute(_candidate_blocks=bitmap)
                mask_arr = mask[:] if isinstance(mask, blosc2.NDArray) else np.asarray(mask)
                if mask_arr.dtype != np.bool_:
                    return None
                return np.flatnonzero(mask_arr).astype(np.int64)

            nrows = int(shape[0])
            mask = expr_result.compute(_candidate_blocks=bitmap, _prune_chunks=True)
            if not isinstance(mask, blosc2.NDArray) or mask.dtype != np.bool_:
                return None
            if mask.schunk.nchunks != len(chunk_mask):
                return None  # unexpected geometry; bail rather than read garbage
            parts = []
            for chunk_idx in np.flatnonzero(chunk_mask):
                chunk_idx = int(chunk_idx)
                lo = chunk_idx * chunk_len
                vlen = min(chunk_len, nrows - lo)  # trim block-alignment padding
                sub = np.frombuffer(mask.schunk.decompress_chunk(chunk_idx), dtype=np.bool_)[:vlen]
                rel = np.flatnonzero(sub)
                if len(rel):
                    parts.append(rel.astype(np.int64) + lo)
            return np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)
        except Exception:
            return None

    def _try_index_where(self, expr_result: blosc2.LazyExpr) -> np.ndarray | None:  # noqa: C901
        """Attempt to resolve *expr_result* via a column index.

        Returns a 1-D int64 array of physical row positions that satisfy the
        predicate, or ``None`` if no usable index was found (caller falls back
        to a full scan).
        """
        from blosc2.indexing import (
            _IN_MEMORY_INDEXES,
            _PERSISTENT_INDEXES,
            _array_key,
            _default_index_store,
            _is_persistent_array,
            evaluate_bucket_query,
            evaluate_segment_query,
            plan_query,
        )

        root = self._root_table
        catalog = root._get_index_catalog()
        if not catalog:
            return None

        positions = self._try_expression_index_where(expr_result, catalog)
        if positions is not None:
            return positions

        expression = expr_result.expression
        operands = dict(expr_result.operands)

        indexed_columns = self._find_indexed_columns(root._cols, catalog, operands)
        if not indexed_columns:
            return None

        primary_col_name, primary_col_arr, _ = indexed_columns[0]
        nullable_indexed = [
            name
            for name, _arr, _descriptor in indexed_columns
            if getattr(root._schema.columns_by_name[name].spec, "null_value", None) is not None
        ]
        # A NaN null sentinel can never satisfy a comparison (every comparison
        # with NaN is False), so the predicate itself already drops those rows;
        # only non-NaN sentinels (which *can* match a predicate) need explicit
        # position-based exclusion.  This lets the fast mask-direct path apply to
        # NaN-nullable columns instead of falling back to the positions path.
        nullable_needs_exclude = []
        for name in nullable_indexed:
            nv = getattr(root._schema.columns_by_name[name].spec, "null_value", None)
            if not (isinstance(nv, float) and np.isnan(nv)):
                nullable_needs_exclude.append(name)

        # Global null post-filtering is not correct for OR expressions.
        if nullable_indexed and ("|" in expr_result.expression or " or " in expr_result.expression):
            return None

        # Inject every usable table-owned descriptor so plan_query can combine them.
        # In .b2z read mode all columns share the same urlpath, so _array_key()
        # returns the same key for every column — causing _SIDECAR_HANDLE_CACHE
        # collisions across queries.  Clear stale handles before each injection so
        # the upcoming query always loads the correct sidecar for this column.
        from blosc2.indexing import _clear_cached_data, _register_descriptor_owner

        # Build column-name → array mapping so the planner can resolve the
        # correct per-column token instead of the generic "__self__".
        array_to_col = {}
        for _col_name, col_arr, descriptor in indexed_columns:
            array_to_col[id(col_arr)] = _col_name
            arr_key = _array_key(col_arr)
            if _is_persistent_array(col_arr):
                store = _PERSISTENT_INDEXES.get(arr_key) or _default_index_store()
                if store["indexes"].get(_col_name) is not descriptor:
                    _clear_cached_data(col_arr, _col_name)
                store["indexes"][_col_name] = descriptor
                _PERSISTENT_INDEXES[arr_key] = store
            else:
                store = _IN_MEMORY_INDEXES.get(id(col_arr)) or _default_index_store()
                store["indexes"][_col_name] = descriptor
                _IN_MEMORY_INDEXES[id(col_arr)] = store
            # Record the owning array so a sibling column sharing this urlpath
            # (and, after column alignment, the same shape/chunks) cannot match
            # this descriptor in _descriptor_for_target.
            _register_descriptor_owner(col_arr, _col_name)

        where_dict = {"_where_x": primary_col_arr}
        merged_operands = {**operands, "_where_x": primary_col_arr}

        plan = plan_query(expression, merged_operands, where_dict, array_to_col=array_to_col)
        if not plan.usable:
            return None

        def _exclude_null_positions(positions):
            positions = np.asarray(positions, dtype=np.int64)
            for name in nullable_indexed:
                col = root._schema.columns_by_name[name]
                raw = root._cols[name][positions]
                nv = getattr(col.spec, "null_value", None)
                if isinstance(nv, float) and np.isnan(nv):
                    keep = ~np.isnan(raw)
                else:
                    keep = raw != nv
                positions = positions[keep]
            return positions

        if plan.exact_positions is not None:
            return _exclude_null_positions(plan.exact_positions)

        if plan.partial_exact_positions is not None:
            # Cross-column refinement: the FULL index on one column gave us
            # exact positions, but the expression has additional predicates on
            # other columns.  Refinement reads every operand column at those
            # candidate positions using sparse/fancy indexing.  For compressed
            # columns this can touch many chunks and be slower than the regular
            # sequential miniexpr scan, which is very fast for simple predicates.
            # Use a cost model to compare refinement vs full scan.
            candidates = np.asarray(plan.partial_exact_positions, dtype=np.int64)
            n_candidates = len(candidates)
            n_operands = len(expr_result.operands)
            target_len = len(root._valid_rows)

            estimated_refine_ms = (
                (n_candidates / 1000.0) * self._GATHER_COST_MS_PER_1K_ITEMS_PER_OP * n_operands
            )
            estimated_scan_ms = (target_len / 1_000_000.0) * self._SCAN_COST_MS_PER_1M_ROWS
            if estimated_refine_ms > estimated_scan_ms:
                return None

            # Read the primary column once and reuse for both null filtering
            # and refinement, avoiding a second sparse gather later.
            primary_op_name = next(
                (vn for vn, va in expr_result.operands.items() if va is primary_col_arr), None
            )
            prefetched = None
            if nullable_indexed and primary_op_name is not None:
                raw = primary_col_arr[candidates]
                raw = np.asarray(raw) if hasattr(raw, "__array__") else raw
                pos = candidates
                for name in nullable_indexed:
                    if name == primary_col_name:
                        nv = getattr(root._schema.columns_by_name[name].spec, "null_value", None)
                        if isinstance(nv, float) and np.isnan(nv):
                            keep = ~np.isnan(raw)
                        else:
                            keep = raw != nv
                        pos = pos[keep]
                        raw = raw[keep]  # already filtered for refinement reuse
                    else:
                        col = root._schema.columns_by_name[name]
                        vals = root._cols[name][pos]
                        nv = getattr(col.spec, "null_value", None)
                        if isinstance(nv, float) and np.isnan(nv):
                            keep = ~np.isnan(vals)
                        else:
                            keep = vals != nv
                        pos = pos[keep]
                candidates = pos
                prefetched = {primary_op_name: raw}
            else:
                candidates = _exclude_null_positions(candidates)

            restricted = self._evaluate_expression_at(expr_result, candidates, prefetched=prefetched)
            if restricted is not None and restricted.dtype == np.bool_:
                refined = candidates[np.asarray(restricted, dtype=bool)]
                return _exclude_null_positions(refined)
            # Fall through to full scan if refinement fails

        if plan.bucket_masks is not None:
            # When bucket pruning covers all units (100 % of chunks are
            # candidates), the per‑chunk evaluation overhead outweighs the
            # benefit over a plain scan.  Fall back to the scan path.
            if plan.total_units > 0 and plan.selected_units >= plan.total_units:
                return None
            _, positions = evaluate_bucket_query(
                expression, merged_operands, {}, where_dict, plan, return_positions=True
            )
            return _exclude_null_positions(positions)

        if plan.candidate_units is not None and plan.segment_len is not None:
            # Cheap pre-gate on the primary column's pruning (the plan is already
            # computed).  When the primary alone cannot prune enough, a plain scan
            # is cheaper, so fall back before doing any extra (per-column) work.
            n_ops = max(1, len(expr_result.operands))
            segment_len = int(plan.segment_len)
            selected = int(plan.selected_units)
            total_rows = int(plan.total_units) * segment_len
            selected_rows = selected * segment_len
            scan_ms = (total_rows / 1e6) * self._SEGMENT_SCAN_MS_PER_1M_ROWS_PER_OP * n_ops
            index_ms = (
                selected * n_ops * self._SEGMENT_OVERHEAD_MS_PER_OP
                + (selected_rows / 1e6) * self._SEGMENT_SCAN_MS_PER_1M_ROWS_PER_OP * n_ops
            )
            if selected > 0 and index_ms >= scan_ms:
                return None
            # Preferred path: evaluate the predicate through miniexpr with a
            # candidate-block bitmap, so non-surviving blocks are skipped inside
            # the fast_eval/miniexpr engine (no per-segment Python loop).  For a
            # conjunction, every column with a SUMMARY index contributes a block
            # mask and they are ANDed, pruning more than the primary alone.  The
            # masked result equals a full scan's (pruned blocks have no matches),
            # so positions are exact.  Falls back to the per-segment query if the
            # bitmap path is unavailable.
            bitmap = self._combined_summary_block_bitmap(expr_result, indexed_columns, primary_col_arr)
            if bitmap is None:
                bitmap = self._summary_candidate_block_bitmap(primary_col_arr, plan)
            if bitmap is not None:
                # Mask-direct: SUMMARY evaluates the predicate through miniexpr
                # with the candidate-block bitmap, so the result is already a
                # compressed boolean mask identical to a full scan's.  Return it
                # as-is (no positions round-trip) when no null post-filtering is
                # needed; the caller views it directly and extracts positions
                # lazily.  Only columns whose null sentinel can satisfy the
                # predicate (non-NaN) require the positions fall-back; NaN
                # sentinels are already excluded by the predicate itself.
                if not nullable_needs_exclude:
                    try:
                        mask = expr_result.compute(_candidate_blocks=bitmap)
                    except Exception:
                        mask = None
                    if isinstance(mask, blosc2.NDArray) and mask.dtype == np.bool_:
                        return mask
                positions = self._block_pruned_positions(expr_result, bitmap, primary_col_arr)
                if positions is not None:
                    return _exclude_null_positions(positions)
            _, positions = evaluate_segment_query(
                expression, merged_operands, {}, where_dict, plan, return_positions=True
            )
            return _exclude_null_positions(positions)

        return None
