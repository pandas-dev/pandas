#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

"""Group-by support for :class:`blosc2.CTable`.

This module contains the Phase-1, NumPy-based implementation.  It is deliberately
chunked and columnar: only grouping columns, aggregation columns, and the
live-row mask are read from the source table.
"""

from __future__ import annotations

import copy
import dataclasses
import math
from collections.abc import Callable, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from blosc2.dsl_kernel import DSLKernel
from blosc2.schema import DictionarySpec, NDArraySpec, SchemaSpec, float64, int64
from blosc2.schema import bool as b2_bool
from blosc2.schema import field as b2_field
from blosc2.schema_compiler import _validate_column_name

if TYPE_CHECKING:  # pragma: no cover
    from blosc2.ctable import CTable


AggName = Literal["size", "count", "sum", "mean", "min", "max", "argmin", "argmax", "udf"]

_NAN_KEY = ("__blosc2_groupby_nan__",)


@dataclasses.dataclass
class _AggSpec:
    input_col: str | None
    op: AggName
    output_col: str
    udf: Callable | None = None
    explicit_dtype: SchemaSpec | None = None


@dataclasses.dataclass
class _AggState:
    op: AggName
    value: Any = None
    count: int = 0


@dataclasses.dataclass
class _Utf8KeyChunk:
    """A utf8 key-column chunk, factorized to chunk-local integer codes.

    ``codes[i]`` indexes ``uniques`` (a ``StringDType`` array sorted
    ascending), so null detection, live-row masking, and per-chunk
    ``np.unique`` all run on int64 codes; only the (few) distinct strings are
    ever decoded.  Produced by :meth:`CTableGroupBy._read_key_chunk` via
    ``Utf8Array.factorizer``.
    """

    codes: np.ndarray
    uniques: np.ndarray

    def __len__(self) -> int:
        return len(self.codes)

    def take(self, mask: np.ndarray) -> _Utf8KeyChunk:
        return _Utf8KeyChunk(self.codes[mask], self.uniques)

    def code_of(self, value: str) -> int:
        """Code of *value* in this chunk, or -1 when absent (uniques are sorted)."""
        i = int(np.searchsorted(self.uniques, value))
        if i < len(self.uniques) and self.uniques[i] == value:
            return i
        return -1


def _is_column_like(value: Any) -> bool:
    return isinstance(getattr(value, "_col_name", None), str)


def _column_name(value: Any) -> str:
    name = getattr(value, "_col_name", value)
    if not isinstance(name, str):
        raise TypeError(f"Expected a column name or Column object, got {type(value)!r}")
    return name


_OP_ALIAS_CACHE: dict | None = None


def _op_alias_map() -> dict:
    """Map blosc2 reduction *function objects* to group-by op names, by identity.

    Built lazily (and cached) to avoid importing the top-level package at module
    load time.  Only the functions that have a matching group-by op are included;
    matching is by object identity, so a user function merely named ``sum`` does
    not collide with :func:`blosc2.sum`.
    """
    global _OP_ALIAS_CACHE
    if _OP_ALIAS_CACHE is None:
        import blosc2

        mapping = {}
        for op in ("sum", "mean", "min", "max", "argmin", "argmax"):
            fn = getattr(blosc2, op, None)
            if fn is not None:
                mapping[fn] = op
        _OP_ALIAS_CACHE = mapping
    return _OP_ALIAS_CACHE


class CTableGroupBy:
    """Deferred group-by operation returned by :meth:`CTable.group_by`.

    The object stores the source table, grouping keys, and execution options.
    It is not a :class:`CTable` view and does not materialize grouped data until
    a terminal method such as :meth:`size`, :meth:`count`, or :meth:`agg` is
    called.
    """

    def __init__(
        self,
        table: CTable,
        keys: str | Sequence[str],
        *,
        sort: bool | None = None,
        dropna: bool = True,
        engine: str = "auto",
        chunk_size: int | None = None,
    ) -> None:
        if isinstance(keys, str) or _is_column_like(keys):
            keys = [_column_name(keys)]
        else:
            keys = [_column_name(k) for k in keys]
        if not keys:
            raise ValueError("group_by() requires at least one key column")

        self.table = table
        self.keys = [table._logical_to_physical_name(k) for k in keys]
        # Tri-state ordering request resolved per execution path by
        # _resolve_sort(): True forces a key sort, False emits first-appearance
        # order, None (the default) sorts only when the path can do so cheaply.
        self.sort = sort
        self.dropna = bool(dropna)
        self.engine = engine
        self.chunk_size = chunk_size
        # Per-key incremental Utf8Factorizer instances, shared across the
        # chunk loop so the string vocabulary is built once (see
        # _read_key_chunk).
        self._utf8_factorizers: dict[str, Any] = {}

        for name in self.keys:
            if name in table._computed_cols:
                raise NotImplementedError("group_by() over computed columns is not supported yet")
            if name not in table._cols:
                raise KeyError(f"No column named {name!r}. Available: {table.col_names}")
            table._ensure_generated_column_not_stale(name)
            col_info = table._schema.columns_by_name[name]
            if isinstance(col_info.spec, NDArraySpec):
                raise TypeError(
                    f"Cannot group by ndarray column {name!r} with per-row shape {col_info.spec.item_shape}. "
                    "Materialize a scalar generated column first, e.g. embedding_norm or embedding_max."
                )
            if table._is_list_column(col_info) or (
                table._is_varlen_scalar_column(col_info) and not table._is_utf8_column(col_info)
            ):
                raise TypeError(f"Cannot group by variable-length/list column {name!r} in Phase 1")

    def size(self, *, urlpath: str | None = None):
        """Return row counts per group as a new :class:`CTable`.

        This is equivalent to SQL ``COUNT(*)``: it counts rows in each group and
        is independent of null values in non-key columns.  If *urlpath* is
        provided, the result is written as a persistent CTable at that path.
        """
        return self._execute([_AggSpec(None, "size", "size")], urlpath=urlpath)

    def count(self, column: str, *, urlpath: str | None = None):
        """Return non-null value counts for *column* per group.

        This is equivalent to SQL ``COUNT(column)`` and to
        ``group_by(...).agg({column: "count"})``.
        """
        col = self.table._logical_to_physical_name(_column_name(column))
        return self._execute([_AggSpec(col, "count", f"{col}_count")], urlpath=urlpath)

    def sum(self, column: str, *, urlpath: str | None = None):
        """Return sums of *column* per group.

        This is equivalent to ``group_by(...).agg({column: "sum"})``.
        """
        return self.agg({_column_name(column): "sum"}, urlpath=urlpath)

    def mean(self, column: str, *, urlpath: str | None = None):
        """Return means of *column* per group.

        This is equivalent to ``group_by(...).agg({column: "mean"})``.
        """
        return self.agg({_column_name(column): "mean"}, urlpath=urlpath)

    def min(self, column: str, *, urlpath: str | None = None):
        """Return minimum values of *column* per group.

        This is equivalent to ``group_by(...).agg({column: "min"})``.
        """
        return self.agg({_column_name(column): "min"}, urlpath=urlpath)

    def max(self, column: str, *, urlpath: str | None = None):
        """Return maximum values of *column* per group.

        This is equivalent to ``group_by(...).agg({column: "max"})``.
        """
        return self.agg({_column_name(column): "max"}, urlpath=urlpath)

    def argmin(self, column: str, *, urlpath: str | None = None):
        """Return logical row positions of minimum non-null *column* values per group.

        Ties keep the first row in the grouped input table or view.  Groups with
        no non-null values for *column* receive ``-1``.
        """
        return self.agg({_column_name(column): "argmin"}, urlpath=urlpath)

    def argmax(self, column: str, *, urlpath: str | None = None):
        """Return logical row positions of maximum non-null *column* values per group.

        Ties keep the first row in the grouped input table or view.  Groups with
        no non-null values for *column* receive ``-1``.
        """
        return self.agg({_column_name(column): "argmax"}, urlpath=urlpath)

    def agg(
        self,
        aggregations: Mapping[str, str | Sequence[str]]
        | Sequence[tuple[Any, str | Sequence[str]]]
        | None = None,
        *,
        urlpath: str | None = None,
        **named: tuple[str, str],
    ):
        """Aggregate value columns per group.

        Three ways to specify aggregations and name the output columns, which may
        be combined:

        * **Auto-named mapping** -- pass *aggregations* as a mapping from input
          column name to an aggregation name (or list of names).  Each result
          column is named ``"<column>_<op>"`` (e.g. ``price_sum``).  Compact and
          collision-safe when a column is aggregated several ways.
        * **Auto-named list of pairs** -- pass *aggregations* as a list of
          ``(column, op-or-ops)`` pairs.  Same ``"<column>_<op>"`` naming, but
          unlike the mapping form it accepts :class:`~blosc2.ctable.Column`
          objects (e.g. ``t.price``), which cannot be dict keys.
        * **Explicitly named** -- pass ``output_name=(column, op)`` keyword
          arguments (pandas-style named aggregation), giving the result column
          exactly the name you want.  This is also the only form that accepts a
          **custom UDF aggregation**: ``output_name=(column, callable)``, or
          ``output_name=(column, callable, dtype)`` to give the output column's
          schema spec explicitly instead of inferring it from the callable's
          results.

        Parameters
        ----------
        aggregations:
            Either a mapping ``{column: op-or-ops}`` or a list of
            ``(column, op-or-ops)`` pairs.  Columns may be name strings or
            :class:`~blosc2.ctable.Column` objects (list/named forms only).
            Supported operations are ``"count"``, ``"sum"``, ``"mean"``,
            ``"min"``, ``"max"``, ``"argmin"``, ``"argmax"`` and the special
            row-count spelling ``"*": "size"``.  An op may also be given as the
            corresponding blosc2 reduction *function* (``blosc2.sum``, ``mean``,
            ``min``, ``max``, ``argmin``, ``argmax``), matched by identity; this
            is a naming shorthand, not a UDF mechanism.  Result columns are
            named ``"<column>_<op>"``.  Custom callables are only accepted via
            the named-aggregation form (see below).
        urlpath:
            If given, write the result as a persistent CTable at that path.
        **named:
            Named aggregations as ``output_name=(column, op)`` pairs, or
            ``output_name=(column, callable[, dtype])`` for a custom UDF
            aggregation.  Use ``("*", "size")`` for a row count.

            A UDF callable receives a 1-D NumPy array of the group's live,
            non-null values (nulls are pre-filtered, same as the built-in
            aggregations) and returns a scalar.  It is called once per group
            with a plain Python loop -- there is no acceleration yet, this is
            the "slow but correct" baseline and semantics oracle for any
            future JIT path.  The output dtype is inferred from the results
            across *every* group (raising a clear error if they disagree) unless
            *dtype* is given explicitly as a blosc2 schema spec.

        Examples
        --------
        >>> g = t.group_by("city")  # doctest: +SKIP
        >>> # Auto-named mapping: columns become sales_sum / sales_mean.
        >>> g.agg({"sales": ["sum", "mean"]})  # doctest: +SKIP
        >>> # Auto-named list of pairs: same names, but accepts Column objects
        >>> # and blosc2 reduction functions as ops.
        >>> g.agg([(t.sales, [blosc2.sum, "mean"])])  # doctest: +SKIP
        >>> # Explicitly named: columns become revenue / avg_sale.
        >>> g.agg(revenue=("sales", "sum"), avg_sale=("sales", "mean"))  # doctest: +SKIP
        >>> # Forms combine, e.g. a list of pairs plus a named row count.
        >>> g.agg([(t.sales, "sum")], n=("*", "size"))  # doctest: +SKIP
        >>> # Custom UDF aggregation (named form only).
        >>> g.agg(sales_range=("sales", lambda a: a.max() - a.min()))  # doctest: +SKIP
        """
        specs = self._normalize_aggs(aggregations, named)
        return self._execute(specs, urlpath=urlpath)

    def _resolve_op(self, op):
        """Resolve an aggregation *op* to its name string.

        Accepts a name string, or one of blosc2's own reduction *functions*
        (:func:`blosc2.sum`, ``mean``, ``min``, ``max``, ``argmin``, ``argmax``)
        matched **by identity** -- so a user function that merely shares a name
        (e.g. a UDF called ``sum``) is *not* silently accepted here.  Custom UDF
        aggregations are only accepted via the named form (``output_name=(col,
        callable)``, see :meth:`agg`); this method only ever sees a callable
        when it fell through that path (auto-named mapping/list forms, which
        cannot derive an output column name for an arbitrary callable).
        """
        if isinstance(op, str):
            return op
        if callable(op):
            alias = _op_alias_map().get(op)
            if alias is not None:
                return alias
            raise ValueError(
                f"Unsupported aggregation function {getattr(op, '__name__', op)!r}.  Pass a "
                f"string op name (e.g. 'sum') or a blosc2 reduction function "
                f"(blosc2.sum/mean/min/max/argmin/argmax).  Custom UDF aggregations are "
                f"supported only via the named form, e.g. "
                f"g.agg(my_range=(col, {getattr(op, '__name__', 'my_func')}))."
            )
        raise ValueError(f"Aggregation op must be a string or a blosc2 reduction function, got {op!r}")

    def _build_agg_spec(
        self, col_name, op, output_col: str | None = None, *, dtype: SchemaSpec | None = None
    ) -> _AggSpec:
        """Validate a single (column, op) pair and build its :class:`_AggSpec`.

        ``output_col`` overrides the default ``"<column>_<op>"`` name; ``"*"`` as
        *col_name* (only with ``op="size"``) yields a row count.  *col_name* may
        be a column name string or a :class:`~blosc2.ctable.Column` object; *op*
        may be an op name string, a blosc2 reduction function (see
        :meth:`_resolve_op`), or an arbitrary callable UDF aggregation (see
        :meth:`agg`), in which case *dtype* optionally gives the output
        column's schema spec (e.g. ``blosc2.float64()``) instead of inferring
        it from the UDF's first result.
        """
        if output_col is not None and callable(op) and op not in _op_alias_map():
            # Only the named form (output_col given) may carry a custom UDF;
            # the auto-named mapping/list forms cannot derive a column name
            # for an arbitrary callable, so they fall through to _resolve_op()
            # below and get its "unsupported aggregation function" error.
            physical = self.table._logical_to_physical_name(_column_name(col_name))
            self._validate_value_column(physical)
            return _AggSpec(physical, "udf", output_col, udf=op, explicit_dtype=dtype)
        op = self._resolve_op(op)
        # Guard the "*" check with isinstance: a Column object overloads __eq__
        # to build an expression, so a bare ``col_name == "*"`` would not return
        # a bool.  Only the literal string "*" is the row-count sentinel.
        if isinstance(col_name, str) and col_name == "*":
            if op != "size":
                raise ValueError("Only the 'size' aggregation is supported for '*' input")
            return _AggSpec(None, "size", output_col or "size")
        physical = self.table._logical_to_physical_name(_column_name(col_name))
        self._validate_value_column(physical)
        if op not in {"count", "sum", "mean", "min", "max", "argmin", "argmax"}:
            raise ValueError(f"Unsupported aggregation {op!r}")
        self._validate_agg_for_column(physical, op)
        return _AggSpec(physical, op, output_col or f"{physical}_{op}")

    def _expand_ops(self, col_name, ops) -> list[_AggSpec]:
        """Expand a (column, op-or-ops) entry into one auto-named spec per op."""
        # A single op may be a string or a reduction function (both non-iterable
        # in the sense we want); only a real sequence of ops is expanded.
        op_list = [ops] if isinstance(ops, str) or callable(ops) else list(ops)
        if not op_list:
            raise ValueError(f"No aggregations specified for column {col_name!r}")
        return [self._build_agg_spec(col_name, op) for op in op_list]

    def _normalize_aggs(
        self,
        aggregations: Mapping[str, str | Sequence[str]] | Sequence[tuple[Any, str | Sequence[str]]] | None,
        named: Mapping[str, tuple[str, str]] | None = None,
    ) -> list[_AggSpec]:
        if not aggregations and not named:
            raise ValueError("agg() requires a mapping/list of pairs and/or named aggregations")
        specs: list[_AggSpec] = []
        if aggregations:
            if isinstance(aggregations, Mapping):
                entries = list(aggregations.items())
            elif isinstance(aggregations, Sequence) and not isinstance(aggregations, (str, bytes)):
                # List of (column, ops) pairs -- lets you use Column objects,
                # which cannot be dict keys (Column is unhashable).
                entries = []
                for pair in aggregations:
                    if not (isinstance(pair, (tuple, list)) and len(pair) == 2):
                        raise ValueError(
                            f"agg() positional list must contain (column, ops) pairs, got {pair!r}"
                        )
                    entries.append(pair)
            else:
                raise ValueError(
                    "agg() positional argument must be a mapping or a list of (column, ops) pairs"
                )
            for col_name, ops in entries:
                specs.extend(self._expand_ops(col_name, ops))
        for out_name, value in (named or {}).items():
            if not (isinstance(value, (tuple, list)) and len(value) in (2, 3)):
                raise ValueError(
                    f"Named aggregation {out_name!r} must be a (column, op) or "
                    f"(column, op, dtype) tuple, got {value!r}"
                )
            col_name, op, *rest = value
            dtype = rest[0] if rest else None
            if dtype is not None and not (callable(op) and op not in _op_alias_map()):
                raise ValueError(
                    f"Named aggregation {out_name!r}: an explicit dtype is only supported "
                    "for callable UDF aggregations"
                )
            if isinstance(op, (tuple, list, set)):
                raise ValueError(
                    f"Named aggregation {out_name!r} takes a single op, got {op!r}.  "
                    f"Each named output maps to one (column, op); use the mapping "
                    f"form agg({{column: [...]}}) for several ops, or give each its "
                    f"own name (e.g. {out_name}_sum=(col, 'sum'))."
                )
            specs.append(self._build_agg_spec(col_name, op, output_col=out_name, dtype=dtype))
        output_names = [s.output_col for s in specs]
        if len(output_names) != len(set(output_names)):
            raise ValueError("Aggregation output column names must be unique")
        return specs

    def _validate_agg_for_column(self, name: str, op: str) -> None:
        dtype = getattr(self.table._schema.columns_by_name[name].spec, "dtype", None)
        if op in {"sum", "mean"} and dtype is not None and dtype.kind not in "biuf":
            raise TypeError(f"Aggregation {op!r} is not supported for column {name!r} with dtype {dtype}")
        if op in {"min", "max", "argmin", "argmax"} and dtype is not None and dtype.kind == "c":
            raise TypeError(f"Aggregation {op!r} is not supported for complex column {name!r}")

    def _validate_value_column(self, name: str) -> None:
        if name in self.table._computed_cols:
            raise NotImplementedError("group_by() aggregations over computed columns are not supported yet")
        if name not in self.table._cols:
            raise KeyError(f"No column named {name!r}. Available: {self.table.col_names}")
        self.table._ensure_generated_column_not_stale(name)
        col_info = self.table._schema.columns_by_name[name]
        if self.table._is_list_column(col_info) or self.table._is_varlen_scalar_column(col_info):
            raise TypeError(f"Cannot aggregate variable-length/list column {name!r} in Phase 1")
        if isinstance(col_info.spec, NDArraySpec):
            raise TypeError(
                f"Cannot aggregate ndarray column {name!r} with per-row shape {col_info.spec.item_shape}. "
                "Materialize a scalar generated column first."
            )
        if self.table._is_dictionary_column(col_info):
            raise TypeError(f"Cannot aggregate dictionary column {name!r} in Phase 1")

    def _execute(self, specs: list[_AggSpec], *, urlpath: str | None = None):
        self._validate_output_names(specs)
        old_result_urlpath = getattr(self, "_result_urlpath", None)
        self._result_urlpath = urlpath
        try:
            return self._execute_with_result_target(specs)
        finally:
            self._result_urlpath = old_result_urlpath

    def _try_fast_paths(self, specs: list[_AggSpec], use_arg_positions: bool):
        """Dispatch to the available fast paths; return a result CTable or None.

        argmin/argmax can only use the dense single-key path (it tracks row
        positions); the Cython kernels do not, so they are skipped for them.
        UDF aggregations always fall through to the generic
        chunked path below, which is the only one that accumulates raw
        per-group values instead of a mergeable scalar state.
        """
        if any(s.op == "udf" for s in specs):
            return None
        if not use_arg_positions:
            for attempt in (
                self._try_execute_cython_dense_int_key,
                self._try_execute_cython_two_int_key_hash,
                self._try_execute_cython_i32_f64_sum,
                self._try_execute_cython_float_integral_key_f64_sum,
            ):
                fast = attempt(specs)
                if fast is not None:
                    return fast
        # Dense single-key path also covers integral-valued float keys with a
        # compact non-negative range (e.g. float32 second/id columns), so try it
        # before the generic float hash path, which is markedly slower.  Unlike
        # the Cython kernels above, it tracks row positions, so it also serves
        # argmin/argmax — keeping them off the slow generic hash+merge path.
        fast = self._try_execute_dense_single_int_key(specs)
        if fast is not None:
            return fast
        if not use_arg_positions:
            return self._try_execute_cython_float_hash(specs)
        return None

    def _execute_with_result_target(self, specs: list[_AggSpec]):
        use_arg_positions = any(s.op in {"argmin", "argmax"} for s in specs)
        fast = self._try_fast_paths(specs, use_arg_positions)
        if fast is not None:
            return fast

        acc: dict[Any, dict[str, _AggState]] = {}
        key_values: dict[Any, tuple[Any, ...]] = {}

        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()
        value_cols = sorted({s.input_col for s in specs if s.input_col is not None})

        logical_seen = 0
        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            logical_positions = logical_seen + np.cumsum(valid, dtype=np.int64) - 1
            logical_seen += int(np.count_nonzero(valid))
            if not np.any(valid):
                continue

            raw_keys = [self._read_key_chunk(name, start, stop) for name in self.keys]
            live_mask = valid.copy()
            if self.dropna:
                for name, values in zip(self.keys, raw_keys, strict=True):
                    live_mask &= ~self._null_mask(name, values, is_key=True)
            if not np.any(live_mask):
                continue

            keys_live = [
                values.take(live_mask)
                if isinstance(values, _Utf8KeyChunk)
                else np.asarray(values)[live_mask]
                for values in raw_keys
            ]
            n_live = len(keys_live[0])
            if n_live == 0:
                continue

            unique_keys, inverse = self._factorize_keys(keys_live)
            value_chunks = {
                name: np.asarray(self.table._cols[name][start:stop])[live_mask] for name in value_cols
            }

            partials = self._compute_partials(
                specs, unique_keys, inverse, value_chunks, logical_positions[live_mask]
            )
            display_keys = self._display_keys(unique_keys)
            normalized_keys = self._normalized_keys(display_keys)
            self._merge_partials(acc, key_values, normalized_keys, display_keys, partials, specs)

        rows = self._final_rows(acc, key_values, specs)
        return self._build_result(rows, specs)

    def _try_execute_cython_two_int_key_hash(self, specs: list[_AggSpec]):  # noqa: C901
        """Cython hash path for two integer/dictionary-code keys."""
        if len(self.keys) != 2:
            return None

        key_arrays = []
        key_is_dict = []
        key_nulls = []
        skip_key_nulls = []
        for key_name in self.keys:
            key_info = self.table._schema.columns_by_name[key_name]
            if self.table._is_dictionary_column(key_info):
                key_arrays.append(self.table._cols[key_name].codes)
                key_is_dict.append(True)
                key_nulls.append(int(key_info.spec.null_code))
                skip_key_nulls.append(self.dropna)
                continue
            key_dtype = getattr(key_info.spec, "dtype", None)
            if key_dtype is None or np.dtype(key_dtype).kind not in "biu":
                return None
            null_value = getattr(key_info.spec, "null_value", None)
            if null_value is not None and not self.dropna:
                return None
            key_arrays.append(self.table._cols[key_name])
            key_is_dict.append(False)
            key_nulls.append(0 if null_value is None else int(null_value))
            skip_key_nulls.append(self.dropna and null_value is not None)

        value_cols = {s.input_col for s in specs if s.input_col is not None}
        if len(value_cols) > 1:
            return None
        value_col = next(iter(value_cols), None)
        if value_col is not None and any(s.op in {"sum", "mean", "min", "max"} for s in specs):
            value_info = self.table._schema.columns_by_name[value_col]
            value_dtype = getattr(value_info.spec, "dtype", None)
            if value_dtype is None or np.dtype(value_dtype).kind != "f":
                return None
            null_value = getattr(value_info.spec, "null_value", None)
            if null_value is not None and not (isinstance(null_value, float) and math.isnan(null_value)):
                return None

        try:
            from blosc2 import groupby_ext
        except ImportError:
            return None
        kernel = getattr(groupby_ext, "groupby_hash_i64x2_f64", None)
        if kernel is None:
            return None

        acc: dict[Any, dict[str, _AggState]] = {}
        key_values: dict[Any, tuple[Any, ...]] = {}
        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()

        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if not np.any(valid):
                continue
            key_chunks = [np.asarray(arr[start:stop], dtype=np.int64) for arr in key_arrays]
            live = valid.copy()
            for key_chunk, skip_null, null_value in zip(key_chunks, skip_key_nulls, key_nulls, strict=True):
                if skip_null:
                    live &= key_chunk != null_value
            if not np.any(live):
                continue

            if value_col is None:
                values = np.empty(len(valid), dtype=np.float64)
                values_valid = np.zeros(len(valid), dtype=bool)
                has_values = False
            else:
                raw_values = np.asarray(self.table._cols[value_col][start:stop])
                values = np.ascontiguousarray(raw_values.astype(np.float64, copy=False))
                values_valid = np.ascontiguousarray(~self._null_mask(value_col, raw_values, is_key=False))
                has_values = True

            (
                out_k0,
                out_k1,
                row_counts,
                value_counts,
                sums,
                mins,
                maxs,
                has_value,
            ) = kernel(
                np.ascontiguousarray(key_chunks[0]),
                np.ascontiguousarray(key_chunks[1]),
                values,
                np.ascontiguousarray(live),
                values_valid,
                has_values,
            )

            for i, (code0, code1) in enumerate(zip(out_k0, out_k1, strict=True)):
                display = []
                norm_parts = []
                for key_pos, code in enumerate((int(code0), int(code1))):
                    if key_is_dict[key_pos]:
                        value = self.table._cols[self.keys[key_pos]].decode(code)
                    else:
                        value = code
                    display.append(value)
                    norm_parts.append(_normalize_key_part(value))
                norm_key = tuple(norm_parts)
                states = acc.setdefault(norm_key, {})
                key_values.setdefault(norm_key, tuple(display))
                for spec in specs:
                    state = states.setdefault(spec.output_col, _AggState(spec.op))
                    if spec.op == "size":
                        state.value = (0 if state.value is None else state.value) + int(row_counts[i])
                    elif spec.op == "count":
                        state.value = (0 if state.value is None else state.value) + int(value_counts[i])
                    elif spec.op in {"sum", "mean"}:
                        if has_value[i]:
                            state.value = (0.0 if state.value is None else state.value) + float(sums[i])
                            state.count += int(value_counts[i])
                    elif spec.op == "min":
                        if has_value[i]:
                            value = float(mins[i])
                            if state.count == 0 or value < state.value:
                                state.value = value
                            state.count += 1
                    elif spec.op == "max" and has_value[i]:
                        value = float(maxs[i])
                        if state.count == 0 or value > state.value:
                            state.value = value
                        state.count += 1

        rows = self._final_rows(acc, key_values, specs)
        return self._build_result(rows, specs)

    def _try_execute_cython_dense_int_key(self, specs: list[_AggSpec]):  # noqa: C901
        """Cython fast path for one compact integer/dictionary key and dense aggregations."""
        if len(self.keys) != 1:
            return None
        key_name = self.keys[0]
        key_info = self.table._schema.columns_by_name[key_name]
        key_is_dict = self.table._is_dictionary_column(key_info)
        if key_is_dict:
            key_arr = self.table._cols[key_name].codes
            key_dtype = np.dtype(np.int32)
            skip_key_null = self.dropna
            key_null = int(key_info.spec.null_code)
        else:
            key_arr = self.table._cols[key_name]
            key_dtype = getattr(key_info.spec, "dtype", None)
            if key_dtype is None:
                return None
            key_dtype = np.dtype(key_dtype)
            if key_dtype.kind not in "biu":
                return None
            key_null_value = getattr(key_info.spec, "null_value", None)
            skip_key_null = self.dropna and key_null_value is not None
            key_null = 0 if key_null_value is None else int(key_null_value)

        try:
            from blosc2 import groupby_ext
        except ImportError:
            return None

        descriptors = []
        for spec in specs:
            desc: dict[str, Any] = {"spec": spec, "op": spec.op}
            if spec.op == "size":
                kernel = getattr(groupby_ext, "groupby_dense_int_size_checked", None)
                if kernel is None:
                    return None
                desc.update({"kernel": kernel, "state_kind": "counts"})
                descriptors.append(desc)
                continue

            if spec.input_col is None:
                return None
            value_info = self.table._schema.columns_by_name[spec.input_col]
            value_dtype = getattr(value_info.spec, "dtype", None)
            if value_dtype is None:
                return None
            value_dtype = np.dtype(value_dtype)
            null_value = getattr(value_info.spec, "null_value", None)

            if spec.op == "count":
                kernel = getattr(groupby_ext, "groupby_dense_int_count_checked", None)
                if kernel is None:
                    return None
                desc.update({"kernel": kernel, "state_kind": "counts", "value_dtype": value_dtype})
            elif spec.op in {"sum", "mean", "min", "max"}:
                if value_dtype.kind == "f":
                    skip_nan = isinstance(null_value, float) and math.isnan(null_value)
                    if null_value is not None and not skip_nan:
                        return None
                    suffix = "sum" if spec.op == "sum" else spec.op
                    kernel = getattr(groupby_ext, f"groupby_dense_int_f64_{suffix}_checked", None)
                    if kernel is None:
                        return None
                    desc.update(
                        {
                            "kernel": kernel,
                            "value_dtype": np.float64,
                            "value_kind": "f64",
                            "skip_nan": skip_nan,
                        }
                    )
                elif value_dtype.kind in "biu":
                    if null_value is not None:
                        return None
                    if spec.op == "mean":
                        kernel = getattr(groupby_ext, "groupby_dense_int_f64_mean_checked", None)
                        if kernel is None:
                            return None
                        desc.update(
                            {
                                "kernel": kernel,
                                "value_dtype": np.float64,
                                "value_kind": "f64",
                                "skip_nan": False,
                            }
                        )
                    else:
                        kernel = getattr(groupby_ext, f"groupby_dense_int_i64_{spec.op}_checked", None)
                        if kernel is None:
                            return None
                        desc.update(
                            {
                                "kernel": kernel,
                                "value_dtype": np.int64,
                                "value_kind": "i64",
                                "skip_nan": False,
                            }
                        )
                else:
                    return None
                if spec.op in {"sum", "min", "max"}:
                    desc["state_kind"] = "value_present" if spec.op == "sum" else "extreme"
                elif spec.op == "mean":
                    desc["state_kind"] = "mean"
            else:
                return None
            descriptors.append(desc)

        compact_limit = 10_000_000
        keys_present = np.zeros(0, dtype=bool)
        states: dict[str, Any] = {}
        for desc in descriptors:
            spec = desc["spec"]
            if desc["state_kind"] == "counts":
                states[spec.output_col] = np.zeros(0, dtype=np.int64)
            elif desc["state_kind"] == "mean":
                states[spec.output_col] = (np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.int64))
            elif desc["state_kind"] == "value_present" or desc["state_kind"] == "extreme":
                dtype = np.float64 if desc["value_kind"] == "f64" else np.int64
                states[spec.output_col] = (np.zeros(0, dtype=dtype), np.zeros(0, dtype=bool))

        def ensure_size(size: int) -> bool:
            nonlocal keys_present, states
            if size > compact_limit:
                return False
            if size <= len(keys_present):
                return True
            old = len(keys_present)
            keys_present = np.pad(keys_present, (0, size - old), constant_values=False)
            for desc in descriptors:
                spec = desc["spec"]
                state = states[spec.output_col]
                if desc["state_kind"] == "counts":
                    states[spec.output_col] = np.pad(state, (0, size - old), constant_values=0)
                else:
                    first, second = state
                    states[spec.output_col] = (
                        np.pad(first, (0, size - old), constant_values=0),
                        np.pad(
                            second, (0, size - old), constant_values=False if second.dtype == np.bool_ else 0
                        ),
                    )
            return True

        def call_checked(kernel, *args) -> bool:
            return int(kernel(*args)) == 0

        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()
        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if not np.any(valid):
                continue
            keys = np.asarray(key_arr[start:stop], dtype=np.int8 if key_dtype.kind == "b" else key_dtype)
            keys = np.ascontiguousarray(keys)
            valid = np.ascontiguousarray(valid)
            live = valid.copy()
            if skip_key_null:
                live &= keys != key_null
            if not np.any(live):
                continue
            live_keys = keys[live]
            if np.min(live_keys) < 0:
                return None
            max_key = int(np.max(live_keys))
            if not ensure_size(max_key + 1):
                return None

            for desc in descriptors:
                spec = desc["spec"]
                state = states[spec.output_col]
                if spec.op == "size":
                    if not call_checked(
                        desc["kernel"], keys, valid, state, keys_present, skip_key_null, key_null
                    ):
                        return None
                elif spec.op == "count":
                    values = np.asarray(self.table._cols[spec.input_col][start:stop])
                    values_valid = np.ascontiguousarray(
                        ~self._null_mask(spec.input_col, values, is_key=False)
                    )
                    if not call_checked(
                        desc["kernel"],
                        keys,
                        valid,
                        values_valid,
                        state,
                        keys_present,
                        skip_key_null,
                        key_null,
                    ):
                        return None
                elif spec.op == "sum":
                    values = np.asarray(
                        self.table._cols[spec.input_col][start:stop], dtype=desc["value_dtype"]
                    )
                    values = np.ascontiguousarray(values)
                    sums, value_present = state
                    args = (
                        keys,
                        values,
                        valid,
                        sums,
                        value_present,
                        keys_present,
                        skip_key_null,
                        key_null,
                    )
                    if desc["value_kind"] == "f64":
                        args = (*args, desc["skip_nan"])
                    if not call_checked(desc["kernel"], *args):
                        return None
                elif spec.op == "mean":
                    values = np.asarray(
                        self.table._cols[spec.input_col][start:stop], dtype=desc["value_dtype"]
                    )
                    values = np.ascontiguousarray(values)
                    sums, counts = state
                    if not call_checked(
                        desc["kernel"],
                        keys,
                        values,
                        valid,
                        sums,
                        counts,
                        keys_present,
                        skip_key_null,
                        key_null,
                        desc["skip_nan"],
                    ):
                        return None
                elif spec.op in {"min", "max"}:
                    values = np.asarray(
                        self.table._cols[spec.input_col][start:stop], dtype=desc["value_dtype"]
                    )
                    values = np.ascontiguousarray(values)
                    extremes, has_value = state
                    args = (
                        keys,
                        values,
                        valid,
                        extremes,
                        has_value,
                        keys_present,
                        skip_key_null,
                        key_null,
                    )
                    if desc["value_kind"] == "f64":
                        args = (*args, desc["skip_nan"])
                    if not call_checked(desc["kernel"], *args):
                        return None

        # For integer keys np.nonzero already yields ascending key order.  When
        # sorting is requested (sort=True, or sort=None since the dict lexsort is
        # cheap) dictionary keys are reordered by decoded string; otherwise they
        # stay in first-appearance (code-assignment) order.
        group_codes = np.nonzero(keys_present)[0]
        if key_is_dict and self._resolve_sort(cheap=True):
            group_codes = self._sort_dict_group_codes(key_name, group_codes)

        key_values = (
            self.table._cols[key_name].decode_batch(group_codes)
            if key_is_dict
            else [_python_scalar(code) for code in group_codes]
        )
        rows = []
        for code, key_value in zip(group_codes, key_values, strict=True):
            row = {key_name: key_value}
            for desc in descriptors:
                spec = desc["spec"]
                state = states[spec.output_col]
                if spec.op in {"size", "count"}:
                    row[spec.output_col] = int(state[code])
                elif spec.op == "sum":
                    sums, value_present = state
                    row[spec.output_col] = (
                        _python_scalar(sums[code])
                        if value_present[code]
                        else _null_output_value(self._result_spec_for_agg(spec))
                    )
                elif spec.op == "mean":
                    sums, counts = state
                    row[spec.output_col] = (
                        math.nan if counts[code] == 0 else float(sums[code]) / int(counts[code])
                    )
                elif spec.op in {"min", "max"}:
                    extremes, has_value = state
                    row[spec.output_col] = (
                        _python_scalar(extremes[code])
                        if has_value[code]
                        else _null_output_value(self._result_spec_for_agg(spec))
                    )
            rows.append(row)
        return self._build_result(rows, specs)

    def _try_execute_cython_i32_f64_sum(self, specs: list[_AggSpec]):  # noqa: C901
        """Cython fast path for one int32 key and one non-null float64 sum."""
        if len(self.keys) != 1 or len(specs) != 1 or specs[0].op != "sum":
            return None
        spec = specs[0]
        if spec.input_col is None:
            return None
        key_name = self.keys[0]
        key_info = self.table._schema.columns_by_name[key_name]
        value_info = self.table._schema.columns_by_name[spec.input_col]
        if self.table._is_dictionary_column(key_info):
            key_arr = self.table._cols[key_name].codes
            key_is_dict = True
            key_null = int(key_info.spec.null_code)
            skip_key_null = self.dropna
        else:
            key_arr = self.table._cols[key_name]
            key_is_dict = False
            key_dtype = getattr(key_info.spec, "dtype", None)
            if key_dtype != np.dtype(np.int32):
                return None
            key_null_value = getattr(key_info.spec, "null_value", None)
            skip_key_null = self.dropna and key_null_value is not None
            key_null = 0 if key_null_value is None else int(key_null_value)
        value_dtype = getattr(value_info.spec, "dtype", None)
        if value_dtype != np.dtype(np.float64) or getattr(value_info.spec, "null_value", None) is not None:
            return None
        try:
            from blosc2 import groupby_ext
        except ImportError:
            return None
        kernel = getattr(groupby_ext, "groupby_dense_i32_f64_sum_checked", None)
        if kernel is None:
            return None

        compact_limit = 10_000_000
        sums = np.zeros(0, dtype=np.float64)
        present = np.zeros(0, dtype=bool)

        def ensure_size(size: int) -> bool:
            nonlocal sums, present
            if size > compact_limit:
                return False
            if size <= len(sums):
                return True
            old = len(sums)
            sums = np.pad(sums, (0, size - old), constant_values=0)
            present = np.pad(present, (0, size - old), constant_values=False)
            return True

        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()
        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if not np.any(valid):
                continue
            keys = np.asarray(key_arr[start:stop], dtype=np.int32)
            values = np.asarray(self.table._cols[spec.input_col][start:stop], dtype=np.float64)
            status = int(kernel(keys, values, valid, sums, present, skip_key_null, key_null, False))
            if status == -1:
                return None
            if status > 0:
                if not ensure_size(status):
                    return None
                status = int(kernel(keys, values, valid, sums, present, skip_key_null, key_null, False))
                if status != 0:
                    return None

        # np.nonzero gives ascending integer-key order; dict keys (which can
        # reach this path when the broader dense_int_key path declined) are
        # reordered by decoded string when a cheap sort is requested, matching
        # the other dense paths.
        present_codes = np.nonzero(present)[0]
        if key_is_dict and self._resolve_sort(cheap=True):
            present_codes = self._sort_dict_group_codes(key_name, present_codes)
        key_values = (
            self.table._cols[key_name].decode_batch(present_codes)
            if key_is_dict
            else [int(code) for code in present_codes]
        )
        rows = [
            {key_name: key_value, spec.output_col: float(sums[code])}
            for code, key_value in zip(present_codes, key_values, strict=True)
        ]
        return self._build_result(rows, specs)

    def _try_execute_cython_float_hash(self, specs: list[_AggSpec]):  # noqa: C901
        """Cython hash path for one arbitrary float key.

        This covers float32/float64 keys that are not suitable for dense
        integral-key indexing.  It currently supports float value columns for
        value reductions and falls back for unsupported mixed/multi-column cases.
        """
        if len(self.keys) != 1:
            return None
        key_name = self.keys[0]
        key_info = self.table._schema.columns_by_name[key_name]
        if self.table._is_dictionary_column(key_info):
            return None
        key_dtype = getattr(key_info.spec, "dtype", None)
        if key_dtype not in {np.dtype(np.float32), np.dtype(np.float64)}:
            return None

        value_cols = {s.input_col for s in specs if s.input_col is not None}
        if len(value_cols) > 1:
            return None
        value_col = next(iter(value_cols), None)
        value_dtype = None
        nullable_nan_value = False
        if value_col is not None:
            value_info = self.table._schema.columns_by_name[value_col]
            value_dtype = getattr(value_info.spec, "dtype", None)
            # Count can operate on any fixed-width value column via values_valid,
            # but other reductions in this hash kernel normalize values to f64.
            if any(s.op in {"sum", "mean", "min", "max"} for s in specs):
                if value_dtype is None or np.dtype(value_dtype).kind != "f":
                    return None
                null_value = getattr(value_info.spec, "null_value", None)
                nullable_nan_value = isinstance(null_value, float) and math.isnan(null_value)
                if null_value is not None and not nullable_nan_value:
                    return None

        try:
            from blosc2 import groupby_ext
        except ImportError:
            return None
        kernel = getattr(groupby_ext, "groupby_hash_f64_f64", None)
        if kernel is None:
            return None

        acc: dict[Any, dict[str, _AggState]] = {}
        key_values: dict[Any, tuple[Any, ...]] = {}
        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()

        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if not np.any(valid):
                continue
            keys = np.ascontiguousarray(np.asarray(self.table._cols[key_name][start:stop], dtype=np.float64))
            if value_col is None:
                values = np.empty(len(keys), dtype=np.float64)
                values_valid = np.zeros(len(keys), dtype=bool)
                has_values = False
            else:
                raw_values = np.asarray(self.table._cols[value_col][start:stop])
                if any(s.op in {"sum", "mean", "min", "max"} for s in specs):
                    values = np.ascontiguousarray(raw_values.astype(np.float64, copy=False))
                else:
                    values = np.empty(len(keys), dtype=np.float64)
                values_valid = np.ascontiguousarray(~self._null_mask(value_col, raw_values, is_key=False))
                has_values = True

            (
                chunk_keys,
                row_counts,
                value_counts,
                sums,
                mins,
                maxs,
                has_value,
            ) = kernel(keys, values, np.ascontiguousarray(valid), values_valid, has_values, self.dropna)

            for i, key in enumerate(chunk_keys):
                key_scalar = np.asarray(key, dtype=key_dtype).item()
                norm_key = _normalize_key_part(float(key_scalar))
                states = acc.setdefault(norm_key, {})
                key_values.setdefault(norm_key, (key_scalar,))
                for spec in specs:
                    state = states.setdefault(spec.output_col, _AggState(spec.op))
                    if spec.op == "size":
                        state.value = (0 if state.value is None else state.value) + int(row_counts[i])
                    elif spec.op == "count":
                        state.value = (0 if state.value is None else state.value) + int(value_counts[i])
                    elif spec.op == "sum" or spec.op == "mean":
                        if has_value[i]:
                            state.value = (0.0 if state.value is None else state.value) + float(sums[i])
                            state.count += int(value_counts[i])
                    elif spec.op == "min":
                        if has_value[i]:
                            value = float(mins[i])
                            if state.count == 0 or value < state.value:
                                state.value = value
                            state.count += 1
                    elif spec.op == "max" and has_value[i]:
                        value = float(maxs[i])
                        if state.count == 0 or value > state.value:
                            state.value = value
                        state.count += 1

        # list(acc) follows the Cython kernel's emission order (hash-bucket
        # order) -- deterministic for a given table but unspecified.  Key-sorting
        # it is a Python list.sort over all groups, so it only runs when
        # requested (sort=True; sort=None leaves it unsorted since this path is
        # not cheap to sort).
        ordered_keys = list(acc)
        if self._resolve_sort(cheap=False):
            ordered_keys.sort(
                key=lambda k: tuple(
                    (1, "") if isinstance(v, float) and math.isnan(v) else (0, v) for v in key_values[k]
                )
            )
        rows = []
        for norm_key in ordered_keys:
            row = dict(zip(self.keys, key_values[norm_key], strict=True))
            states = acc[norm_key]
            for spec in specs:
                state = states[spec.output_col]
                if spec.op == "mean":
                    row[spec.output_col] = math.nan if state.count == 0 else state.value / state.count
                elif spec.op in {"sum", "min", "max"} and state.count == 0:
                    row[spec.output_col] = _null_output_value(self._result_spec_for_agg(spec))
                else:
                    row[spec.output_col] = 0 if state.value is None else state.value
            rows.append(row)
        return self._build_result(rows, specs)

    def _try_execute_cython_float_integral_key_f64_sum(self, specs: list[_AggSpec]):  # noqa: C901
        """Cython fast path for integral float32/float64 keys and one non-null float64 sum."""
        if len(self.keys) != 1 or len(specs) != 1 or specs[0].op != "sum":
            return None
        spec = specs[0]
        if spec.input_col is None:
            return None
        key_name = self.keys[0]
        key_info = self.table._schema.columns_by_name[key_name]
        value_info = self.table._schema.columns_by_name[spec.input_col]
        key_dtype = getattr(key_info.spec, "dtype", None)
        value_dtype = getattr(value_info.spec, "dtype", None)
        if key_dtype not in {np.dtype(np.float32), np.dtype(np.float64)} or value_dtype != np.dtype(
            np.float64
        ):
            return None
        if getattr(value_info.spec, "null_value", None) is not None:
            return None
        # The fast path can skip NaNs.  If dropna=False and NaNs are present,
        # the Cython kernel reports unsupported and we fall back to generic
        # grouping, which can materialize a NaN group.
        skip_key_nan = self.dropna
        try:
            from blosc2 import groupby_ext
        except ImportError:
            return None
        kernel_name = (
            "groupby_dense_f32_integral_key_f64_sum_checked"
            if key_dtype == np.dtype(np.float32)
            else "groupby_dense_f64_integral_key_f64_sum_checked"
        )
        kernel = getattr(groupby_ext, kernel_name, None)
        if kernel is None:
            return None

        compact_limit = 10_000_000
        sums = np.zeros(0, dtype=np.float64)
        present = np.zeros(0, dtype=bool)

        def ensure_size(size: int) -> bool:
            nonlocal sums, present
            if size > compact_limit:
                return False
            if size <= len(sums):
                return True
            old = len(sums)
            sums = np.pad(sums, (0, size - old), constant_values=0)
            present = np.pad(present, (0, size - old), constant_values=False)
            return True

        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()
        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if not np.any(valid):
                continue
            keys = np.asarray(self.table._cols[key_name][start:stop], dtype=key_dtype)
            values = np.asarray(self.table._cols[spec.input_col][start:stop], dtype=np.float64)
            status = int(kernel(keys, values, valid, sums, present, skip_key_nan, False))
            if status == -1:
                return None
            if status > 0:
                if not ensure_size(status):
                    return None
                status = int(kernel(keys, values, valid, sums, present, skip_key_nan, False))
                if status != 0:
                    return None

        rows = [
            {key_name: float(code), spec.output_col: float(sums[code])} for code in np.nonzero(present)[0]
        ]
        return self._build_result(rows, specs)

    def _try_execute_dense_single_int_key(self, specs: list[_AggSpec]):  # noqa: C901
        """Fast path for one dense integer/dictionary-code key.

        This avoids per-chunk ``np.unique`` and Python dictionary merging.  It is
        intentionally conservative: keys must be non-negative and the observed
        key range must stay reasonably compact.
        """
        if len(self.keys) != 1:
            return None
        key_name = self.keys[0]
        key_info = self.table._schema.columns_by_name[key_name]
        key_is_dict = self.table._is_dictionary_column(key_info)
        key_dtype = np.dtype(np.int32) if key_is_dict else getattr(key_info.spec, "dtype", None)
        if key_dtype is None or key_dtype.kind not in "biuf":
            return None
        # Float keys are accepted only when their values are integral and fit a
        # compact non-negative range; per-chunk casting verifies integrality and
        # bails (falling back to the hash path) on the first fractional value.
        key_is_integral_float = (not key_is_dict) and key_dtype.kind == "f"
        extreme_ops = {"min", "max", "argmin", "argmax"}
        if any(spec.op in extreme_ops and spec.input_col is not None for spec in specs):
            for spec in specs:
                if spec.op in extreme_ops and spec.input_col is not None:
                    dtype = getattr(self.table._schema.columns_by_name[spec.input_col].spec, "dtype", None)
                    if dtype is None or np.dtype(dtype).kind not in "biufmM":
                        return None
        need_positions = any(spec.op in {"argmin", "argmax"} for spec in specs)

        compact_limit = 10_000_000
        present = np.zeros(0, dtype=bool)
        states: dict[str, Any] = {}
        for spec in specs:
            if spec.op in {"size", "count"}:
                states[spec.output_col] = np.zeros(0, dtype=np.int64)
            elif spec.op == "sum":
                out_dtype = np.int64
                if spec.input_col is not None:
                    dtype = np.dtype(self.table._schema.columns_by_name[spec.input_col].spec.dtype)
                    out_dtype = np.float64 if dtype.kind == "f" else np.int64
                states[spec.output_col] = (np.zeros(0, dtype=out_dtype), np.zeros(0, dtype=np.int64))
            elif spec.op == "mean":
                states[spec.output_col] = (np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.int64))
            elif spec.op in {"min", "max"}:
                assert spec.input_col is not None
                dtype = np.dtype(self.table._schema.columns_by_name[spec.input_col].spec.dtype)
                identity = _max_identity(dtype) if spec.op == "min" else _min_identity(dtype)
                states[spec.output_col] = (np.full(0, identity, dtype=dtype), np.zeros(0, dtype=bool))
            elif spec.op in {"argmin", "argmax"}:
                assert spec.input_col is not None
                dtype = np.dtype(self.table._schema.columns_by_name[spec.input_col].spec.dtype)
                identity = _max_identity(dtype) if spec.op == "argmin" else _min_identity(dtype)
                states[spec.output_col] = (
                    np.full(0, identity, dtype=dtype),  # best value seen per group
                    np.full(0, -1, dtype=np.int64),  # first row position attaining it
                    np.zeros(0, dtype=bool),  # group has any non-null value
                )

        def ensure_size(size: int) -> bool:
            nonlocal present, states
            if size > compact_limit:
                return False
            if size <= len(present):
                return True
            old = len(present)
            present = np.pad(present, (0, size - old), constant_values=False)
            for spec in specs:
                state = states[spec.output_col]
                if spec.op in {"size", "count"}:
                    states[spec.output_col] = np.pad(state, (0, size - old), constant_values=0)
                elif spec.op in {"sum", "mean"}:
                    sums, counts = state
                    states[spec.output_col] = (
                        np.pad(sums, (0, size - old), constant_values=0),
                        np.pad(counts, (0, size - old), constant_values=0),
                    )
                elif spec.op in {"min", "max"}:
                    values, has = state
                    dtype = values.dtype
                    identity = _max_identity(dtype) if spec.op == "min" else _min_identity(dtype)
                    states[spec.output_col] = (
                        np.pad(values, (0, size - old), constant_values=identity),
                        np.pad(has, (0, size - old), constant_values=False),
                    )
                elif spec.op in {"argmin", "argmax"}:
                    values, positions, has = state
                    dtype = values.dtype
                    identity = _max_identity(dtype) if spec.op == "argmin" else _min_identity(dtype)
                    states[spec.output_col] = (
                        np.pad(values, (0, size - old), constant_values=identity),
                        np.pad(positions, (0, size - old), constant_values=-1),
                        np.pad(has, (0, size - old), constant_values=False),
                    )
            return True

        phys_len = len(self.table._valid_rows)
        chunk_size = self._chunk_size()
        value_cols = sorted({s.input_col for s in specs if s.input_col is not None})
        logical_seen = 0  # running count of live rows, for argmin/argmax positions
        for start in range(0, phys_len, chunk_size):
            stop = min(start + chunk_size, phys_len)
            valid = np.asarray(self.table._valid_rows[start:stop], dtype=bool)
            if need_positions:
                logical_chunk = logical_seen + np.cumsum(valid, dtype=np.int64) - 1
                logical_seen += int(np.count_nonzero(valid))
            if not np.any(valid):
                continue
            raw_keys = self._read_key_chunk(key_name, start, stop)
            live_mask = valid.copy()
            if self.dropna:
                live_mask &= ~self._null_mask(key_name, raw_keys, is_key=True)
            if not np.any(live_mask):
                continue
            keys = np.asarray(raw_keys[live_mask])
            if keys.dtype.kind == "b":
                keys = keys.astype(np.int8, copy=False)
            elif key_is_integral_float:
                # Non-finite keys (NaN/inf, e.g. a retained NaN group when
                # dropna=False) and fractional keys both make the dense integer
                # mapping invalid; defer to the generic/hash fallback.
                if not np.all(np.isfinite(keys)):
                    return None
                keys_int = keys.astype(np.int64)
                if not np.array_equal(keys_int, keys):
                    return None
                keys = keys_int
            if len(keys) == 0:
                continue
            min_key = int(np.min(keys))
            if min_key < 0:
                return None
            max_key = int(np.max(keys))
            if not ensure_size(max_key + 1):
                return None
            present[keys] = True
            value_chunks = {
                name: np.asarray(self.table._cols[name][start:stop])[live_mask] for name in value_cols
            }
            row_positions = logical_chunk[live_mask] if need_positions else None

            for spec in specs:
                if spec.op == "size":
                    states[spec.output_col] += np.bincount(keys, minlength=len(present)).astype(np.int64)
                    continue
                assert spec.input_col is not None
                values = value_chunks[spec.input_col]
                non_null = ~self._null_mask(spec.input_col, values, is_key=False)
                if spec.op == "count":
                    states[spec.output_col] += np.bincount(
                        keys, weights=non_null.astype(np.int64), minlength=len(present)
                    ).astype(np.int64)
                elif spec.op == "sum":
                    sums, counts = states[spec.output_col]
                    if values.dtype.kind in "biu":
                        np.add.at(sums, keys[non_null], values[non_null].astype(np.int64, copy=False))
                    else:
                        sums += np.bincount(
                            keys, weights=np.where(non_null, values, 0), minlength=len(present)
                        ).astype(sums.dtype, copy=False)
                    counts += np.bincount(
                        keys, weights=non_null.astype(np.int64), minlength=len(present)
                    ).astype(np.int64)
                elif spec.op == "mean":
                    sums, counts = states[spec.output_col]
                    sums += np.bincount(keys, weights=np.where(non_null, values, 0), minlength=len(present))
                    counts += np.bincount(
                        keys, weights=non_null.astype(np.int64), minlength=len(present)
                    ).astype(np.int64)
                elif spec.op in {"min", "max"}:
                    values_state, has_state = states[spec.output_col]
                    if spec.op == "min":
                        np.minimum.at(values_state, keys[non_null], values[non_null])
                    else:
                        np.maximum.at(values_state, keys[non_null], values[non_null])
                    has_state[keys[non_null]] = True
                elif spec.op in {"argmin", "argmax"}:
                    best_state, pos_state, has_state = states[spec.output_col]
                    nstates = len(present)
                    k_nn = keys[non_null]
                    v_nn = values[non_null]
                    p_nn = row_positions[non_null]
                    # This chunk's per-group extreme and the first position attaining it.
                    if spec.op == "argmin":
                        chunk_best = np.full(
                            nstates, _max_identity(best_state.dtype), dtype=best_state.dtype
                        )
                        np.minimum.at(chunk_best, k_nn, v_nn)
                    else:
                        chunk_best = np.full(
                            nstates, _min_identity(best_state.dtype), dtype=best_state.dtype
                        )
                        np.maximum.at(chunk_best, k_nn, v_nn)
                    chunk_has = np.zeros(nstates, dtype=bool)
                    chunk_has[k_nn] = True
                    attains = v_nn == chunk_best[k_nn]
                    chunk_pos = np.full(nstates, np.iinfo(np.int64).max, dtype=np.int64)
                    np.minimum.at(chunk_pos, k_nn[attains], p_nn[attains])
                    # Merge: a strictly better value replaces the position; an equal
                    # value keeps the earlier (lower) position, so ties keep the first row.
                    if spec.op == "argmin":
                        better = chunk_has & (~has_state | (chunk_best < best_state))
                    else:
                        better = chunk_has & (~has_state | (chunk_best > best_state))
                    best_state[better] = chunk_best[better]
                    pos_state[better] = chunk_pos[better]
                    has_state |= chunk_has

        # np.nonzero gives ascending integer-key order.  When sorting is
        # requested (sort=True, or sort=None since the dict lexsort is cheap)
        # dictionary keys are reordered by decoded string; otherwise they stay in
        # first-appearance (code-assignment) order.
        group_codes = np.nonzero(present)[0]
        if key_is_dict and self._resolve_sort(cheap=True):
            group_codes = self._sort_dict_group_codes(key_name, group_codes)
        if key_is_dict:
            key_values = self.table._cols[key_name].decode_batch(group_codes)
        elif key_is_integral_float:
            key_values = [float(code) for code in group_codes]
        else:
            key_values = [_python_scalar(code) for code in group_codes]
        rows = []
        for code, key_value in zip(group_codes, key_values, strict=True):
            row = {key_name: key_value}
            for spec in specs:
                state = states[spec.output_col]
                if spec.op == "mean":
                    sums, counts = state
                    row[spec.output_col] = (
                        math.nan if counts[code] == 0 else float(sums[code]) / int(counts[code])
                    )
                elif spec.op == "sum":
                    sums, counts = state
                    row[spec.output_col] = (
                        _python_scalar(sums[code])
                        if counts[code] > 0
                        else _null_output_value(self._result_spec_for_agg(spec))
                    )
                elif spec.op in {"min", "max"}:
                    values_state, has_state = state
                    row[spec.output_col] = (
                        _python_scalar(values_state[code])
                        if has_state[code]
                        else _null_output_value(self._result_spec_for_agg(spec))
                    )
                elif spec.op in {"argmin", "argmax"}:
                    _best_state, pos_state, has_state = state
                    row[spec.output_col] = (
                        int(pos_state[code])
                        if has_state[code]
                        else _null_output_value(self._result_spec_for_agg(spec))
                    )
                else:
                    row[spec.output_col] = _python_scalar(state[code])
            rows.append(row)
        return self._build_result(rows, specs)

    def _chunk_size(self) -> int:
        if self.chunk_size is not None:
            if self.chunk_size <= 0:
                raise ValueError("chunk_size must be positive")
            return int(self.chunk_size)
        target = 1 << 20
        chunks = getattr(self.table._valid_rows, "chunks", None)
        if not chunks:
            return target
        base = max(int(chunks[0]), 1)
        if base >= target:
            return base
        # Batching the loop at the raw validity chunk shape can be pathological
        # (in-memory tables created small and grown by resize keep their tiny
        # initial chunk shape, e.g. 64 rows), and even 64 Ki-row batches leave
        # the per-batch Python bookkeeping (group display/merge) visible for
        # multi-key groupbys.  Scale up to ~1 Mi rows while staying
        # chunk-aligned; per-column batch memory stays modest (8 MB per int64
        # column).
        return base * -(-target // base)

    def _read_key_chunk(self, name: str, start: int, stop: int) -> np.ndarray:
        col_info = self.table._schema.columns_by_name[name]
        if self.table._is_dictionary_column(col_info):
            return np.asarray(self.table._cols[name].codes[start:stop], dtype=np.int32)
        if self.table._is_utf8_column(col_info):
            # Factorize the chunk from its raw offsets/bytes buffers: no row
            # is decoded, only the distinct values (codes flow through the
            # rest of the pipeline).  The factorizer is shared across chunks
            # so values seen before are hash-matched instead of re-sorted.
            # Utf8Array is sized to the logical row count, not the physical
            # valid_rows capacity, so a chunk boundary can run past its end;
            # rows beyond it are never live (the row can't have been written
            # without this column), so the padding code is never read live.
            col = self.table._cols[name]
            fact = self._utf8_factorizers.get(name)
            if fact is None:
                fact = self._utf8_factorizers[name] = col.factorizer()
            n = len(col)
            codes = fact.codes_for_span(start, min(stop, n)) if start < n else np.empty(0, dtype=np.int64)
            uniques = fact.uniques()
            # Rank-normalize so this chunk keeps the np.unique contract the
            # pipeline expects (uniques ascending, codes = string ranks).
            order = np.argsort(uniques, kind="stable")
            rank = np.empty(len(order), dtype=np.int64)
            rank[order] = np.arange(len(order))
            codes = rank[codes] if len(order) else codes
            if stop > n:
                codes = np.concatenate([codes, np.zeros(stop - max(start, n), dtype=np.int64)])
            return _Utf8KeyChunk(codes, uniques[order])
        return np.asarray(self.table._cols[name][start:stop])

    def _factorize_keys(
        self, keys_live: list[np.ndarray]
    ) -> tuple[np.ndarray | list[np.ndarray], np.ndarray]:
        if len(keys_live) == 1:
            arr = keys_live[0]
            if isinstance(arr, _Utf8KeyChunk):
                # The chunk is already factorized to dense string-rank codes;
                # dedupe them with an O(n) bincount instead of a sort.  The
                # np.unique contract — uniques ascending by string — holds
                # because the codes are ranks into the sorted uniques.
                present = np.bincount(arr.codes, minlength=len(arr.uniques)) > 0
                remap = np.cumsum(present) - 1
                return arr.uniques[present], remap[arr.codes]
            if arr.dtype.kind in ("U", "S") and arr.dtype.itemsize % 4 == 0 and arr.dtype.itemsize:
                return _factorize_fixed_width_str(arr)
            unique, inverse = np.unique(arr, return_inverse=True)
            return unique, inverse

        # utf8 key chunks pack as their int64 codes (codes are string-rank
        # within the chunk, so the packed sort order matches the string sort
        # order); the codes in the deduped rows are mapped back to strings
        # below, into object fields (StringDType cannot be a structured-array
        # field).
        pack_arrs = [arr.codes if isinstance(arr, _Utf8KeyChunk) else arr for arr in keys_live]
        composite = self._composite_int_factorize(pack_arrs)
        if composite is not None:
            unique_fields, inverse = composite
        else:
            dtype = [(f"k{i}", arr.dtype) for i, arr in enumerate(pack_arrs)]
            packed = np.empty(len(pack_arrs[0]), dtype=dtype)
            for i, arr in enumerate(pack_arrs):
                packed[f"k{i}"] = arr
            packed_unique, inverse = np.unique(packed, return_inverse=True)
            unique_fields = [packed_unique[f"k{i}"] for i in range(len(pack_arrs))]
        out_dtype = [
            (f"k{i}", object if isinstance(arr, _Utf8KeyChunk) else pack_arrs[i].dtype)
            for i, arr in enumerate(keys_live)
        ]
        unique = np.empty(len(unique_fields[0]), dtype=out_dtype)
        for i, arr in enumerate(keys_live):
            field = unique_fields[i]
            unique[f"k{i}"] = arr.uniques[field] if isinstance(arr, _Utf8KeyChunk) else field
        return unique, inverse

    @staticmethod
    def _composite_int_factorize(
        pack_arrs: list[np.ndarray],
    ) -> tuple[list[np.ndarray], np.ndarray] | None:
        """Dedup rows of all-integer key columns via one combined int64 key.

        ``np.unique`` over a structured dtype compares rows field-by-field
        through void comparisons and is ~an order of magnitude slower than
        over a plain int64 array, so when every key column is integral and
        the product of the per-column value ranges fits int64, combine the
        columns into a single integer (Horner over the zero-based fields —
        the combined sort order equals the structured lexicographic order)
        and dedup that.  Returns ``(per-field unique values, inverse)``, or
        ``None`` when the keys don't fit this scheme.
        """
        if not all(arr.dtype.kind in "biu" for arr in pack_arrs):
            return None
        if len(pack_arrs[0]) == 0:
            return None
        mins = [int(arr.min()) for arr in pack_arrs]
        spans = [int(arr.max()) - mn + 1 for arr, mn in zip(pack_arrs, mins, strict=True)]
        if math.prod(spans) >= 1 << 62:
            return None
        combined = np.zeros(len(pack_arrs[0]), dtype=np.int64)
        for arr, mn, span in zip(pack_arrs, mins, spans, strict=True):
            if arr.dtype.kind == "u":
                # Zero-base within the unsigned dtype first: the raw values may
                # not fit int64, but value - min always fits (span is checked).
                zero_based = (arr - arr.dtype.type(mn)).astype(np.int64, copy=False)
            else:
                zero_based = arr.astype(np.int64, copy=False) - mn
            combined *= span
            combined += zero_based
        unique_c, inverse = np.unique(combined, return_inverse=True)
        unique_fields: list[np.ndarray] = [None] * len(pack_arrs)  # type: ignore[list-item]
        rem = unique_c
        for i in range(len(pack_arrs) - 1, -1, -1):
            unique_fields[i] = (rem % spans[i] + mins[i]).astype(pack_arrs[i].dtype)
            rem = rem // spans[i]
        return unique_fields, inverse

    def _resolve_sort(self, *, cheap: bool) -> bool:
        """Resolve the tri-state ``self.sort`` request for the current path.

        ``cheap`` is declared by the call site: ``True`` for vectorized/free
        ordering (``np.nonzero`` ascending or the vectorized
        :meth:`_sort_dict_group_codes` lexsort), ``False`` for a Python
        ``list.sort`` over all groups.  ``None`` (auto) sorts only cheap paths.
        """
        return _resolve_sort(self.sort, cheap=cheap)

    def _sort_dict_group_codes(self, key_name: str, group_codes: np.ndarray) -> np.ndarray:
        """Reorder dictionary *group_codes* so groups come out sorted by string.

        Vectorized: the dictionary store holds only the unique category values
        (cardinality ``D``, not the ``N`` rows), so decompressing it whole into a
        NumPy array and ``np.lexsort``-ing the present codes is cheap — far
        cheaper than a Python ``sorted`` with a per-code ``decode`` call.  Null
        groups (``null_code``) sort first, matching ``_sortable_key_part``'s
        ordering of ``None`` before real values.
        """
        col = self.table._cols[key_name]
        null_code = int(col._spec.null_code)
        all_strings = np.asarray(col._dict_store[:])  # D unique values, no nulls
        is_non_null = group_codes != null_code
        decoded = np.empty(len(group_codes), dtype=all_strings.dtype if all_strings.size else "<U1")
        # Null groups keep a placeholder; the primary lexsort key sorts them
        # first regardless, so their string value is never compared.
        decoded[is_non_null] = all_strings[group_codes[is_non_null]]
        decoded[~is_non_null] = ""
        order = np.lexsort((decoded, is_non_null.astype(np.int8)))
        return group_codes[order]

    def _display_keys(self, unique_keys: np.ndarray | list[np.ndarray]) -> list[tuple[Any, ...]]:
        if len(self.keys) == 1:
            name = self.keys[0]
            col_info = self.table._schema.columns_by_name[name]
            unique_arr = np.asarray(unique_keys)
            if self.table._is_dictionary_column(col_info):
                decoded = self.table._cols[name].decode_batch(unique_arr)
                return [(value,) for value in decoded]
            return [(_python_scalar(value),) for value in unique_arr]

        result = []
        assert isinstance(unique_keys, np.ndarray)
        for row in unique_keys:
            vals = []
            for i, name in enumerate(self.keys):
                value = row[f"k{i}"]
                col_info = self.table._schema.columns_by_name[name]
                if self.table._is_dictionary_column(col_info):
                    vals.append(self.table._cols[name].decode(int(value)))
                else:
                    vals.append(_python_scalar(value))
            result.append(tuple(vals))
        return result

    def _normalized_keys(self, display_keys: list[tuple[Any, ...]]) -> list[Any]:
        normalized = []
        for key in display_keys:
            norm = tuple(_normalize_key_part(v) for v in key)
            normalized.append(norm[0] if len(norm) == 1 else norm)
        return normalized

    def _compute_partials(
        self,
        specs: list[_AggSpec],
        unique_keys: np.ndarray | list[np.ndarray],
        inverse: np.ndarray,
        value_chunks: dict[str, np.ndarray],
        row_positions: np.ndarray,
    ) -> dict[str, Any]:
        n_groups = len(unique_keys)
        partials: dict[str, Any] = {}
        for spec in specs:
            if spec.op == "size":
                partials[spec.output_col] = np.bincount(inverse, minlength=n_groups).astype(np.int64)
                continue

            assert spec.input_col is not None
            values = value_chunks[spec.input_col]
            non_null = ~self._null_mask(spec.input_col, values, is_key=False)

            if spec.op == "count":
                partials[spec.output_col] = np.bincount(
                    inverse, weights=non_null.astype(np.int64), minlength=n_groups
                ).astype(np.int64)
            elif spec.op in {"sum", "mean"}:
                counts = np.bincount(inverse, weights=non_null.astype(np.int64), minlength=n_groups).astype(
                    np.int64
                )
                if spec.op == "sum" and values.dtype.kind in "biu":
                    sums = np.zeros(n_groups, dtype=np.int64)
                    np.add.at(sums, inverse[non_null], values[non_null].astype(np.int64, copy=False))
                else:
                    weights = np.where(non_null, values, 0)
                    sums = np.bincount(inverse, weights=weights, minlength=n_groups)
                partials[spec.output_col] = (sums, counts)
            elif spec.op in {"min", "max"}:
                partials[spec.output_col] = self._minmax_partials(
                    spec.op, inverse, values, non_null, n_groups
                )
            elif spec.op in {"argmin", "argmax"}:
                partials[spec.output_col] = self._argminmax_partials(
                    spec.op, inverse, values, non_null, row_positions, n_groups
                )
            elif spec.op == "udf":
                partials[spec.output_col] = self._udf_value_partials(inverse, values, non_null, n_groups)
        return partials

    def _udf_value_partials(
        self, inverse: np.ndarray, values: np.ndarray, non_null: np.ndarray, n_groups: int
    ) -> list[np.ndarray]:
        """Split this chunk's non-null values by group, for UDF aggregations.

        Unlike the built-in ops, an arbitrary Python callable cannot be
        incrementally merged across chunks, so each chunk instead contributes
        its raw per-group values; :meth:`_merge_partials` collects these into
        a growing list per group, and :meth:`_final_rows` concatenates and
        calls the UDF once, after all chunks have been read.
        """
        groups = inverse[non_null]
        vals = values[non_null]
        order = np.argsort(groups, kind="stable")
        sorted_groups = groups[order]
        sorted_vals = vals[order]
        boundaries = np.searchsorted(sorted_groups, np.arange(n_groups + 1))
        return [sorted_vals[boundaries[g] : boundaries[g + 1]] for g in range(n_groups)]

    def _minmax_partials(
        self, op: AggName, inverse: np.ndarray, values: np.ndarray, non_null: np.ndarray, n_groups: int
    ) -> tuple[np.ndarray, np.ndarray]:
        if values.dtype.kind in "biufcmM":
            if op == "min":
                identity = _max_identity(values.dtype)
                out = np.full(n_groups, identity, dtype=values.dtype)
                np.minimum.at(out, inverse[non_null], values[non_null])
            else:
                identity = _min_identity(values.dtype)
                out = np.full(n_groups, identity, dtype=values.dtype)
                np.maximum.at(out, inverse[non_null], values[non_null])
        else:
            out = np.empty(n_groups, dtype=values.dtype)
            has = np.zeros(n_groups, dtype=bool)
            for group, value, ok in zip(inverse, values, non_null, strict=True):
                if not ok:
                    continue
                if not has[group] or (value < out[group] if op == "min" else value > out[group]):
                    out[group] = value
                    has[group] = True
            return out, has
        has_value = np.bincount(inverse, weights=non_null.astype(np.int64), minlength=n_groups) > 0
        return out, has_value

    def _argminmax_partials(
        self,
        op: AggName,
        inverse: np.ndarray,
        values: np.ndarray,
        non_null: np.ndarray,
        row_positions: np.ndarray,
        n_groups: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Reduce the per-group extreme value with the (vectorized) min/max path,
        # then pick the *first* row position that attains it — vectorized, instead
        # of the former per-row Python loop (~30x faster on millions of rows).
        mm_op = "min" if op == "argmin" else "max"
        best_values, has_value = self._minmax_partials(mm_op, inverse, values, non_null, n_groups)
        best_positions = np.full(n_groups, -1, dtype=np.int64)
        # Rows in a value-bearing group whose value equals their group's extreme.
        # On ties np.minimum.at keeps the smallest position, matching "first row".
        attains = non_null & has_value[inverse] & (values == best_values[inverse])
        if np.any(attains):
            int_max = np.iinfo(np.int64).max
            pos_best = np.full(n_groups, int_max, dtype=np.int64)
            np.minimum.at(pos_best, inverse[attains], np.asarray(row_positions)[attains])
            best_positions = np.where(pos_best != int_max, pos_best, -1)
        return best_values, best_positions, has_value

    def _merge_partials(  # noqa: C901
        self,
        acc: dict[Any, dict[str, _AggState]],
        key_values: dict[Any, tuple[Any, ...]],
        normalized_keys: list[Any],
        display_keys: list[tuple[Any, ...]],
        partials: dict[str, Any],
        specs: list[_AggSpec],
    ) -> None:
        for i, norm_key in enumerate(normalized_keys):
            states = acc.setdefault(norm_key, {})
            key_values.setdefault(norm_key, display_keys[i])
            for spec in specs:
                state = states.setdefault(spec.output_col, _AggState(spec.op))
                partial = partials[spec.output_col]
                if spec.op in {"size", "count"}:
                    state.value = (0 if state.value is None else state.value) + int(partial[i])
                elif spec.op == "sum":
                    sums, counts = partial
                    if counts[i] > 0:
                        state.value = (0 if state.value is None else state.value) + _python_scalar(sums[i])
                        state.count += int(counts[i])
                elif spec.op == "mean":
                    sums, counts = partial
                    if counts[i] > 0:
                        state.value = (0.0 if state.value is None else state.value) + float(sums[i])
                        state.count += int(counts[i])
                elif spec.op in {"min", "max"}:
                    values, has_value = partial
                    if has_value[i]:
                        value = _python_scalar(values[i])
                        if (
                            state.count == 0
                            or (spec.op == "min" and value < state.value)
                            or (spec.op == "max" and value > state.value)
                        ):
                            state.value = value
                        state.count += 1
                elif spec.op in {"argmin", "argmax"}:
                    values, positions, has_value = partial
                    if has_value[i]:
                        value = _python_scalar(values[i])
                        if (
                            state.count == 0
                            or (spec.op == "argmin" and value < state.value[0])
                            or (spec.op == "argmax" and value > state.value[0])
                        ):
                            state.value = (value, int(positions[i]))
                        state.count += 1
                elif spec.op == "udf":
                    chunk_values = partial[i]
                    if state.value is None:
                        state.value = []
                    if len(chunk_values):
                        state.value.append(chunk_values)

    def _final_rows(  # noqa: C901
        self,
        acc: dict[Any, dict[str, _AggState]],
        key_values: dict[Any, tuple[Any, ...]],
        specs: list[_AggSpec],
    ) -> list[dict[str, Any]]:
        # This path handles keys that miss the vectorized fast paths (fixed-width
        # strings, multi-key, exotic dtypes); ordering is a Python list.sort over
        # all groups, so it only runs when requested.  list(acc) is deterministic
        # first-appearance order (sort=False, or sort=None since this path is not
        # cheap to sort).
        keys = list(acc)
        if self._resolve_sort(cheap=False):
            keys.sort(key=lambda k: tuple(_sortable_key_part(v) for v in key_values[k]))

        # Sentinel marking a UDF aggregation group with no real result yet --
        # either it had zero non-null input values (the UDF was never
        # called), or the UDF itself returned None to signal a null result.
        # Both cases are excluded from dtype inference and patched to the
        # output dtype's null value once it is known, below.
        _empty_udf_group = object()

        udf_results: dict[str, list] = {spec.output_col: [] for spec in specs if spec.op == "udf"}
        rows = []
        for norm_key in keys:
            row = dict(zip(self.keys, key_values[norm_key], strict=True))
            states = acc[norm_key]
            for spec in specs:
                state = states[spec.output_col]
                if spec.op == "mean":
                    row[spec.output_col] = math.nan if state.count == 0 else state.value / state.count
                elif spec.op == "udf":
                    chunks = state.value
                    if not chunks:
                        # No non-null values for this group/column; matches
                        # the sum/min/max convention of a null result rather
                        # than calling the UDF with an empty array.
                        row[spec.output_col] = _empty_udf_group
                    else:
                        group_values = np.concatenate(chunks)
                        # A @blosc2.dsl_kernel-decorated UDF is a DSLKernel
                        # instance whose __call__ expects the array-kernel
                        # calling convention (inputs_tuple, output, offset),
                        # not this "one array in, one scalar out" aggregation
                        # convention -- call the wrapped plain function instead.
                        udf_callable = spec.udf.func if isinstance(spec.udf, DSLKernel) else spec.udf
                        try:
                            result = _python_scalar(udf_callable(group_values))
                        except Exception as exc:
                            raise RuntimeError(
                                f"UDF aggregation {spec.output_col!r} raised for group "
                                f"{key_values[norm_key]!r}: {exc}"
                            ) from exc
                        if result is None:
                            # The UDF itself signaled a null result; treat it
                            # like an empty group rather than feeding None
                            # into dtype inference (which would otherwise see
                            # it as a genuinely inconsistent result type).
                            row[spec.output_col] = _empty_udf_group
                        else:
                            row[spec.output_col] = result
                            udf_results[spec.output_col].append(result)
                elif spec.op in {"sum", "min", "max", "argmin", "argmax"} and state.count == 0:
                    row[spec.output_col] = _null_output_value(self._result_spec_for_agg(spec))
                elif spec.op in {"argmin", "argmax"}:
                    row[spec.output_col] = state.value[1]
                else:
                    row[spec.output_col] = 0 if state.value is None else state.value
            rows.append(row)
        for spec in specs:
            if spec.op != "udf":
                continue
            if spec.explicit_dtype is None:
                spec.explicit_dtype = self._infer_udf_spec(udf_results[spec.output_col], spec.output_col)
            null_value = _null_output_value(spec.explicit_dtype)
            for row in rows:
                if row[spec.output_col] is _empty_udf_group:
                    row[spec.output_col] = null_value
        return rows

    @staticmethod
    def _infer_udf_spec(results: list, name: str) -> SchemaSpec:
        """Infer a CTable schema spec from a UDF aggregation's collected results.

        Probing every group's result (not just the first) means a UDF that
        returns inconsistent types (e.g. int for one group, a string for
        another) is caught here with a clear error, rather than surfacing as
        an opaque failure while building the result table.
        """
        if not results:
            # No group ever produced a value (e.g. an empty table, or every
            # group is all-null so the UDF was never called) -- there is
            # nothing to infer a dtype from.
            raise ValueError(
                f"Cannot infer a CTable dtype for UDF aggregation {name!r}: it was never "
                f"called (empty table, or every group had no non-null values). Pass an "
                f"explicit dtype in the named-agg tuple, e.g. "
                f"g.agg({name}=(col, fn, blosc2.float64()))."
            )
        try:
            arr = np.asarray(results)
        except ValueError as exc:
            # Ragged/inhomogeneous results (e.g. a UDF returning a list for
            # one group and a scalar for another) raise here rather than
            # producing an object array.
            raise ValueError(
                f"UDF aggregation {name!r} produced inconsistent or unsupported types "
                f"across groups: {results!r} ({exc}). Pass an explicit dtype in the "
                f"named-agg tuple, e.g. g.agg({name}=(col, fn, blosc2.float64()))."
            ) from exc
        if arr.dtype == object:
            raise ValueError(
                f"UDF aggregation {name!r} produced inconsistent or unsupported types "
                f"across groups: {results!r}. Pass an explicit dtype in the named-agg "
                f"tuple, e.g. g.agg({name}=(col, fn, blosc2.float64()))."
            )
        if arr.dtype.kind == "b":
            return b2_bool()
        if arr.dtype.kind in "iu":
            return int64()
        if arr.dtype.kind == "f":
            return float64()
        raise ValueError(
            f"Cannot infer a CTable dtype for UDF aggregation {name!r} result dtype "
            f"{arr.dtype!r}. Pass an explicit dtype in the named-agg tuple, e.g. "
            f"g.agg({name}=(col, fn, blosc2.string(max_length=32)))."
        )

    def _build_result(self, rows: list[dict[str, Any]], specs: list[_AggSpec]):
        from blosc2.ctable import CTable

        columns = self.keys + [spec.output_col for spec in specs]
        schema_specs = {name: self._result_spec_for_key(name) for name in self.keys}
        for spec in specs:
            schema_specs[spec.output_col] = self._result_spec_for_agg(spec)

        # CTable construction is schema-from-dataclass based, and dataclass field
        # names must be Python identifiers.  Build with temporary identifier
        # aliases, then rename the result columns back to their requested logical
        # names.  This keeps nested names such as ``trip.sec`` in the resulting
        # table while reusing the normal CTable initialization path.
        alias_by_name = self._result_aliases(columns)
        fields = []
        for name in columns:
            alias = alias_by_name[name]
            fields.append((alias, _python_type_for_spec(schema_specs[name]), b2_field(schema_specs[name])))
        row_type = dataclasses.make_dataclass("CTableGroupByRow", fields)
        data = {alias_by_name[name]: [row[name] for row in rows] for name in columns}
        urlpath = getattr(self, "_result_urlpath", None)
        kwargs = {"urlpath": str(urlpath), "mode": "w"} if urlpath is not None else {}
        out = CTable(row_type, new_data=data, expected_size=max(len(rows), 1), validate=False, **kwargs)
        for name in columns:
            alias = alias_by_name[name]
            if alias != name:
                out.rename_column(alias, name)
        return out

    def _result_aliases(self, names: Sequence[str]) -> dict[str, str]:
        final_names = set(names)
        aliases: dict[str, str] = {}
        used: set[str] = set()
        for i, name in enumerate(names):
            alias = f"groupby_result_col_{i}"
            suffix = 0
            while alias in final_names or alias in used:
                suffix += 1
                alias = f"groupby_result_col_{i}_{suffix}"
            aliases[name] = alias
            used.add(alias)
        return aliases

    def _validate_output_names(self, specs: list[_AggSpec]) -> None:
        names = self.keys + [s.output_col for s in specs]
        for name in names:
            _validate_column_name(name)
        if len(names) != len(set(names)):
            raise ValueError("Group-by result column names would not be unique")

    def _result_spec_for_key(self, name: str) -> SchemaSpec:
        return copy.deepcopy(self.table._schema.columns_by_name[name].spec)

    def _result_spec_for_agg(self, spec: _AggSpec) -> SchemaSpec:
        if spec.op in {"size", "count"}:
            return int64()
        if spec.op in {"argmin", "argmax"}:
            return int64(null_value=-1)
        if spec.op == "mean":
            return float64()
        if spec.op == "udf":
            # _final_rows() always sets this (inferred, unless the user gave
            # an explicit dtype) before _build_result() reads it.
            assert spec.explicit_dtype is not None
            return spec.explicit_dtype
        assert spec.input_col is not None
        input_spec = self.table._schema.columns_by_name[spec.input_col].spec
        dtype = getattr(input_spec, "dtype", None)
        if spec.op == "sum":
            if dtype is not None and dtype.kind in "iu":
                return int64()
            if dtype is not None and dtype.kind == "b":
                return int64()
            if dtype is not None and dtype.kind == "f":
                return float64()
        return copy.deepcopy(input_spec)

    def _null_mask(self, name: str, values: np.ndarray, *, is_key: bool) -> np.ndarray:
        col_info = self.table._schema.columns_by_name[name]
        spec = col_info.spec
        if isinstance(values, _Utf8KeyChunk):
            null_value = getattr(spec, "null_value", None)
            if null_value is None:
                return np.zeros(len(values), dtype=bool)
            return values.codes == values.code_of(null_value)
        if isinstance(spec, DictionarySpec):
            mask = values == np.int32(spec.null_code)
            return mask if is_key or getattr(spec, "nullable", False) else np.zeros(len(values), dtype=bool)
        null_value = getattr(spec, "null_value", None)
        mask = np.zeros(len(values), dtype=bool)
        # For keys, treat all NaNs as missing so dropna behaves predictably.
        # For values, only nullable NaN sentinels are skipped.
        if values.dtype.kind == "f" and (
            is_key or (isinstance(null_value, float) and math.isnan(null_value))
        ):
            mask |= np.isnan(values)
        if null_value is not None and not (isinstance(null_value, float) and math.isnan(null_value)):
            mask |= values == null_value
        return mask


_HASH_MIX = np.uint64(0x9E3779B97F4A7C15)


def _factorize_fixed_width_str(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact ``np.unique(arr, return_inverse=True)`` for fixed-width string
    keys, ~2x faster.

    Argsorting N strings is the string-key groupby bottleneck (a ``U8`` key
    compares 32 bytes of UTF-32 per element), so instead: hash each row's raw
    bytes into one uint64, factorize the *integers*, and recover the string
    for each group from one representative row.  A vectorized verify pass
    keeps it exact — on a hash collision (different strings, same hash) fall
    back to plain ``np.unique``.  Output contract (uniques sorted ascending,
    inverse indices into them) is identical to ``np.unique``, so callers
    cannot tell the difference.

    ponytail: per-chunk cost is now the int64 unique-argsort; a cross-chunk
    vocabulary cache (searchsorted against known hashes) could roughly halve
    it again if string-key groupby speed ever matters more.
    """
    words = arr.view(np.uint32).reshape(len(arr), arr.dtype.itemsize // 4)
    h = words[:, 0].astype(np.uint64)
    for i in range(1, words.shape[1]):
        h = (h * _HASH_MIX) ^ words[:, i]
    hash_uniques, inverse = np.unique(h, return_inverse=True)
    representative = np.empty(len(hash_uniques), dtype=np.int64)
    representative[inverse] = np.arange(len(arr))  # any row of the group will do
    reps = arr[representative]
    if not (arr == reps[inverse]).all():
        return np.unique(arr, return_inverse=True)  # collision: exact fallback
    # np.unique contract: uniques ascending by *string*, not by hash.
    order = np.argsort(reps, kind="stable")
    rank = np.empty(len(order), dtype=inverse.dtype)
    rank[order] = np.arange(len(order), dtype=inverse.dtype)
    return reps[order], rank[inverse]


def _normalize_key_part(value: Any) -> Any:
    if isinstance(value, float) and math.isnan(value):
        return _NAN_KEY
    return value


def _sortable_key_part(value: Any) -> tuple[int, Any]:
    if value is None:
        return (0, "")
    if isinstance(value, float) and math.isnan(value):
        return (0, "")
    return (1, value)


def _python_scalar(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    return value


def _python_type_for_spec(spec: SchemaSpec):
    if isinstance(spec, DictionarySpec):
        return str
    if isinstance(spec, b2_bool):
        return bool
    dtype = getattr(spec, "dtype", None)
    if dtype is not None:
        if dtype.kind in "iu":
            return int
        if dtype.kind == "f":
            return float
        if dtype.kind == "b":
            return bool
        if dtype.kind in "US":
            return str if dtype.kind == "U" else bytes
    return getattr(spec, "python_type", object)


def _max_identity(dtype: np.dtype):
    dtype = np.dtype(dtype)
    if dtype.kind in "iu":
        return np.iinfo(dtype).max
    if dtype.kind == "f":
        return np.inf
    if dtype.kind in "mM":
        return np.iinfo(np.int64).max
    return None


def _min_identity(dtype: np.dtype):
    dtype = np.dtype(dtype)
    if dtype.kind in "iu":
        return np.iinfo(dtype).min
    if dtype.kind == "f":
        return -np.inf
    if dtype.kind in "mM":
        return np.iinfo(np.int64).min
    return None


def _null_output_value(spec: SchemaSpec):
    dtype = getattr(spec, "dtype", None)
    null_value = getattr(spec, "null_value", None)
    if null_value is not None:
        return null_value
    if dtype is not None and dtype.kind == "f":
        return math.nan
    if dtype is not None and dtype.kind in "iu":
        return 0
    if dtype is not None and dtype.kind == "b":
        return False
    if dtype is not None and dtype.kind == "U":
        return ""
    if dtype is not None and dtype.kind == "S":
        return b""
    return None


# ----------------------------------------------------------------------
# Public array-oriented grouped reductions
# ----------------------------------------------------------------------


def group_reduce(keys, values=None, op: AggName = "size", *, sort: bool | None = None, dropna: bool = True):
    """Group *keys* and reduce *values* with *op*.

    This is a lower-level, array-oriented grouped reduction primitive.  It exposes
    Blosc2's optimized group-reduce kernels for one-dimensional array-like inputs,
    including NumPy arrays and :class:`blosc2.NDArray`, without requiring a
    :class:`blosc2.CTable`.

    Parameters
    ----------
    keys : array-like
        One-dimensional grouping keys.
    values : array-like, optional
        One-dimensional values to reduce.  Required for ``"count"``, ``"sum"``,
        ``"mean"``, ``"min"`` and ``"max"``.  Ignored for ``"size"``.
    op : {"size", "count", "sum", "mean", "min", "max"}, default: "size"
        Reduction operation.  ``"size"`` counts rows per group, while
        ``"count"`` counts non-NaN values per group.
    sort : bool or None, default: None
        Output group ordering.  ``True`` always sorts groups by key; ``False``
        never sorts (order is implementation dependent but deterministic for a
        given input); ``None`` (the default) sorts only when cheap -- the
        vectorized integer and float kernels sort (``np.argsort``), while the
        generic Python fallback is left unsorted to avoid an O(G log G) Python
        sort that can rival the grouping cost on high-cardinality keys.
    dropna : bool, default: True
        If true, skip NaN float keys.  If false, all NaN keys form one group.

    Returns
    -------
    groups, result : numpy.ndarray, numpy.ndarray
        Group keys and reduced values.

    Examples
    --------
    >>> import blosc2
    >>> keys = blosc2.array([1, 2, 1, 2, 1])
    >>> values = blosc2.array([10., 20., 30., 40., 50.])
    >>> groups, sums = blosc2.group_reduce(keys, values, op="sum", sort=True)
    >>> groups
    array([1, 2])
    >>> sums
    array([90., 60.])
    """
    if op not in {"size", "count", "sum", "mean", "min", "max"}:
        raise ValueError(f"unsupported group_reduce operation {op!r}")

    keys_arr = np.asarray(keys)
    if keys_arr.ndim != 1:
        raise ValueError("keys must be a 1-D array")

    if op == "size":
        values_arr = None
    else:
        if values is None:
            raise ValueError(f"values are required for group_reduce op {op!r}")
        values_arr = np.asarray(values)
        if values_arr.ndim != 1:
            raise ValueError("values must be a 1-D array")
        if len(values_arr) != len(keys_arr):
            raise ValueError("keys and values must have the same length")

    if len(keys_arr) == 0:
        return keys_arr.copy(), np.empty(0, dtype=_result_dtype(values_arr, op))

    # The dense integer and float-hash kernels sort via vectorized np.argsort
    # (cheap); the generic Python fallback sorts via list.sort (expensive).
    # Resolve the tri-state per path's cost so None auto-sorts only the cheap ones.
    fast = _try_dense_integer(keys_arr, values_arr, op, sort=_resolve_sort(sort, cheap=True))
    if fast is not None:
        return fast

    fast = _try_float_hash(keys_arr, values_arr, op, sort=_resolve_sort(sort, cheap=True), dropna=dropna)
    if fast is not None:
        return fast

    return _group_reduce_numpy(
        keys_arr, values_arr, op, sort=_resolve_sort(sort, cheap=False), dropna=dropna
    )


def _try_dense_integer(keys: np.ndarray, values: np.ndarray | None, op: str, *, sort: bool):  # noqa: C901
    key_dtype = np.dtype(keys.dtype)
    if key_dtype.kind == "b":
        keys = keys.astype(np.int8, copy=False)
    elif key_dtype.kind not in "iu":
        return None
    keys = np.ascontiguousarray(keys)
    if len(keys) == 0:
        return None
    if np.min(keys) < 0:
        return None
    max_key = int(np.max(keys))
    if max_key + 1 > 10_000_000:
        return None

    try:
        from blosc2 import groupby_ext
    except ImportError:
        return None

    valid = np.ones(len(keys), dtype=bool)
    keys_present = np.zeros(max_key + 1, dtype=bool)

    if op == "size":
        kernel = getattr(groupby_ext, "groupby_dense_int_size_checked", None)
        if kernel is None:
            return None
        counts = np.zeros(max_key + 1, dtype=np.int64)
        kernel(keys, valid, counts, keys_present, False, 0)
        groups = np.nonzero(keys_present)[0].astype(key_dtype if key_dtype.kind != "b" else np.bool_)
        result = counts[np.nonzero(keys_present)[0]]
        return _maybe_sort(groups, result, sort)

    assert values is not None
    value_dtype = np.dtype(values.dtype)
    if op == "count":
        kernel = getattr(groupby_ext, "groupby_dense_int_count_checked", None)
        if kernel is None:
            return None
        counts = np.zeros(max_key + 1, dtype=np.int64)
        values_valid = _values_valid(values)
        kernel(keys, valid, np.ascontiguousarray(values_valid), counts, keys_present, False, 0)
        codes = np.nonzero(keys_present)[0]
        return _maybe_sort(
            codes.astype(key_dtype if key_dtype.kind != "b" else np.bool_), counts[codes], sort
        )

    if op == "mean" or value_dtype.kind == "f":
        vals = np.ascontiguousarray(values.astype(np.float64, copy=False))
        skip_nan = value_dtype.kind == "f"
        if op == "sum":
            kernel = getattr(groupby_ext, "groupby_dense_int_f64_sum_checked", None)
            if kernel is None:
                return None
            sums = np.zeros(max_key + 1, dtype=np.float64)
            present = np.zeros(max_key + 1, dtype=bool)
            kernel(keys, vals, valid, sums, present, keys_present, False, 0, skip_nan)
            codes = np.nonzero(keys_present)[0]
            result = sums[codes]
            result[~present[codes]] = np.nan
        elif op == "mean":
            kernel = getattr(groupby_ext, "groupby_dense_int_f64_mean_checked", None)
            if kernel is None:
                return None
            sums = np.zeros(max_key + 1, dtype=np.float64)
            counts = np.zeros(max_key + 1, dtype=np.int64)
            kernel(keys, vals, valid, sums, counts, keys_present, False, 0, skip_nan)
            codes = np.nonzero(keys_present)[0]
            result = np.full(len(codes), np.nan, dtype=np.float64)
            ok = counts[codes] > 0
            result[ok] = sums[codes][ok] / counts[codes][ok]
        elif op in {"min", "max"}:
            state = np.zeros(max_key + 1, dtype=np.float64)
            has_value = np.zeros(max_key + 1, dtype=bool)
            kernel = getattr(groupby_ext, f"groupby_dense_int_f64_{op}_checked", None)
            if kernel is None:
                return None
            kernel(keys, vals, valid, state, has_value, keys_present, False, 0, skip_nan)
            codes = np.nonzero(keys_present)[0]
            result = state[codes]
            result[~has_value[codes]] = np.nan
        else:  # pragma: no cover
            return None
        return _maybe_sort(codes.astype(key_dtype if key_dtype.kind != "b" else np.bool_), result, sort)

    if value_dtype.kind not in "biu":
        return None
    vals_i64 = np.ascontiguousarray(values.astype(np.int64, copy=False))
    state = np.zeros(max_key + 1, dtype=np.int64)
    present = np.zeros(max_key + 1, dtype=bool)
    kernel = getattr(groupby_ext, f"groupby_dense_int_i64_{op}_checked", None)
    if kernel is None:
        return None
    kernel(keys, vals_i64, valid, state, present, keys_present, False, 0)
    codes = np.nonzero(keys_present)[0]
    return _maybe_sort(codes.astype(key_dtype if key_dtype.kind != "b" else np.bool_), state[codes], sort)


def _try_float_hash(keys: np.ndarray, values: np.ndarray | None, op: str, *, sort: bool, dropna: bool):
    key_dtype = np.dtype(keys.dtype)
    if key_dtype.kind != "f":
        return None
    if values is not None and np.dtype(values.dtype).kind != "f" and op != "count":
        return None
    try:
        from blosc2 import groupby_ext
    except ImportError:
        return None

    keys_f64 = np.ascontiguousarray(keys.astype(np.float64, copy=False))
    valid = np.ones(len(keys_f64), dtype=bool)
    if values is None:
        values_f64 = np.empty(len(keys_f64), dtype=np.float64)
        values_valid = np.zeros(len(keys_f64), dtype=bool)
        has_values = False
    else:
        values_f64 = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        values_valid = np.ascontiguousarray(_values_valid(values))
        has_values = True

    kernel = getattr(groupby_ext, "groupby_hash_f64_f64", None)
    if kernel is None:
        return None
    groups, row_counts, value_counts, sums, mins, maxs, has_value = kernel(
        keys_f64, values_f64, valid, values_valid, has_values, dropna
    )
    groups = groups.astype(key_dtype, copy=False)
    if op == "size":
        result = row_counts
    elif op == "count":
        result = value_counts
    elif op == "sum":
        result = sums.copy()
        result[~has_value] = np.nan
    elif op == "mean":
        result = np.full(len(groups), np.nan, dtype=np.float64)
        ok = value_counts > 0
        result[ok] = sums[ok] / value_counts[ok]
    elif op == "min":
        result = mins.copy()
        result[~has_value] = np.nan
    elif op == "max":
        result = maxs.copy()
        result[~has_value] = np.nan
    else:  # pragma: no cover
        return None
    return _maybe_sort(groups, result, sort)


def _group_reduce_numpy(  # noqa: C901
    keys: np.ndarray, values: np.ndarray | None, op: str, *, sort: bool, dropna: bool
):
    acc: dict[object, list] = {}
    display: dict[object, object] = {}
    for i, key in enumerate(keys):
        key_item = _python_scalar(key)
        if isinstance(key_item, float) and math.isnan(key_item):
            if dropna:
                continue
            norm_key = _NAN_KEY
        else:
            norm_key = key_item
        display.setdefault(norm_key, key_item)
        state = acc.setdefault(norm_key, [0, 0, 0.0, None, None])
        state[0] += 1
        if values is None:
            continue
        value = _python_scalar(values[i])
        if isinstance(value, float) and math.isnan(value):
            continue
        state[1] += 1
        if op in {"sum", "mean"}:
            state[2] += value
        elif op == "min" and (state[3] is None or value < state[3]):
            state[3] = value
        elif op == "max" and (state[4] is None or value > state[4]):
            state[4] = value

    order = list(acc)
    if sort:
        order.sort(key=lambda k: _group_reduce_sort_key(display[k]))
    groups = np.asarray([display[k] for k in order], dtype=keys.dtype)
    result = []
    for k in order:
        rows, count, total, min_value, max_value = acc[k]
        if op == "size":
            result.append(rows)
        elif op == "count":
            result.append(count)
        elif op == "sum":
            result.append(total if count else _null_value_for(values))
        elif op == "mean":
            result.append(math.nan if count == 0 else total / count)
        elif op == "min":
            result.append(min_value if count else _null_value_for(values))
        elif op == "max":
            result.append(max_value if count else _null_value_for(values))
    return groups, np.asarray(result, dtype=_result_dtype(values, op))


def _group_reduce_sort_key(value: Any) -> tuple[int, Any]:
    """Sort group_reduce keys with None first and NaN groups last."""
    if value is None:
        return (0, "")
    if isinstance(value, float) and math.isnan(value):
        return (2, "")
    return (1, value)


def _resolve_sort(sort: bool | None, *, cheap: bool) -> bool:
    """Resolve a tri-state ``sort`` request for a given execution path.

    ``sort`` is ``True`` (always sort), ``False`` (never sort), or ``None``
    (auto: sort only when this path can do so cheaply).  ``cheap`` is declared
    by the call site: ``True`` for vectorized/free ordering (``np.nonzero``
    ascending or a vectorized lexsort/argsort), ``False`` for a Python
    ``list.sort`` over all groups.
    """
    if sort is None:
        return cheap
    return sort


def _maybe_sort(groups: np.ndarray, result: np.ndarray, sort: bool):
    if sort and len(groups):
        order = np.argsort(groups, kind="stable")
        return groups[order], result[order]
    return groups, result


def _values_valid(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind == "f":
        return ~np.isnan(values)
    return np.ones(len(values), dtype=bool)


def _result_dtype(values: np.ndarray | None, op: str):
    if op in {"size", "count"}:
        return np.int64
    if op == "mean" or values is None:
        return np.float64
    dtype = np.dtype(values.dtype)
    if op == "sum" and dtype.kind in "biu":
        return np.int64
    return dtype


def _null_value_for(values: np.ndarray | None):
    if values is not None and np.dtype(values.dtype).kind in "iu":
        return 0
    return math.nan
