#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Storage backends for CTable.

Two concrete backends:

* :class:`InMemoryTableStorage` — all arrays live in RAM (default when
  ``urlpath`` is not provided).
* :class:`FileTableStorage` — arrays are stored inside a :class:`blosc2.TreeStore`
  rooted at ``urlpath``; logical object metadata lives in ``/_meta`` and table
  data lives under ``/_valid_rows`` and ``/_cols/<name>``.
"""

from __future__ import annotations

import contextlib
import copy
import json
import os
from typing import TYPE_CHECKING, Any

import numpy as np

import blosc2
from blosc2.batch_array import BatchArray
from blosc2.dictionary_column import DictionaryColumn
from blosc2.list_array import ListArray
from blosc2.scalar_array import (
    _make_persistent_backend,
    _open_persistent_backend,
    _ScalarVarLenArray,
    _validate_role_metadata,
)
from blosc2.schema import Utf8Spec
from blosc2.schunk import process_opened_object
from blosc2.utf8_array import Utf8Array, _new_backend_arrays

if TYPE_CHECKING:
    from blosc2.schema import ListSpec

# Directory inside the table root that holds per-column index sidecar files.
_INDEXES_DIR = "_indexes"


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------


class TableStorage:
    """Interface that CTable uses to create/open its backing arrays."""

    def create_column(
        self,
        name: str,
        *,
        dtype: np.dtype,
        shape: tuple[int, ...],
        chunks: tuple[int, ...],
        blocks: tuple[int, ...],
        cparams: dict[str, Any] | None,
        dparams: dict[str, Any] | None,
    ) -> blosc2.NDArray:
        raise NotImplementedError

    def install_column(self, name: str, ndarray: blosc2.NDArray) -> blosc2.NDArray:
        """Store a pre-built NDArray as column *name*, preserving its storage config.

        Faster than create_column + fill when the caller already has the fully
        compressed array (e.g. from NDArray.copy with new block settings).
        Subclasses should override with a path that avoids double recompression.
        """
        raise NotImplementedError

    def open_column(self, name: str) -> blosc2.NDArray:
        raise NotImplementedError

    def create_list_column(
        self,
        name: str,
        *,
        spec: ListSpec,
        cparams: dict[str, Any] | None,
        dparams: dict[str, Any] | None,
    ) -> ListArray:
        raise NotImplementedError

    def install_list_column(self, name: str, list_array: ListArray) -> ListArray:
        """Store a pre-built ListArray as column *name*.

        Faster than create_list_column + extend when the caller already has the
        fully populated array (e.g. from ListArray.copy with chunk_copy).
        Subclasses should override with a path that skips element-wise iteration.
        """
        raise NotImplementedError

    def open_list_column(self, name: str) -> ListArray:
        raise NotImplementedError

    def create_varlen_scalar_column(
        self,
        name: str,
        *,
        spec,
        cparams=None,
        dparams=None,
    ) -> _ScalarVarLenArray:
        raise NotImplementedError

    def open_varlen_scalar_column(self, name: str, spec) -> _ScalarVarLenArray:
        raise NotImplementedError

    def create_dictionary_column(
        self,
        name: str,
        *,
        spec,
        cparams: dict[str, Any] | None = None,
        dparams: dict[str, Any] | None = None,
    ) -> DictionaryColumn:
        raise NotImplementedError

    def open_dictionary_column(self, name: str, spec) -> DictionaryColumn:
        raise NotImplementedError

    def create_valid_rows(
        self,
        *,
        shape: tuple[int, ...],
        chunks: tuple[int, ...],
        blocks: tuple[int, ...],
    ) -> blosc2.NDArray:
        raise NotImplementedError

    def open_valid_rows(self) -> blosc2.NDArray:
        raise NotImplementedError

    def save_schema(self, schema_dict: dict[str, Any]) -> None:
        raise NotImplementedError

    def load_schema(self) -> dict[str, Any] | None:
        raise NotImplementedError

    def table_exists(self) -> bool:
        raise NotImplementedError

    def is_read_only(self) -> bool:
        raise NotImplementedError

    def open_mode(self) -> str | None:
        raise NotImplementedError

    def delete_column(self, name: str) -> None:
        raise NotImplementedError

    def rename_column(self, old: str, new: str) -> blosc2.NDArray:
        raise NotImplementedError

    def close(self) -> None:
        raise NotImplementedError

    def discard(self) -> None:
        """Clean up resources without persisting changes back to the archive."""
        self.close()

    # -- Index catalog and epoch helpers -------------------------------------

    def load_index_catalog(self) -> dict:
        """Return the current index catalog (column_name → descriptor dict)."""
        raise NotImplementedError

    def save_index_catalog(self, catalog: dict) -> None:
        """Persist *catalog* (column_name → descriptor dict)."""
        raise NotImplementedError

    def index_catalog_revision(self) -> int:
        """Return a process-local revision for cache invalidation."""
        return int(getattr(self, "_index_catalog_revision", 0))

    def _bump_index_catalog_revision(self) -> None:
        self._index_catalog_revision = self.index_catalog_revision() + 1

    def get_epoch_counters(self) -> tuple[int, int]:
        """Return ``(value_epoch, visibility_epoch)``."""
        raise NotImplementedError

    def bump_value_epoch(self) -> int:
        """Increment and return the value epoch (data values changed)."""
        raise NotImplementedError

    def bump_visibility_epoch(self) -> int:
        """Increment and return the visibility epoch (row set changed by delete)."""
        raise NotImplementedError

    def index_anchor_path(self, col_name: str) -> str | None:
        """Return the urlpath used as the anchor for index sidecar naming.

        Returns *None* for in-memory storage.  For file-backed storage returns
        a path of the form ``<root>/_indexes/<col_name>/_anchor``.
        """
        raise NotImplementedError


# ---------------------------------------------------------------------------
# In-memory backend
# ---------------------------------------------------------------------------


class InMemoryTableStorage(TableStorage):
    """All arrays are plain in-memory blosc2.NDArray objects."""

    def __init__(self) -> None:
        self._index_catalog: dict = {}
        self._value_epoch: int = 0
        self._visibility_epoch: int = 0

    def create_column(self, name, *, dtype, shape, chunks, blocks, cparams, dparams):
        kwargs: dict[str, Any] = {"chunks": chunks, "blocks": blocks}
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        return blosc2.zeros(shape, dtype=dtype, **kwargs)

    def install_column(self, name, ndarray: blosc2.NDArray) -> blosc2.NDArray:
        """Store a pre-built NDArray as column *name* (skips the zeros+fill pattern)."""
        return ndarray

    def open_column(self, name):
        raise RuntimeError("In-memory tables have no on-disk representation to open.")

    def create_list_column(self, name, *, spec, cparams, dparams):
        kwargs = {}
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        return ListArray(spec=spec, **kwargs)

    def install_list_column(self, name, list_array: ListArray) -> ListArray:
        """Store a pre-built ListArray (in-memory: chunk_copy without urlpath)."""
        return list_array.copy()

    def open_list_column(self, name):
        raise RuntimeError("In-memory tables have no on-disk representation to open.")

    def create_varlen_scalar_column(self, name, *, spec, cparams=None, dparams=None):
        if isinstance(spec, Utf8Spec):
            offsets, data = _new_backend_arrays(cparams, dparams)
            return Utf8Array(spec, offsets, data)
        return _ScalarVarLenArray(spec)

    def open_varlen_scalar_column(self, name, spec):
        raise RuntimeError("In-memory tables have no on-disk representation to open.")

    def create_dictionary_column(
        self,
        name,
        *,
        spec,
        cparams=None,
        dparams=None,
        codes_shape=(4096,),
        codes_chunks=(4096,),
        codes_blocks=(256,),
    ):
        from blosc2.schema import VLStringSpec

        codes = blosc2.zeros(codes_shape, dtype=np.int32, chunks=codes_chunks, blocks=codes_blocks)
        dict_store = _ScalarVarLenArray(VLStringSpec(nullable=False))
        return DictionaryColumn(spec, codes, dict_store)

    def open_dictionary_column(self, name, spec):
        raise RuntimeError("In-memory tables have no on-disk representation to open.")

    def create_valid_rows(self, *, shape, chunks, blocks):
        return blosc2.zeros(shape, dtype=np.bool_, chunks=chunks, blocks=blocks)

    def open_valid_rows(self):
        raise RuntimeError("In-memory tables have no on-disk representation to open.")

    def save_schema(self, schema_dict):
        pass  # nothing to persist

    def load_schema(self):
        return None

    def table_exists(self):
        return False

    def is_read_only(self):
        return False

    def open_mode(self) -> str | None:
        return None

    def delete_column(self, name):
        raise RuntimeError("In-memory tables have no on-disk representation to mutate.")

    def rename_column(self, old: str, new: str):
        raise RuntimeError("In-memory tables have no on-disk representation to mutate.")

    def close(self):
        pass

    # -- Index catalog and epoch helpers -------------------------------------

    def load_index_catalog(self) -> dict:
        return copy.deepcopy(self._index_catalog)

    def save_index_catalog(self, catalog: dict) -> None:
        self._index_catalog = copy.deepcopy(catalog)
        self._bump_index_catalog_revision()

    def get_epoch_counters(self) -> tuple[int, int]:
        return self._value_epoch, self._visibility_epoch

    def bump_value_epoch(self) -> int:
        self._value_epoch += 1
        return self._value_epoch

    def bump_visibility_epoch(self) -> int:
        self._visibility_epoch += 1
        return self._visibility_epoch

    def index_anchor_path(self, col_name: str) -> str | None:
        return None


# ---------------------------------------------------------------------------
# File-backed backend
# ---------------------------------------------------------------------------

_META_KEY = "/_meta"
_VALID_ROWS_KEY = "/_valid_rows"
_COLS_DIR = "_cols"
_VLMETA_KEY = "/_vlmeta"


def split_field_path(path: str) -> tuple[str, ...]:
    """Split a dotted logical field path into segments.

    A backslash escapes separator characters, so ``"a\\.b.c"`` means the
    two-segment path ``("a.b", "c")``.  The empty string is the canonical root.
    """
    if path == "":
        return ()
    parts: list[str] = []
    buf: list[str] = []
    escaped = False
    for ch in path:
        if escaped:
            buf.append(ch)
            escaped = False
        elif ch == "\\":
            escaped = True
        elif ch == ".":
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if escaped:
        buf.append("\\")
    parts.append("".join(buf))
    return tuple(parts)


def join_field_path(parts: tuple[str, ...] | list[str]) -> str:
    """Join logical path segments using dot syntax with backslash escaping."""
    escaped_parts = []
    for part in parts:
        buf: list[str] = []
        for ch in part:
            if ch in {"\\", ".", "/"}:
                buf.append("\\")
            buf.append(ch)
        escaped_parts.append("".join(buf))
    return ".".join(escaped_parts)


def _encode_storage_segment(segment: str) -> str:
    """Percent-encode characters that are structural in logical/storage paths."""
    return segment.replace("%", "%25").replace("/", "%2F").replace(".", "%2E").replace("\\", "%5C")


def _column_name_to_relpath(name: str) -> str:
    """Map a logical column name to a hierarchical path under ``_cols``.

    Unescaped dots are interpreted as nested path separators
    (``a.b.c`` -> ``a/b/c``). Literal dots/slashes/backslashes in field names
    can be represented with :func:`join_field_path` and are percent-encoded in
    the physical storage path.
    """
    return "/".join(_encode_storage_segment(part) for part in split_field_path(name))


# Suffix used to store a dictionary column's value store alongside its codes.
_DICT_SUFFIX = "_dict"

# Key suffix for the byte-blob companion of a utf8 column (its offsets array
# lives at the plain column key).  The literal dot cannot collide with a user
# column: unescaped dots in logical names become path separators and escaped
# dots are percent-encoded, so no column name maps to a key containing ".".
_UTF8_DATA_SUFFIX = ".utf8"


class EmbedStoreTableStorage(TableStorage):
    """Read-only :class:`CTable` storage backed by an in-memory :class:`blosc2.EmbedStore`.

    This is the reconstruction side of :meth:`blosc2.CTable.to_cframe` /
    :func:`blosc2.ctable_from_cframe`.  A cframe deserializes into an
    :class:`blosc2.EmbedStore` (a dict of blosc2 objects keyed by ``/``-paths);
    this storage exposes those entries through the :class:`TableStorage`
    interface so that :meth:`CTable._open_from_storage` can build a normal
    in-memory :class:`CTable` from them.

    The layout mirrors :class:`FileTableStorage`::

        /_meta            SChunk  (vlmeta: kind, version, schema)
        /_valid_rows      NDArray (bool)
        /_vlmeta          SChunk   (user vlmeta, optional)
        /_cols/<name>     NDArray | ListArray | BatchArray
        /_cols/<name>_dict  BatchArray (dictionary column value store)

    The store is read-only; mutation methods raise.
    """

    def __init__(self, estore: blosc2.EmbedStore) -> None:
        self._estore = estore
        self._meta: blosc2.SChunk | None = None
        self._vlmeta: blosc2.SChunk | None = None

    # -- key helpers ------------------------------------------------------

    @staticmethod
    def _col_key(name: str) -> str:
        return f"/{_COLS_DIR}/{_column_name_to_relpath(name)}"

    # -- meta / schema ----------------------------------------------------

    def _open_meta(self) -> blosc2.SChunk:
        if self._meta is None:
            opened = self._estore[_META_KEY]
            if not isinstance(opened, blosc2.SChunk):
                raise ValueError("CTable manifest '/_meta' must be an SChunk.")
            self._meta = opened
        return self._meta

    def check_kind(self) -> None:
        kind = self._open_meta().vlmeta["kind"]
        if isinstance(kind, bytes):
            kind = kind.decode()
        if kind != "ctable":
            raise ValueError(f"Object is not a CTable (kind={kind!r})")

    def load_schema(self) -> dict[str, Any]:
        raw = self._open_meta().vlmeta["schema"]
        if isinstance(raw, bytes):
            raw = raw.decode()
        return json.loads(raw)

    def _open_vlmeta(self) -> blosc2.SChunk | None:
        if self._vlmeta is not None:
            return self._vlmeta
        if _VLMETA_KEY in self._estore:
            opened = self._estore[_VLMETA_KEY]
            if isinstance(opened, blosc2.SChunk):
                self._vlmeta = opened
                return opened
        return None

    # -- columns ----------------------------------------------------------

    def open_column(self, name: str) -> blosc2.NDArray:
        return self._estore[self._col_key(name)]

    def open_list_column(self, name: str) -> ListArray:
        return self._estore[self._col_key(name)]

    def open_varlen_scalar_column(self, name: str, spec) -> _ScalarVarLenArray:
        if isinstance(spec, Utf8Spec):
            offsets = self._estore[self._col_key(name)]
            data = self._estore[self._col_key(name) + _UTF8_DATA_SUFFIX]
            return Utf8Array(spec, offsets, data)
        backend = self._estore[self._col_key(name)]
        return _ScalarVarLenArray(spec, backend)

    def open_dictionary_column(self, name: str, spec) -> DictionaryColumn:
        from blosc2.schema import VLStringSpec

        codes = self._estore[self._col_key(name)]
        dict_backend = self._estore[self._col_key(name) + _DICT_SUFFIX]
        dict_spec = VLStringSpec(nullable=False)
        dict_store = _ScalarVarLenArray(dict_spec, dict_backend)
        return DictionaryColumn(spec, codes, dict_store)

    # -- valid rows -------------------------------------------------------

    def open_valid_rows(self) -> blosc2.NDArray:
        return self._estore[_VALID_ROWS_KEY]

    # -- status -----------------------------------------------------------

    def table_exists(self) -> bool:
        return True

    def is_read_only(self) -> bool:
        return True

    def open_mode(self) -> str | None:
        return "r"

    def close(self) -> None:
        pass

    # -- mutation: not supported (transport is read-only) -----------------

    def _not_supported(self, *args, **kwargs):
        raise RuntimeError("EmbedStoreTableStorage is read-only (cframe transport).")

    create_column = _not_supported
    install_column = _not_supported
    create_list_column = _not_supported
    install_list_column = _not_supported
    create_varlen_scalar_column = _not_supported
    create_dictionary_column = _not_supported
    create_valid_rows = _not_supported
    save_schema = _not_supported
    save_vlmeta = _not_supported
    delete_column = _not_supported
    rename_column = _not_supported

    # -- index catalog (in-memory defaults; indexes are not transported) ---

    def load_index_catalog(self) -> dict:
        return {}

    def save_index_catalog(self, catalog: dict) -> None:
        pass

    def get_epoch_counters(self) -> tuple[int, int]:
        return 0, 0

    def bump_value_epoch(self) -> int:
        return 0

    def bump_visibility_epoch(self) -> int:
        return 0

    def index_anchor_path(self, col_name: str) -> str | None:
        return None


class FileTableStorage(TableStorage):
    """Arrays stored as TreeStore leaves inside *urlpath*.

    Parameters
    ----------
    urlpath:
        Path to the backing TreeStore (typically ``.b2d`` or ``.b2z``).
    mode:
        ``'w'`` — create (overwrite existing files).
        ``'a'`` — open existing or create new.
        ``'r'`` — open existing read-only.
    mmap_mode:
        ``'r'`` to memory-map the backing store (requires ``mode='r'``);
        members are then read from mapped pages — for ``.b2z``, in place at
        their offsets inside the single mapped container file.
    """

    def __init__(
        self,
        urlpath: str,
        mode: str,
        store: blosc2.TreeStore | None = None,
        mmap_mode: str | None = None,
    ) -> None:
        if mode not in ("r", "a", "w"):
            raise ValueError(f"mode must be 'r', 'a', or 'w'; got {mode!r}")
        if mmap_mode is not None and mode != "r":
            raise ValueError("mmap_mode requires mode='r'")
        self._root = urlpath
        self._mode = mode
        self._mmap_mode = mmap_mode
        self._meta: blosc2.SChunk | None = None
        self._vlmeta: blosc2.SChunk | None = None
        # CTable internals must always use external-file storage (never the
        # embed store) so that small SChunk overwrites (e.g. _meta with
        # nbytes=0) are reliably persisted.  Normalise a pre-existing store
        # that was opened by generic dispatch without this setting.
        if store is not None and store.threshold != 0:
            store.threshold = 0
        self._store: blosc2.TreeStore | None = store
        self._registered_sidecar_paths: list[str] = []

    # ------------------------------------------------------------------
    # Key helpers
    # ------------------------------------------------------------------

    @property
    def _meta_path(self) -> str:
        return self._key_to_path(_META_KEY)

    @property
    def _valid_rows_path(self) -> str:
        return self._key_to_path(_VALID_ROWS_KEY)

    @property
    def _vlmeta_path(self) -> str:
        return self._key_to_path(_VLMETA_KEY)

    def _col_path(self, name: str) -> str:
        return self._key_to_path(self._col_key(name))

    def _list_col_path(self, name: str) -> str:
        rel_key = self._col_key(name).lstrip("/")
        # Use working_dir so .b2z stores write into their temp dir (gets zipped on close).
        # For .b2d, working_dir == self._root, so behaviour is unchanged.
        return os.path.join(self._open_store().working_dir, rel_key + ".b2b")

    def _dict_col_path(self, name: str) -> str:
        """Path for the dictionary values store of a dictionary column."""
        rel_key = self._col_key(name).lstrip("/")
        return os.path.join(self._open_store().working_dir, rel_key + "_dict.b2b")

    def _col_key(self, name: str) -> str:
        return f"/{_COLS_DIR}/{_column_name_to_relpath(name)}"

    def _key_to_path(self, key: str) -> str:
        rel_key = key.lstrip("/")
        suffix = ".b2f" if key in (_META_KEY, _VLMETA_KEY) else ".b2nd"
        if self._root.endswith(".b2d"):
            return os.path.join(self._root, rel_key + suffix)
        return os.path.join(self._root, rel_key + suffix)

    def _open_store(self) -> blosc2.TreeStore:
        if self._store is None:
            kwargs: dict[str, Any] = {"mode": self._mode}
            if self._mode != "r":
                # Force table internals to be stored as proper external leaves so
                # reopened arrays stay live and mutable through the TreeStore.
                kwargs["threshold"] = 0
            elif self._mmap_mode is not None:
                # Memory-map the read-only store: members are read straight
                # from mapped pages (for .b2z, at their offsets inside the
                # single mapped container file, shared across readers).
                kwargs["mmap_mode"] = self._mmap_mode
            self._store = blosc2.TreeStore(self._root, **kwargs)
        return self._store

    # ------------------------------------------------------------------
    # TableStorage interface
    # ------------------------------------------------------------------

    def table_exists(self) -> bool:
        return os.path.exists(self._root)

    def is_read_only(self) -> bool:
        return self._mode == "r"

    def open_mode(self) -> str | None:
        return self._mode

    def create_column(self, name, *, dtype, shape, chunks, blocks, cparams, dparams):
        kwargs: dict[str, Any] = {
            "chunks": chunks,
            "blocks": blocks,
        }
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        col = blosc2.zeros(shape, dtype=dtype, **kwargs)
        store = self._open_store()
        store[self._col_key(name)] = col
        return store[self._col_key(name)]

    def install_column(self, name, ndarray: blosc2.NDArray) -> blosc2.NDArray:
        """Store a pre-built NDArray as column *name* (skips the zeros+fill pattern)."""
        store = self._open_store()
        store[self._col_key(name)] = ndarray
        return store[self._col_key(name)]

    def open_column(self, name: str) -> blosc2.NDArray:
        return self._open_store()[self._col_key(name)]

    def create_list_column(self, name, *, spec, cparams, dparams):
        kwargs: dict[str, Any] = {"urlpath": self._list_col_path(name), "mode": "w", "contiguous": True}
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        os.makedirs(os.path.dirname(self._list_col_path(name)), exist_ok=True)
        return ListArray(spec=spec, **kwargs)

    def install_list_column(self, name, list_array: ListArray) -> ListArray:
        """Bulk-copy a pre-built ListArray to the column path (chunk-level transfer)."""
        dest_path = self._list_col_path(name)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        result = list_array.copy(urlpath=dest_path, mode="w", contiguous=True)
        result.flush()
        return result

    def open_list_column(self, name: str) -> ListArray:
        store = self._open_store()
        if store.is_zip_store and self._mode == "r":
            # In read mode, .b2z is never extracted — read the member at its zip offset directly.
            rel = self._col_key(name).lstrip("/") + ".b2b"
            if rel not in store.offsets:
                raise KeyError(f"List column {name!r} not found in {self._root!r}")
            opened = blosc2.blosc2_ext.open(
                store.b2z_path, mode="r", offset=store.offsets[rel]["offset"], mmap_mode=store.mmap_mode
            )
            return process_opened_object(opened)
        return blosc2.open(self._list_col_path(name), mode=self._mode)

    def create_varlen_scalar_column(self, name, *, spec, cparams=None, dparams=None) -> _ScalarVarLenArray:
        if isinstance(spec, Utf8Spec):
            offsets, data = _new_backend_arrays(cparams, dparams)
            store = self._open_store()
            key = self._col_key(name)
            data_key = key + _UTF8_DATA_SUFFIX
            store[key] = offsets
            store[data_key] = data
            return Utf8Array(spec, store[key], store[data_key])
        urlpath = self._list_col_path(name)
        backend = _make_persistent_backend(spec, urlpath, "w", cparams=cparams, dparams=dparams)
        return _ScalarVarLenArray(spec, backend)

    def open_varlen_scalar_column(self, name: str, spec) -> _ScalarVarLenArray:
        if isinstance(spec, Utf8Spec):
            store = self._open_store()
            offsets = store[self._col_key(name)]
            data = store[self._col_key(name) + _UTF8_DATA_SUFFIX]
            return Utf8Array(spec, offsets, data)
        store = self._open_store()
        path = self._list_col_path(name)
        if store.is_zip_store and self._mode == "r":
            rel = self._col_key(name).lstrip("/") + ".b2b"
            if rel not in store.offsets:
                raise KeyError(f"Varlen scalar column {name!r} not found in {self._root!r}")
            backend = BatchArray(
                _from_schunk=blosc2.blosc2_ext.open(
                    store.b2z_path, mode="r", offset=store.offsets[rel]["offset"], mmap_mode=store.mmap_mode
                )
            )
        else:
            backend = _open_persistent_backend(path, self._mode, spec=spec)
        _validate_role_metadata(backend, spec)
        return _ScalarVarLenArray(spec, backend)

    def create_dictionary_column(
        self,
        name,
        *,
        spec,
        cparams=None,
        dparams=None,
        codes_shape=(4096,),
        codes_chunks=(4096,),
        codes_blocks=(256,),
    ) -> DictionaryColumn:
        from blosc2.schema import VLStringSpec

        # Codes: stored as a regular NDArray under _cols/name
        codes = self.create_column(
            name,
            dtype=np.int32,
            shape=codes_shape,
            chunks=codes_chunks,
            blocks=codes_blocks,
            cparams=cparams,
            dparams=dparams,
        )
        # Dictionary values: stored as a varlen scalar (vlstring) at name_dict.b2b
        dict_spec = VLStringSpec(nullable=False)
        dict_path = self._dict_col_path(name)
        dict_backend = _make_persistent_backend(dict_spec, dict_path, "w")
        dict_store = _ScalarVarLenArray(dict_spec, dict_backend)
        return DictionaryColumn(spec, codes, dict_store)

    def open_dictionary_column(self, name: str, spec) -> DictionaryColumn:
        from blosc2.schema import VLStringSpec

        codes = self.open_column(name)
        dict_spec = VLStringSpec(nullable=False)
        store = self._open_store()
        dict_path = self._dict_col_path(name)
        if store.is_zip_store and self._mode == "r":
            rel = self._col_key(name).lstrip("/") + "_dict.b2b"
            if rel not in store.offsets:
                raise KeyError(f"Dictionary column dict store {name!r} not found in {self._root!r}")
            dict_backend = BatchArray(
                _from_schunk=blosc2.blosc2_ext.open(
                    store.b2z_path, mode="r", offset=store.offsets[rel]["offset"], mmap_mode=store.mmap_mode
                )
            )
        else:
            dict_backend = _open_persistent_backend(dict_path, self._mode, spec=dict_spec)
        dict_store = _ScalarVarLenArray(dict_spec, dict_backend)
        return DictionaryColumn(spec, codes, dict_store)

    def create_valid_rows(self, *, shape, chunks, blocks):
        valid_rows = blosc2.zeros(
            shape,
            dtype=np.bool_,
            chunks=chunks,
            blocks=blocks,
        )
        store = self._open_store()
        store[_VALID_ROWS_KEY] = valid_rows
        return store[_VALID_ROWS_KEY]

    def open_valid_rows(self) -> blosc2.NDArray:
        return self._open_store()[_VALID_ROWS_KEY]

    def save_schema(self, schema_dict: dict[str, Any]) -> None:
        """Write *schema_dict* (plus kind/version markers) to ``/_meta``."""
        meta = blosc2.SChunk()
        meta.vlmeta["kind"] = "ctable"
        meta.vlmeta["version"] = 1
        meta.vlmeta["schema"] = json.dumps(schema_dict)
        store = self._open_store()
        store[_META_KEY] = meta
        opened = store[_META_KEY]
        if not isinstance(opened, blosc2.SChunk):
            raise ValueError("CTable manifest '/_meta' must materialize as an SChunk.")
        self._meta = opened

    def save_vlmeta(self, schunk: blosc2.SChunk) -> None:
        """Persist the user vlmeta SChunk to the storage."""
        if self._mode == "r":
            return
        self._vlmeta = schunk
        if self._store is not None:
            self._store[_VLMETA_KEY] = schunk

    def _open_vlmeta(self) -> blosc2.SChunk | None:
        """Open (or return cached) the ``/_vlmeta`` SChunk.

        Returns ``None`` if the file does not exist (read-only open of a
        table that never had user vlmeta written).
        """
        uv = getattr(self, "_vlmeta", None)
        if uv is not None:
            return uv
        # Try TreeStore first
        try:
            opened = self._open_store()[_VLMETA_KEY]
            if isinstance(opened, blosc2.SChunk):
                self._vlmeta = opened
                return opened
        except (KeyError, FileNotFoundError):
            pass
        # Fallback: try opening the filesystem path directly
        uv_path = self._vlmeta_path
        if os.path.exists(uv_path):
            opened = blosc2.open(uv_path, mode="r")
            if isinstance(opened, blosc2.SChunk):
                self._vlmeta = opened
                return opened
        return None

    def _open_meta(self) -> blosc2.SChunk:
        """Open (or return cached) the ``/_meta`` SChunk."""
        if self._meta is None:
            try:
                opened = self._open_store()[_META_KEY]
            except KeyError as exc:
                raise FileNotFoundError(f"No CTable manifest found at {self._root!r}") from exc
            if not isinstance(opened, blosc2.SChunk):
                raise ValueError(f"CTable manifest at {self._root!r} must be an SChunk.")
            self._meta = opened
        return self._meta

    def load_schema(self) -> dict[str, Any]:
        """Read and return the schema dict stored in ``/_meta``."""
        raw = self._open_meta().vlmeta["schema"]
        if isinstance(raw, bytes):
            raw = raw.decode()
        return json.loads(raw)

    def check_kind(self) -> None:
        """Raise :exc:`ValueError` if ``_meta`` does not identify a CTable."""
        kind = self._open_meta().vlmeta["kind"]
        if isinstance(kind, bytes):
            kind = kind.decode()
        if kind != "ctable":
            raise ValueError(f"Path {self._root!r} does not contain a CTable (kind={kind!r}).")

    def column_names_from_schema(self) -> list[str]:
        d = self.load_schema()
        return [c["name"] for c in d["columns"]]

    def delete_column(self, name: str) -> None:
        key = self._col_key(name)
        if key in self._open_store():
            del self._open_store()[key]
            data_key = key + _UTF8_DATA_SUFFIX
            if data_key in self._open_store():
                del self._open_store()[data_key]
            return
        list_path = self._list_col_path(name)
        if os.path.exists(list_path):
            blosc2.remove_urlpath(list_path)
            return
        raise KeyError(name)

    def rename_column(self, old: str, new: str):
        store = self._open_store()
        old_key = self._col_key(old)
        new_key = self._col_key(new)
        if old_key in store:
            store[new_key] = store[old_key]
            del store[old_key]
            old_data_key = old_key + _UTF8_DATA_SUFFIX
            if old_data_key in store:
                store[new_key + _UTF8_DATA_SUFFIX] = store[old_data_key]
                del store[old_data_key]
            return store[new_key]
        old_path = self._list_col_path(old)
        new_path = self._list_col_path(new)
        if os.path.exists(old_path):
            os.makedirs(os.path.dirname(new_path), exist_ok=True)
            os.replace(old_path, new_path)
            return blosc2.open(new_path, mode=self._mode)
        raise KeyError(old)

    def close(self) -> None:
        self._unregister_sidecar_zip_paths()
        self._evict_cached_index_handles()
        if self._store is not None:
            self._store.close()
            self._store = None
        self._meta = None

    def discard(self) -> None:
        """Clean up without repacking the .b2z archive."""
        self._unregister_sidecar_zip_paths()
        self._evict_cached_index_handles()
        if self._store is not None:
            self._store.discard()
            self._store = None
        self._meta = None

    def _evict_cached_index_handles(self) -> None:
        """Release process-global index handle/data caches for this table's
        files, so closing it does not leak a file descriptor per table."""
        from blosc2.indexing import evict_cached_index_handles

        with contextlib.suppress(Exception):
            evict_cached_index_handles(self._root)

    def _unregister_sidecar_zip_paths(self) -> None:
        if not self._registered_sidecar_paths:
            return
        from blosc2.indexing import _SIDECAR_ZIP_REGISTRY

        for path in self._registered_sidecar_paths:
            _SIDECAR_ZIP_REGISTRY.pop(path, None)
        self._registered_sidecar_paths.clear()

    # -- Index catalog and epoch helpers -------------------------------------

    @staticmethod
    def _walk_descriptor_paths(descriptor: dict):
        """Yield (obj, key) for every string value that looks like a file path."""
        stack = [descriptor]
        while stack:
            obj = stack.pop()
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if (k == "path" or k.endswith("_path")) and isinstance(v, str):
                        yield obj, k
                    elif isinstance(v, (dict, list)):
                        stack.append(v)
            elif isinstance(obj, list):
                for item in obj:
                    if isinstance(item, (dict, list)):
                        stack.append(item)

    @staticmethod
    def _relativize_descriptor(descriptor: dict, working_dir: str) -> dict:
        """Replace paths inside *working_dir* with working-dir relative paths."""
        abs_working_dir = os.path.abspath(working_dir)
        prefix = abs_working_dir.rstrip(os.sep) + os.sep
        d = copy.deepcopy(descriptor)
        for obj, key in FileTableStorage._walk_descriptor_paths(d):
            v = obj[key]
            abs_v = v if os.path.isabs(v) else os.path.abspath(v)
            if abs_v.startswith(prefix):
                obj[key] = abs_v[len(prefix) :].replace(os.sep, "/")
        return d

    @staticmethod
    def _absolutize_descriptor(descriptor: dict, working_dir: str) -> dict:
        """Expand working-dir relative paths back to absolute paths."""
        d = copy.deepcopy(descriptor)
        for obj, key in FileTableStorage._walk_descriptor_paths(d):
            v = obj[key]
            if not os.path.isabs(v):
                obj[key] = os.path.abspath(v) if os.path.exists(v) else os.path.join(working_dir, v)
        return d

    def _register_index_zip_paths(self, store, descriptor: dict) -> None:
        """Register sidecar paths from *descriptor* in the zip-offset registry.

        This lets indexing code open sidecar arrays directly at their byte offset
        inside the .b2z archive, avoiding the need to extract them to disk first.
        """
        from blosc2.indexing import _SIDECAR_ZIP_REGISTRY

        working_dir = store.working_dir
        for obj, key in self._walk_descriptor_paths(descriptor):
            abs_path = obj[key]
            if abs_path in _SIDECAR_ZIP_REGISTRY:
                continue
            rel = os.path.relpath(abs_path, working_dir).replace(os.sep, "/")
            info = store.offsets.get(rel)
            if info is None:
                continue
            _SIDECAR_ZIP_REGISTRY[abs_path] = (store.b2z_path, info["offset"])
            self._registered_sidecar_paths.append(abs_path)

    def load_index_catalog(self) -> dict:
        meta = self._open_meta()
        raw = meta.vlmeta.get("index_catalog")
        if not isinstance(raw, dict):
            return {}
        catalog = copy.deepcopy(raw)
        store = self._open_store()
        working_dir = store.working_dir
        # Expand relative paths and, for b2z read mode, register sidecars in the
        # zip-offset registry so indexing code can open them without extraction.
        for col_name, descriptor in catalog.items():
            catalog[col_name] = self._absolutize_descriptor(descriptor, working_dir)
            if store.is_zip_store and self._mode == "r":
                self._register_index_zip_paths(store, catalog[col_name])
        return catalog

    def save_index_catalog(self, catalog: dict) -> None:
        meta = self._open_meta()
        working_dir = self._open_store().working_dir
        relativized = {col: self._relativize_descriptor(desc, working_dir) for col, desc in catalog.items()}
        meta.vlmeta["index_catalog"] = relativized
        self._bump_index_catalog_revision()

    def get_epoch_counters(self) -> tuple[int, int]:
        meta = self._open_meta()
        ve = int(meta.vlmeta.get("value_epoch", 0) or 0)
        vis_e = int(meta.vlmeta.get("visibility_epoch", 0) or 0)
        return ve, vis_e

    def bump_value_epoch(self) -> int:
        meta = self._open_meta()
        ve = int(meta.vlmeta.get("value_epoch", 0) or 0) + 1
        meta.vlmeta["value_epoch"] = ve
        return ve

    def bump_visibility_epoch(self) -> int:
        meta = self._open_meta()
        vis_e = int(meta.vlmeta.get("visibility_epoch", 0) or 0) + 1
        meta.vlmeta["visibility_epoch"] = vis_e
        return vis_e

    def index_anchor_path(self, col_name: str) -> str | None:
        return os.path.join(self._open_store().working_dir, _INDEXES_DIR, col_name, "_anchor")


# ---------------------------------------------------------------------------
# TreeStore-backed backend (inline subtree layout)
# ---------------------------------------------------------------------------


class TreeStoreTableStorage(TableStorage):
    """TableStorage backend that stores a CTable inline inside an outer TreeStore.

    All CTable components are written as normal external leaves under *root_key*
    inside *store*'s working directory.  This avoids nested ZIP files and allows
    the entire bundle to be packed as a flat ``.b2z`` archive.

    Parameters
    ----------
    store:
        The outer :class:`blosc2.TreeStore` (or its subtree view) that will
        hold the CTable internals.
    root_key:
        Full absolute key where the CTable lives, e.g. ``"/table"``.
    mode:
        Open mode (``'r'``, ``'a'``, or ``'w'``).  Should match ``store.mode``.
    owns_store:
        If ``True``, ``close()`` / ``discard()`` will also close / discard
        *store*.  Use ``False`` (default) when *store* is owned by the caller.
    """

    def __init__(
        self,
        store: blosc2.TreeStore,
        root_key: str,
        mode: str,
        owns_store: bool = False,
    ) -> None:
        self._store = store
        self._root_key = root_key.rstrip("/")
        self._mode = mode
        self._owns_store = owns_store
        self._meta: blosc2.SChunk | None = None
        self._vlmeta: blosc2.SChunk | None = None
        self._registered_sidecar_paths: list[str] = []

    # ------------------------------------------------------------------
    # Key / path helpers
    # ------------------------------------------------------------------

    def _table_key(self, logical_key: str) -> str:
        """Translate a CTable-internal logical key to an outer-store absolute key.

        For example, if *root_key* is ``"/table"`` and *logical_key* is
        ``"/_meta"``, the result is ``"/table/_meta"``.
        """
        return self._root_key + logical_key

    def _working_dir(self) -> str:
        return self._store.working_dir

    def _dest_path(self, logical_key: str, ext: str) -> str:
        """Absolute filesystem path for the external leaf file."""
        rel = self._table_key(logical_key).lstrip("/")
        return os.path.join(self._working_dir(), rel + ext)

    def _write_leaf(self, logical_key: str, value: Any, ext: str) -> None:
        """Write *value* as a raw cframe file and register it in the outer
        store's map_tree so DictStore can find it again on open."""
        dest_path = self._dest_path(logical_key, ext)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        if isinstance(value, blosc2.SChunk) or hasattr(value, "to_cframe"):
            with open(dest_path, "wb") as f:
                f.write(value.to_cframe())
        else:
            value.save(urlpath=dest_path, mode="w")
        rel_path = os.path.relpath(dest_path, self._working_dir()).replace(os.sep, "/")
        full_key = self._table_key(logical_key)
        self._store.map_tree[full_key] = rel_path
        self._store._modified = True

    def _open_leaf(self, logical_key: str) -> Any:
        """Open a leaf via the outer store's map_tree / estore logic."""
        from blosc2.dict_store import DictStore

        full_key = self._table_key(logical_key)
        return DictStore.__getitem__(self._store, full_key)

    def _col_logical_key(self, name: str) -> str:
        return f"/{_COLS_DIR}/{_column_name_to_relpath(name)}"

    def _list_col_path(self, name: str) -> str:
        """Filesystem path for a list-style column (``.b2b``)."""
        return self._dest_path(self._col_logical_key(name), ".b2b")

    # ------------------------------------------------------------------
    # TableStorage interface — lifecycle
    # ------------------------------------------------------------------

    def table_exists(self) -> bool:
        full_key = self._table_key("/_meta")
        return full_key in self._store.map_tree or full_key in self._store._estore

    def is_read_only(self) -> bool:
        return self._mode == "r"

    def open_mode(self) -> str | None:
        return self._mode

    def close(self) -> None:
        self._unregister_sidecar_zip_paths()
        self._evict_cached_index_handles()
        if self._owns_store and self._store is not None:
            self._store.close()
            self._store = None
        self._meta = None

    def discard(self) -> None:
        self._unregister_sidecar_zip_paths()
        self._evict_cached_index_handles()
        if self._owns_store and self._store is not None:
            self._store.discard()
            self._store = None
        self._meta = None

    def _evict_cached_index_handles(self) -> None:
        """Release process-global index handle/data caches for this table's
        subtree, so closing it does not leak a file descriptor per table."""
        from blosc2.indexing import evict_cached_index_handles

        with contextlib.suppress(Exception):
            root = os.path.join(self._store.working_dir, self._root_key.lstrip("/"))
            evict_cached_index_handles(root)

    def _unregister_sidecar_zip_paths(self) -> None:
        if not self._registered_sidecar_paths:
            return
        from blosc2.indexing import _SIDECAR_ZIP_REGISTRY

        for path in self._registered_sidecar_paths:
            _SIDECAR_ZIP_REGISTRY.pop(path, None)
        self._registered_sidecar_paths.clear()

    # ------------------------------------------------------------------
    # TableStorage interface — columns and valid_rows
    # ------------------------------------------------------------------

    def create_column(
        self,
        name: str,
        *,
        dtype: np.dtype,
        shape: tuple[int, ...],
        chunks: tuple[int, ...],
        blocks: tuple[int, ...],
        cparams: dict[str, Any] | None,
        dparams: dict[str, Any] | None,
    ) -> blosc2.NDArray:
        kwargs: dict[str, Any] = {"chunks": chunks, "blocks": blocks}
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        dest_path = self._dest_path(self._col_logical_key(name), ".b2nd")
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        col = blosc2.zeros(shape, dtype=dtype, urlpath=dest_path, mode="w", **kwargs)
        rel_path = os.path.relpath(dest_path, self._working_dir()).replace(os.sep, "/")
        self._store.map_tree[self._table_key(self._col_logical_key(name))] = rel_path
        self._store._modified = True
        return col

    def install_column(self, name: str, ndarray: blosc2.NDArray) -> blosc2.NDArray:
        dest_path = self._dest_path(self._col_logical_key(name), ".b2nd")
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        saved = ndarray.copy(urlpath=dest_path)
        rel_path = os.path.relpath(dest_path, self._working_dir()).replace(os.sep, "/")
        self._store.map_tree[self._table_key(self._col_logical_key(name))] = rel_path
        self._store._modified = True
        return saved

    def open_column(self, name: str) -> blosc2.NDArray:
        return self._open_leaf(self._col_logical_key(name))

    def create_list_column(
        self,
        name: str,
        *,
        spec: ListSpec,
        cparams: dict[str, Any] | None,
        dparams: dict[str, Any] | None,
    ) -> ListArray:
        kwargs: dict[str, Any] = {
            "urlpath": self._list_col_path(name),
            "mode": "w",
            "contiguous": True,
        }
        if cparams is not None:
            kwargs["cparams"] = cparams
        if dparams is not None:
            kwargs["dparams"] = dparams
        os.makedirs(os.path.dirname(self._list_col_path(name)), exist_ok=True)
        return ListArray(spec=spec, **kwargs)

    def install_list_column(self, name: str, list_array: ListArray) -> ListArray:
        """Bulk-copy a pre-built ListArray to the column path (chunk-level transfer)."""
        dest_path = self._list_col_path(name)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        result = list_array.copy(urlpath=dest_path, mode="w", contiguous=True)
        result.flush()
        return result

    def open_list_column(self, name: str) -> ListArray:
        if self._store.is_zip_store and self._mode == "r":
            rel = self._table_key(self._col_logical_key(name)).lstrip("/") + ".b2b"
            if rel not in self._store.offsets:
                raise KeyError(f"List column {name!r} not found in {self._store.localpath!r}")
            opened = blosc2.blosc2_ext.open(
                self._store.b2z_path,
                mode="r",
                offset=self._store.offsets[rel]["offset"],
                mmap_mode=self._store.mmap_mode,
            )
            return process_opened_object(opened)
        return blosc2.open(self._list_col_path(name), mode=self._mode)

    def create_varlen_scalar_column(
        self,
        name: str,
        *,
        spec,
        cparams=None,
        dparams=None,
    ) -> _ScalarVarLenArray:
        if isinstance(spec, Utf8Spec):
            logical_key = self._col_logical_key(name)
            offsets_path = self._dest_path(logical_key, ".b2nd")
            data_path = self._dest_path(logical_key + _UTF8_DATA_SUFFIX, ".b2nd")
            os.makedirs(os.path.dirname(offsets_path), exist_ok=True)
            offsets, data = _new_backend_arrays(
                cparams, dparams, offsets_urlpath=offsets_path, data_urlpath=data_path
            )
            for logical, dest_path in (
                (logical_key, offsets_path),
                (logical_key + _UTF8_DATA_SUFFIX, data_path),
            ):
                rel_path = os.path.relpath(dest_path, self._working_dir()).replace(os.sep, "/")
                self._store.map_tree[self._table_key(logical)] = rel_path
            self._store._modified = True
            return Utf8Array(spec, offsets, data)
        urlpath = self._list_col_path(name)
        os.makedirs(os.path.dirname(urlpath), exist_ok=True)
        return _make_persistent_backend(spec, urlpath, "w", cparams=cparams, dparams=dparams)

    def open_varlen_scalar_column(self, name: str, spec) -> _ScalarVarLenArray:
        if isinstance(spec, Utf8Spec):
            logical_key = self._col_logical_key(name)
            offsets = self._open_leaf(logical_key)
            data = self._open_leaf(logical_key + _UTF8_DATA_SUFFIX)
            return Utf8Array(spec, offsets, data)
        if self._store.is_zip_store and self._mode == "r":
            rel = self._table_key(self._col_logical_key(name)).lstrip("/") + ".b2b"
            if rel not in self._store.offsets:
                raise KeyError(f"Varlen scalar column {name!r} not found in {self._store.localpath!r}")
            backend = BatchArray(
                _from_schunk=blosc2.blosc2_ext.open(
                    self._store.b2z_path,
                    mode="r",
                    offset=self._store.offsets[rel]["offset"],
                    mmap_mode=self._store.mmap_mode,
                )
            )
        else:
            backend = _open_persistent_backend(self._list_col_path(name), self._mode, spec=spec)
        _validate_role_metadata(backend, spec)
        return _ScalarVarLenArray(spec, backend)

    def _dict_col_path(self, name: str) -> str:
        """Path for the dictionary values store of a dictionary column."""
        return self._dest_path(self._col_logical_key(name), "_dict.b2b")

    def create_dictionary_column(
        self,
        name: str,
        *,
        spec,
        cparams=None,
        dparams=None,
        codes_shape=(4096,),
        codes_chunks=(4096,),
        codes_blocks=(256,),
    ) -> DictionaryColumn:
        from blosc2.schema import VLStringSpec

        codes = self.create_column(
            name,
            dtype=np.int32,
            shape=codes_shape,
            chunks=codes_chunks,
            blocks=codes_blocks,
            cparams=cparams,
            dparams=dparams,
        )
        dict_spec = VLStringSpec(nullable=False)
        dict_path = self._dict_col_path(name)
        os.makedirs(os.path.dirname(dict_path), exist_ok=True)
        dict_backend = _make_persistent_backend(dict_spec, dict_path, "w")
        dict_store = _ScalarVarLenArray(dict_spec, dict_backend)
        return DictionaryColumn(spec, codes, dict_store)

    def open_dictionary_column(self, name: str, spec) -> DictionaryColumn:
        from blosc2.schema import VLStringSpec

        codes = self.open_column(name)
        dict_spec = VLStringSpec(nullable=False)
        if self._store.is_zip_store and self._mode == "r":
            rel = self._table_key(self._col_logical_key(name)).lstrip("/") + "_dict.b2b"
            if rel not in self._store.offsets:
                raise KeyError(
                    f"Dictionary column dict store {name!r} not found in {self._store.localpath!r}"
                )
            dict_backend = BatchArray(
                _from_schunk=blosc2.blosc2_ext.open(
                    self._store.b2z_path,
                    mode="r",
                    offset=self._store.offsets[rel]["offset"],
                    mmap_mode=self._store.mmap_mode,
                )
            )
        else:
            dict_backend = _open_persistent_backend(self._dict_col_path(name), self._mode, spec=dict_spec)
        dict_store = _ScalarVarLenArray(dict_spec, dict_backend)
        return DictionaryColumn(spec, codes, dict_store)

    def create_valid_rows(
        self,
        *,
        shape: tuple[int, ...],
        chunks: tuple[int, ...],
        blocks: tuple[int, ...],
    ) -> blosc2.NDArray:
        dest_path = self._dest_path("/_valid_rows", ".b2nd")
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        valid_rows = blosc2.zeros(
            shape,
            dtype=np.bool_,
            chunks=chunks,
            blocks=blocks,
            urlpath=dest_path,
            mode="w",
        )
        rel_path = os.path.relpath(dest_path, self._working_dir()).replace(os.sep, "/")
        self._store.map_tree[self._table_key("/_valid_rows")] = rel_path
        self._store._modified = True
        return valid_rows

    def open_valid_rows(self) -> blosc2.NDArray:
        return self._open_leaf("/_valid_rows")

    # ------------------------------------------------------------------
    # TableStorage interface — schema and manifest
    # ------------------------------------------------------------------

    def save_schema(self, schema_dict: dict[str, Any]) -> None:
        meta = blosc2.SChunk()
        meta.vlmeta["kind"] = "ctable"
        meta.vlmeta["version"] = 1
        meta.vlmeta["schema"] = json.dumps(schema_dict)
        self._write_leaf("/_meta", meta, ".b2f")
        opened = self._open_leaf("/_meta")
        if not isinstance(opened, blosc2.SChunk):
            raise ValueError("CTable manifest '/_meta' must materialise as an SChunk.")
        self._meta = opened

    def _open_meta(self) -> blosc2.SChunk:
        if self._meta is None:
            try:
                opened = self._open_leaf("/_meta")
            except KeyError as exc:
                raise FileNotFoundError(f"No CTable manifest found at {self._root_key!r}") from exc
            if not isinstance(opened, blosc2.SChunk):
                raise ValueError(f"CTable manifest at {self._root_key!r} must be an SChunk.")
            self._meta = opened
        return self._meta

    def load_schema(self) -> dict[str, Any]:
        raw = self._open_meta().vlmeta["schema"]
        if isinstance(raw, bytes):
            raw = raw.decode()
        return json.loads(raw)

    def check_kind(self) -> None:
        kind = self._open_meta().vlmeta["kind"]
        if isinstance(kind, bytes):
            kind = kind.decode()
        if kind != "ctable":
            raise ValueError(f"Object at {self._root_key!r} is not a CTable (kind={kind!r})")

    def save_vlmeta(self, schunk: blosc2.SChunk) -> None:
        """Persist the user vlmeta SChunk to the outer TreeStore."""
        if self._mode == "r":
            return
        self._vlmeta = schunk
        self._write_leaf("/_vlmeta", schunk, ".b2f")

    def _open_vlmeta(self) -> blosc2.SChunk | None:
        """Open (or return cached) the ``/_vlmeta`` SChunk.

        Returns ``None`` if the leaf does not exist (read-only open of a
        table that never had user vlmeta written).
        """
        uv = getattr(self, "_vlmeta", None)
        if uv is not None:
            return uv
        try:
            opened = self._open_leaf("/_vlmeta")
        except (KeyError, FileNotFoundError):
            return None
        if not isinstance(opened, blosc2.SChunk):
            return None
        self._vlmeta = opened
        return opened

    def column_names_from_schema(self) -> list[str]:
        return [c["name"] for c in self.load_schema()["columns"]]

    def delete_column(self, name: str) -> None:
        full_key = self._table_key(self._col_logical_key(name))
        if full_key in self._store.map_tree:
            keys = [full_key]
            data_key = self._table_key(self._col_logical_key(name) + _UTF8_DATA_SUFFIX)
            if data_key in self._store.map_tree:
                keys.append(data_key)
            for key in keys:
                filepath = self._store.map_tree.pop(key)
                full_path = os.path.join(self._working_dir(), filepath)
                if os.path.exists(full_path):
                    os.remove(full_path)
            return
        list_path = self._list_col_path(name)
        if os.path.exists(list_path):
            blosc2.remove_urlpath(list_path)
            return
        raise KeyError(name)

    def rename_column(self, old: str, new: str) -> blosc2.NDArray:
        old_key = self._table_key(self._col_logical_key(old))
        if old_key in self._store.map_tree:
            moves = [(old_key, self._col_logical_key(new))]
            old_data_key = self._table_key(self._col_logical_key(old) + _UTF8_DATA_SUFFIX)
            if old_data_key in self._store.map_tree:
                moves.append((old_data_key, self._col_logical_key(new) + _UTF8_DATA_SUFFIX))
            new_dest = None
            for src_key, dst_logical in moves:
                dst_dest = self._dest_path(dst_logical, ".b2nd")
                src_dest = os.path.join(self._working_dir(), self._store.map_tree[src_key])
                os.makedirs(os.path.dirname(dst_dest), exist_ok=True)
                os.replace(src_dest, dst_dest)
                del self._store.map_tree[src_key]
                self._store.map_tree[self._table_key(dst_logical)] = os.path.relpath(
                    dst_dest, self._working_dir()
                ).replace(os.sep, "/")
                if new_dest is None:
                    new_dest = dst_dest
            self._store._modified = True
            return blosc2.open(new_dest, mode=self._mode)
        old_path = self._list_col_path(old)
        new_path = self._list_col_path(new)
        if os.path.exists(old_path):
            os.makedirs(os.path.dirname(new_path), exist_ok=True)
            os.replace(old_path, new_path)
            return blosc2.open(new_path, mode=self._mode)
        raise KeyError(old)

    # ------------------------------------------------------------------
    # TableStorage interface — index catalog and epoch counters
    # ------------------------------------------------------------------

    def load_index_catalog(self) -> dict:
        meta = self._open_meta()
        raw = meta.vlmeta.get("index_catalog")
        if not isinstance(raw, dict):
            return {}
        catalog = copy.deepcopy(raw)
        working_dir = self._working_dir()
        store = self._store
        for col_name, descriptor in catalog.items():
            catalog[col_name] = FileTableStorage._absolutize_descriptor(descriptor, working_dir)
            if store.is_zip_store and self._mode == "r":
                self._register_index_zip_paths(catalog[col_name])
        return catalog

    def _register_index_zip_paths(self, descriptor: dict) -> None:
        """Register sidecar paths from *descriptor* in the zip-offset registry."""
        from blosc2.indexing import _SIDECAR_ZIP_REGISTRY

        store = self._store
        working_dir = self._working_dir()
        for obj, key in FileTableStorage._walk_descriptor_paths(descriptor):
            abs_path = obj[key]
            if abs_path in _SIDECAR_ZIP_REGISTRY:
                continue
            rel = os.path.relpath(abs_path, working_dir).replace(os.sep, "/")
            info = store.offsets.get(rel)
            if info is None:
                continue
            _SIDECAR_ZIP_REGISTRY[abs_path] = (store.b2z_path, info["offset"])
            self._registered_sidecar_paths.append(abs_path)

    def save_index_catalog(self, catalog: dict) -> None:
        meta = self._open_meta()
        working_dir = self._working_dir()
        relativized = {
            col: FileTableStorage._relativize_descriptor(desc, working_dir) for col, desc in catalog.items()
        }
        meta.vlmeta["index_catalog"] = relativized
        self._bump_index_catalog_revision()

    def get_epoch_counters(self) -> tuple[int, int]:
        meta = self._open_meta()
        ve = int(meta.vlmeta.get("value_epoch", 0) or 0)
        vis_e = int(meta.vlmeta.get("visibility_epoch", 0) or 0)
        return ve, vis_e

    def bump_value_epoch(self) -> int:
        meta = self._open_meta()
        ve = int(meta.vlmeta.get("value_epoch", 0) or 0) + 1
        meta.vlmeta["value_epoch"] = ve
        return ve

    def bump_visibility_epoch(self) -> int:
        meta = self._open_meta()
        vis_e = int(meta.vlmeta.get("visibility_epoch", 0) or 0) + 1
        meta.vlmeta["visibility_epoch"] = vis_e
        return vis_e

    def index_anchor_path(self, col_name: str) -> str | None:
        table_rel = self._root_key.lstrip("/")
        return os.path.join(self._working_dir(), table_rel, _INDEXES_DIR, col_name, "_anchor")
