#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Internal variable-length scalar column adapter over BatchArray.

This module is *not* part of the public API.  It provides row-wise scalar
semantics (one str, bytes, struct dict, or schema-less object value per row) backed by batched
msgpack storage via :class:`blosc2.BatchArray`.

Physical layout: each chunk in the backing BatchArray stores a list of
scalar values, e.g. ``["foo", None, "bar", "baz"]``.  Nulls are represented
as native Python ``None`` and never converted to a sentinel value.
"""

from __future__ import annotations

import os
from bisect import bisect_right
from collections import defaultdict
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

from blosc2.batch_array import BatchArray

# Meta tag written on the backing BatchArray's SChunk so the storage layer
# can identify the role of this container on reopen.
_CTABLE_VARLEN_SCALAR_META_KEY = "ctable_varlen_scalar"


def _role_metadata_for_spec(spec) -> dict[str, Any]:
    """Return the fixed metadata role tag for a CTable varlen scalar backend."""
    if spec.python_type is str:
        py_type = "str"
    elif spec.python_type is bytes:
        py_type = "bytes"
    elif spec.python_type is dict:
        py_type = "struct"
    else:
        py_type = "object"
    return {
        "version": 1,
        "py_type": py_type,
        "nullable": bool(getattr(spec, "nullable", False)),
        "batch_rows": getattr(spec, "batch_rows", 2048),
    }


def _storage_with_role_meta(spec, **storage_kwargs: Any):
    import blosc2

    storage = blosc2.Storage(**storage_kwargs)
    fixed_meta = dict(storage.meta or {})
    fixed_meta[_CTABLE_VARLEN_SCALAR_META_KEY] = _role_metadata_for_spec(spec)
    storage.meta = fixed_meta
    return storage


def _validate_role_metadata(backend: BatchArray, spec) -> None:
    """Validate the optional CTable varlen scalar role tag on an opened backend."""
    meta = backend.schunk.meta
    if _CTABLE_VARLEN_SCALAR_META_KEY not in meta:
        # Older local artifacts may only have the BatchArray tag; the CTable schema
        # still identifies the logical role, so keep reopen tolerant.
        return
    role = meta[_CTABLE_VARLEN_SCALAR_META_KEY]
    if spec.python_type is str:
        expected_py_type = "str"
    elif spec.python_type is bytes:
        expected_py_type = "bytes"
    elif spec.python_type is dict:
        expected_py_type = "struct"
    else:
        expected_py_type = "object"
    if role.get("py_type") != expected_py_type:
        raise ValueError(
            f"Varlen scalar backend type mismatch: expected {expected_py_type!r}, "
            f"found {role.get('py_type')!r}."
        )


def _make_backend(spec) -> BatchArray:
    """Create a fresh in-memory BatchArray for a varlen scalar spec."""
    storage = _storage_with_role_meta(spec)
    return BatchArray(
        storage=storage,
        items_per_block=getattr(spec, "items_per_block", None),
        serializer=getattr(spec, "serializer", "msgpack"),
    )


def _make_persistent_backend(spec, urlpath: str, mode: str, *, cparams=None, dparams=None) -> BatchArray:
    """Create or open a persistent BatchArray for a varlen scalar spec."""
    os.makedirs(os.path.dirname(urlpath), exist_ok=True)
    kwargs: dict[str, Any] = {}
    if cparams is not None:
        kwargs["cparams"] = cparams
    if dparams is not None:
        kwargs["dparams"] = dparams
    storage = _storage_with_role_meta(spec, urlpath=urlpath, mode=mode, contiguous=True)
    return BatchArray(
        storage=storage,
        items_per_block=getattr(spec, "items_per_block", None),
        serializer=getattr(spec, "serializer", "msgpack"),
        **kwargs,
    )


def _open_persistent_backend(urlpath: str, mode: str, spec=None) -> BatchArray:
    """Reopen an existing persistent BatchArray (any mode)."""
    backend = BatchArray(urlpath=urlpath, mode=mode)
    if spec is not None:
        _validate_role_metadata(backend, spec)
    return backend


class _ScalarVarLenArray:
    """Row-wise variable-length scalar array backed by a :class:`~blosc2.BatchArray`.

    Provides the same row-oriented interface expected by CTable columns:
    ``append``, ``extend``, ``flush``, ``__len__``, ``__getitem__``, and
    ``__setitem__``.

    This class is internal; do not use it directly.

    Parameters
    ----------
    spec:
        A :class:`~blosc2.schema.VLStringSpec`,
        :class:`~blosc2.schema.VLBytesSpec`,
        :class:`~blosc2.schema.StructSpec`, or
        :class:`~blosc2.schema.ObjectSpec` describing this column.
    backend:
        Pre-constructed :class:`~blosc2.BatchArray`.  If ``None``, a fresh
        in-memory backend is created from *spec*.
    """

    def __init__(self, spec, backend: BatchArray | None = None) -> None:
        from blosc2.schema import ObjectSpec, StructSpec, VLBytesSpec, VLStringSpec

        if not isinstance(spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec)):
            raise TypeError(
                "_ScalarVarLenArray requires a VLStringSpec, VLBytesSpec, StructSpec, or "
                f"ObjectSpec, got {type(spec)!r}"
            )
        self._spec = spec
        self._py_type: type = spec.python_type  # str, bytes, dict, or object
        self._nullable: bool = getattr(spec, "nullable", False)
        self._batch_rows: int = int(getattr(spec, "batch_rows", 2048) or 2048)

        if backend is None:
            backend = _make_backend(spec)
        self._backend: BatchArray = backend

        # Pending rows not yet flushed to the backend.
        self._pending: list[Any] = []
        # Cumulative row count flushed into the backend (sum of all batch lengths).
        self._persisted_row_count: int = self._compute_persisted_rows()
        # Cache for prefix sums over batch lengths (invalidated on flush/setitem).
        self._prefix_cache: list[int] | None = None

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _compute_persisted_rows(self) -> int:
        """Compute the total number of rows persisted in the backend."""
        if len(self._backend) == 0:
            return 0
        return sum(self._backend._load_or_compute_batch_lengths())

    def _persisted_prefix_sums(self) -> list[int]:
        """Return a list of cumulative batch start positions (length = n_batches + 1)."""
        if self._prefix_cache is not None:
            return self._prefix_cache
        lengths = self._backend._load_or_compute_batch_lengths()
        prefix: list[int] = [0]
        total = 0
        for length in lengths:
            total += int(length)
            prefix.append(total)
        self._prefix_cache = prefix
        return prefix

    def _invalidate_prefix_cache(self) -> None:
        self._prefix_cache = None

    def _locate_persisted_row(self, row_index: int) -> tuple[int, int]:
        """Return (batch_index, inner_index) for *row_index* in the persisted region."""
        prefix = self._persisted_prefix_sums()
        batch_index = bisect_right(prefix, row_index) - 1
        inner_index = row_index - prefix[batch_index]
        return batch_index, inner_index

    def _get_batch_items(self, batch_index: int) -> list[Any]:
        """Return the items in batch *batch_index* as a plain Python list."""
        return self._backend[batch_index][:]

    def _coerce(self, value: Any) -> Any:
        """Coerce *value* to the column's Python type, respecting nullability."""
        if value is None:
            if not self._nullable:
                raise TypeError(f"Column {self._py_type.__name__!r} is not nullable; received None.")
            return None
        if self._py_type is str:
            if isinstance(value, str):
                return value
            raise TypeError(f"Expected str for vlstring column, got {type(value).__name__!r}.")
        if self._py_type is bytes:
            if isinstance(value, (bytes, bytearray, memoryview)):
                return bytes(value)
            raise TypeError(f"Expected bytes for vlbytes column, got {type(value).__name__!r}.")
        if self._py_type is dict:
            from blosc2.list_array import _coerce_struct_item

            return _coerce_struct_item(self._spec, value)
        return value

    def _flush_full_batches(self) -> None:
        """Flush as many full batches as possible from _pending."""
        while len(self._pending) >= self._batch_rows:
            batch = self._pending[: self._batch_rows]
            self._backend.append(batch)
            self._pending = self._pending[self._batch_rows :]
            self._persisted_row_count += len(batch)
            self._invalidate_prefix_cache()

    # ------------------------------------------------------------------
    # Public write interface
    # ------------------------------------------------------------------

    def append(self, value: Any) -> None:
        """Append one scalar row."""
        self._pending.append(self._coerce(value))
        self._flush_full_batches()

    def extend(self, values: Iterable[Any]) -> None:
        """Append many scalar rows."""
        for v in values:
            self._pending.append(self._coerce(v))
            if len(self._pending) >= self._batch_rows:
                self._flush_full_batches()

    def flush(self) -> None:
        """Flush any remaining pending rows to the backend as one batch."""
        if self._pending:
            batch = list(self._pending)
            self._backend.append(batch)
            self._persisted_row_count += len(batch)
            self._pending.clear()
            self._invalidate_prefix_cache()

    # ------------------------------------------------------------------
    # Public read interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self._persisted_row_count + len(self._pending)

    def __iter__(self) -> Iterator[Any]:
        yield from self[:]

    def __getitem__(self, index: int | slice | list | tuple) -> Any | list[Any]:
        if isinstance(index, int):
            n = len(self)
            if index < 0:
                index += n
            if not (0 <= index < n):
                raise IndexError("_ScalarVarLenArray index out of range")
            if index >= self._persisted_row_count:
                return self._pending[index - self._persisted_row_count]
            batch_index, inner_index = self._locate_persisted_row(index)
            return self._get_batch_items(batch_index)[inner_index]

        if isinstance(index, slice):
            indices = list(range(*index.indices(len(self))))
            return self._get_many(indices)

        # numpy array or list/tuple of indices
        try:
            import numpy as np

            if isinstance(index, np.ndarray):
                if index.dtype == bool:
                    if len(index) != len(self):
                        raise IndexError(
                            f"Boolean mask length {len(index)} does not match array length {len(self)}"
                        )
                    return self._get_many(np.flatnonzero(index).tolist())
                return self._get_many(index.tolist())
        except ImportError:
            pass
        if isinstance(index, (list, tuple)):
            return self._get_many(list(index))

        raise TypeError(f"_ScalarVarLenArray indices must be int, slice, or array; got {type(index)!r}")

    def __setitem__(self, index: int, value: Any) -> None:
        value = self._coerce(value)
        n = len(self)
        if index < 0:
            index += n
        if not (0 <= index < n):
            raise IndexError("_ScalarVarLenArray index out of range")
        if index >= self._persisted_row_count:
            self._pending[index - self._persisted_row_count] = value
            return
        # Rewrite the persisted batch.
        batch_index, inner_index = self._locate_persisted_row(index)
        items = self._get_batch_items(batch_index)
        items[inner_index] = value
        self._backend[batch_index] = items
        self._invalidate_prefix_cache()

    # ------------------------------------------------------------------
    # Bulk access helpers
    # ------------------------------------------------------------------

    def _get_many_grouped(self, indices: list[int]) -> list[Any]:
        out: list[Any] = [None] * len(indices)
        grouped: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for out_i, index in enumerate(indices):
            if index >= self._persisted_row_count:
                out[out_i] = self._pending[index - self._persisted_row_count]
            else:
                batch_index, inner_index = self._locate_persisted_row(index)
                grouped[batch_index].append((out_i, inner_index))
        for batch_index, refs in grouped.items():
            items = self._get_batch_items(batch_index)
            for out_i, inner_index in refs:
                out[out_i] = items[inner_index]
        return out

    def _get_many_monotonic(self, indices: list[int]) -> list[Any]:
        out: list[Any] = [None] * len(indices)
        prefix = self._persisted_prefix_sums()
        batch_index = 0
        batch_items: list[Any] | None = None

        i = 0
        while i < len(indices):
            index = indices[i]
            if index >= self._persisted_row_count:
                pending_start = index - self._persisted_row_count
                j = i + 1
                while (
                    j < len(indices)
                    and indices[j] >= self._persisted_row_count
                    and indices[j] == indices[j - 1] + 1
                ):
                    j += 1
                span = j - i
                out[i:j] = self._pending[pending_start : pending_start + span]
                i = j
                continue

            while batch_index + 1 < len(prefix) and index >= prefix[batch_index + 1]:
                batch_index += 1
                batch_items = None
            if batch_items is None:
                batch_items = self._get_batch_items(batch_index)

            batch_start = prefix[batch_index]
            batch_end = prefix[batch_index + 1]
            local_start = index - batch_start
            j = i + 1
            while j < len(indices) and indices[j] == indices[j - 1] + 1 and indices[j] < batch_end:
                j += 1
            span = j - i
            out[i:j] = batch_items[local_start : local_start + span]
            i = j

        return out

    def _get_many(self, indices: list[int]) -> list[Any]:
        if len(indices) <= 1:
            return self._get_many_grouped(indices)
        # Check if monotonic (allows faster sequential scan)
        monotonic = all(indices[k] < indices[k + 1] for k in range(len(indices) - 1))
        if monotonic:
            return self._get_many_monotonic(indices)
        return self._get_many_grouped(indices)

    # ------------------------------------------------------------------
    # Properties mirroring NDArray / ListArray interface expected by CTable
    # ------------------------------------------------------------------

    @property
    def dtype(self):
        """Always ``None`` for varlen scalar columns (no fixed NumPy dtype)."""
        return None

    @property
    def schunk(self):
        return self._backend.schunk

    @property
    def urlpath(self) -> str | None:
        return self._backend.urlpath

    @property
    def nbytes(self) -> int:
        return self._backend.nbytes

    @property
    def cbytes(self) -> int:
        return self._backend.cbytes

    @property
    def cratio(self) -> float:
        return self._backend.cratio

    def copy(self, spec=None, **kwargs: Any) -> _ScalarVarLenArray:
        """Return an in-memory copy."""
        if spec is None:
            spec = self._spec
        out = _ScalarVarLenArray(spec)
        out.extend(self)
        out.flush()
        return out
