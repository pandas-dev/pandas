#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Internal variable-length UTF-8 string column backed by offsets + bytes.

This module is *not* part of the public API.  It provides row-wise string
semantics (one ``str`` value per row) stored Arrow-style as two companion
:class:`blosc2.NDArray` objects:

* **offsets** — ``int64``, length ``n + 1`` where ``n`` is the number of
  persisted rows.  ``offsets[0]`` is always ``0`` and ``offsets[i+1]`` is the
  end byte position of row ``i``.
* **data** — ``uint8``, the concatenated UTF-8 encoding of all row values.
  Its length is at least ``offsets[n]`` (one slack byte is kept when empty
  because zero-length NDArrays cannot be created).

Reading rows ``[a, b)`` needs ``offsets[a : b + 1]`` plus
``bytes[offsets[a] : offsets[b]]`` — both plain NDArray slice reads.

Nulls are represented with a per-column sentinel string (like every other
scalar CTable column); ``None`` written to a nullable column is converted to
the sentinel, and reads return the sentinel verbatim.

Bulk reads return :class:`numpy.dtypes.StringDType` arrays (NumPy >= 2.0),
which support vectorized comparison, ordering, and ``np.strings`` functions.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

# Pending rows are flushed to the backing arrays once either bound is hit.
# _FLUSH_ROWS is set well above what per-row overhead would suggest: each
# flush pays a fixed NDArray resize + slice-write cost, so fewer, larger
# flushes amortize that cost over more rows -- a taxi-like ASCII workload
# (1e7 rows, ~26-char average) measured a monotonic ingest speedup from
# 4096 (~3.6s) up to 65536 (~0.9s) with essentially no peak-memory increase,
# then diminishing returns beyond that (largely because _FLUSH_CHARS starts
# binding before _FLUSH_ROWS does).
_FLUSH_ROWS = 65536
_FLUSH_CHARS = 1 << 22  # ~4 Mi characters (>= 4 MiB encoded)

# Storage grids for freshly created backing arrays.  Both arrays are created
# tiny (shape (1,)) and grown by resize; the chunk shape is fixed at creation
# time, so it must be sized for the eventual data, not the initial shape.
_OFFSETS_CHUNKS = (2**17,)  # 1 MiB chunks of int64 row offsets
_DATA_CHUNKS = (2**21,)  # 2 MiB chunks of UTF-8 bytes

# Sparse gathers read the persisted region in clusters; a new cluster starts
# when the gap between consecutive row indices exceeds this many rows.
_GATHER_GAP = 1024

# Multiplier for the byte-hash in factorize_span (same mixer as
# groupby._factorize_fixed_width_str).
_HASH_MIX = np.uint64(0x9E3779B97F4A7C15)


def _factorize_byte_rows(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact factorization of the rows of a ``(k, L)`` uint8 matrix.

    Returns ``(rep_rows, inverse)`` where ``rep_rows`` holds one representative
    row index per distinct byte string (in unspecified group order — the
    caller re-sorts groups) and ``inverse`` maps each row to its group.
    Hashes each row into one uint64 and factorizes the integers, with a
    vectorized verify pass; on a hash collision, falls back to an exact sort
    of the raw rows.
    """
    h = mat[:, 0].astype(np.uint64)
    for i in range(1, mat.shape[1]):
        h = (h * _HASH_MIX) ^ mat[:, i]
    hash_uniques, inverse = np.unique(h, return_inverse=True)
    rep_rows = np.empty(len(hash_uniques), dtype=np.int64)
    rep_rows[inverse] = np.arange(len(mat))
    if not (mat == mat[rep_rows][inverse]).all():  # collision: exact fallback
        as_void = np.ascontiguousarray(mat).view([("", np.uint8, mat.shape[1])]).ravel()
        void_uniques, inverse = np.unique(as_void, return_inverse=True)
        rep_rows = np.empty(len(void_uniques), dtype=np.int64)
        rep_rows[inverse] = np.arange(len(mat))
    return rep_rows, inverse


def have_string_dtype() -> bool:
    """True when the installed NumPy provides ``StringDType`` (NumPy >= 2.0)."""
    return hasattr(np.dtypes, "StringDType")


def string_dtype():
    """Return a ``numpy.dtypes.StringDType`` instance, or raise if unavailable."""
    try:
        return np.dtypes.StringDType()
    except AttributeError:  # pragma: no cover - only on numpy < 2.0
        raise TypeError(
            "utf8 columns require NumPy >= 2.0 (numpy.dtypes.StringDType); "
            f"installed version is {np.__version__}. Use blosc2.vlstring() or "
            "blosc2.string(max_length=...) instead."
        ) from None


def _pack_utf8_kernel():
    """Return the compiled bulk StringDType packer, or ``None``.

    ``None`` means the ``utf8_ext`` extension is unavailable (a source
    install without a C toolchain, a WASM build, or a NumPy without the
    StringDType C API); callers fall back to a pure-Python per-row loop.
    """
    try:
        from blosc2 import utf8_ext
    except ImportError:
        return None
    return getattr(utf8_ext, "pack_utf8_span", None)


def _encode_utf8_kernel():
    """Return the compiled bulk UTF-8 encoder, or ``None``.

    ``None`` means the ``utf8_ext`` extension is unavailable; callers fall
    back to the pure-Python join+encode paths in ``_rewrite_from``.
    """
    try:
        from blosc2 import utf8_ext
    except ImportError:
        return None
    return getattr(utf8_ext, "encode_utf8_span", None)


def _new_backend_arrays(cparams=None, dparams=None, *, offsets_urlpath=None, data_urlpath=None):
    """Create fresh (offsets, data) NDArrays for an empty utf8 column."""
    import blosc2

    kwargs: dict[str, Any] = {}
    if cparams is not None:
        kwargs["cparams"] = cparams
    if dparams is not None:
        kwargs["dparams"] = dparams
    off_kwargs = dict(kwargs)
    data_kwargs = dict(kwargs)
    if offsets_urlpath is not None:
        off_kwargs["urlpath"] = offsets_urlpath
        off_kwargs["mode"] = "w"
    if data_urlpath is not None:
        data_kwargs["urlpath"] = data_urlpath
        data_kwargs["mode"] = "w"
    offsets = blosc2.zeros((1,), dtype=np.int64, chunks=_OFFSETS_CHUNKS, **off_kwargs)
    data = blosc2.zeros((1,), dtype=np.uint8, chunks=_DATA_CHUNKS, **data_kwargs)
    return offsets, data


class Utf8Array:
    """Row-wise variable-length UTF-8 string array over offsets + bytes NDArrays.

    Provides the row-oriented interface expected by CTable columns:
    ``append``, ``extend``, ``flush``, ``__len__``, ``__getitem__``, and
    ``__setitem__``.  Bulk reads return ``StringDType`` NumPy arrays; single
    reads return ``str``.

    In-place assignment to row ``i`` rewrites the byte blob and offsets of
    every row after ``i`` (a new value usually has a different byte length,
    which shifts all subsequent offsets), so ``__setitem__`` costs O(n - i).

    This class is internal; obtain instances via
    ``storage.create_varlen_scalar_column()`` or
    ``storage.open_varlen_scalar_column()`` with a ``Utf8Spec``.

    Parameters
    ----------
    spec:
        The :class:`~blosc2.schema.Utf8Spec` describing this column.
    offsets:
        ``int64`` NDArray of row offsets (length ``n + 1``).  Created fresh
        (in memory) when ``None``.
    data:
        ``uint8`` NDArray with the concatenated UTF-8 bytes.  Created fresh
        (in memory) when ``None``.
    """

    def __init__(self, spec, offsets=None, data=None) -> None:
        from blosc2.schema import Utf8Spec

        if not isinstance(spec, Utf8Spec):
            raise TypeError(f"Utf8Array requires a Utf8Spec, got {type(spec)!r}")
        self._dtype = string_dtype()
        self._spec = spec
        if (offsets is None) != (data is None):
            raise ValueError("offsets and data must be provided together")
        if offsets is None:
            offsets, data = _new_backend_arrays()
        self._offsets = offsets
        self._data = data
        self._persisted_rows: int = int(offsets.shape[0]) - 1
        # End byte position of the persisted region; resolved lazily because it
        # needs a chunk read from the offsets array.
        self._bytes_used_cache: int | None = 0 if self._persisted_rows == 0 else None
        # Rows not yet flushed to the backing arrays (list of str).
        self._pending: list[str] = []
        self._pending_chars: int = 0

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @property
    def _bytes_used(self) -> int:
        if self._bytes_used_cache is None:
            self._bytes_used_cache = int(self._offsets[self._persisted_rows])
        return self._bytes_used_cache

    def _coerce(self, value: Any) -> str:
        """Coerce *value* to ``str``, mapping ``None`` to the null sentinel."""
        if value is None:
            null_value = getattr(self._spec, "null_value", None)
            if null_value is None:
                raise TypeError("Column of utf8 strings is not nullable; received None.")
            return null_value
        if isinstance(value, str):
            return str(value)  # np.str_ instances are normalized to plain str
        raise TypeError(f"Expected str for utf8 column, got {type(value).__name__!r}.")

    def _flush_if_needed(self) -> None:
        if len(self._pending) >= _FLUSH_ROWS or self._pending_chars >= _FLUSH_CHARS:
            self.flush()

    def _read_persisted_span(self, a: int, b: int) -> np.ndarray:
        """Return persisted rows ``[a, b)`` as a StringDType array."""
        n = b - a
        if n <= 0:
            return np.empty(0, dtype=self._dtype)
        offs = np.asarray(self._offsets[a : b + 1], dtype=np.int64)
        start, end = int(offs[0]), int(offs[-1])
        rel = offs - start
        out = np.empty(n, dtype=self._dtype)
        kernel = _pack_utf8_kernel()
        if kernel is not None:
            data = (
                np.ascontiguousarray(self._data[start:end]) if end > start else np.zeros(1, dtype=np.uint8)
            )
            kernel(rel, data, out)
            return out
        blob = np.asarray(self._data[start:end]).tobytes() if end > start else b""
        for i in range(n):
            out[i] = blob[rel[i] : rel[i + 1]].decode("utf-8")
        return out

    def _gather_persisted(self, indices: np.ndarray) -> np.ndarray:
        """Gather persisted rows at *indices* (any order) as a StringDType array.

        Indices are read in sorted clusters of nearby rows so that a sparse
        gather (e.g. a few head and tail rows for display) does not read the
        whole column.
        """
        out = np.empty(len(indices), dtype=self._dtype)
        if len(indices) == 0:
            return out
        order = np.argsort(indices, kind="stable")
        sorted_idx = indices[order]
        cluster_starts = np.flatnonzero(np.diff(sorted_idx) > _GATHER_GAP) + 1
        pos = 0
        for cluster in np.split(sorted_idx, cluster_starts):
            lo, hi = int(cluster[0]), int(cluster[-1])
            span = self._read_persisted_span(lo, hi + 1)
            out[order[pos : pos + len(cluster)]] = span[cluster - lo]
            pos += len(cluster)
        return out

    def _read_span(self, a: int, b: int) -> np.ndarray:
        """Return rows ``[a, b)`` (persisted + pending) as a StringDType array."""
        np_rows = self._persisted_rows
        if b <= np_rows:
            return self._read_persisted_span(a, b)
        persisted = self._read_persisted_span(a, min(b, np_rows)) if a < np_rows else None
        pending = np.array(self._pending[max(0, a - np_rows) : b - np_rows], dtype=self._dtype)
        if persisted is None:
            return pending
        return np.concatenate([persisted, pending])

    def _get_many(self, indices: np.ndarray) -> np.ndarray:
        indices = np.asarray(indices)
        if indices.ndim != 1:
            indices = indices.ravel()
        n = len(self)
        indices = np.where(indices < 0, indices + n, indices).astype(np.int64, copy=False)
        m = len(indices)
        if m and (indices.min() < 0 or indices.max() >= n):
            raise IndexError("Utf8Array index out of range")
        if m and indices[-1] - indices[0] == m - 1 and bool((np.diff(indices) == 1).all()):
            # A contiguous ascending run (e.g. a full-column read routed here
            # via an index array rather than a step-1 slice) is just a span
            # read -- skip the sort/cluster/scatter-copy gather below, which
            # pays a full StringDType element-wise copy for no reason here.
            return self._read_span(int(indices[0]), int(indices[-1]) + 1)
        out = np.empty(m, dtype=self._dtype)
        np_rows = self._persisted_rows
        pending_mask = indices >= np_rows
        if pending_mask.any():
            out[pending_mask] = [self._pending[i - np_rows] for i in indices[pending_mask]]
        if not pending_mask.all():
            persisted_mask = ~pending_mask
            out[persisted_mask] = self._gather_persisted(indices[persisted_mask])
        return out

    def _rewrite_from(self, pos: int, values: list[str]) -> None:
        """Replace persisted rows ``pos ..`` with *values*, shifting offsets.

        ``pos + len(values)`` becomes the new persisted row count; the byte
        blob and offsets after ``pos`` are rewritten.
        """
        kernel = _encode_utf8_kernel()
        if kernel is not None and values:
            data_arr, lengths = kernel(values)
        elif values and all(v.isascii() for v in values):
            data_arr = np.frombuffer("".join(values).encode("ascii"), dtype=np.uint8)
            lengths = np.fromiter(map(len, values), dtype=np.int64, count=len(values))
        else:
            encoded = [v.encode("utf-8") for v in values]
            data_arr = np.frombuffer(b"".join(encoded), dtype=np.uint8)
            lengths = np.fromiter((len(e) for e in encoded), dtype=np.int64, count=len(encoded))
        if pos == 0:
            start = 0
        elif pos == self._persisted_rows:
            start = self._bytes_used
        else:
            start = int(self._offsets[pos])
        new_used = start + len(data_arr)
        new_rows = pos + len(values)
        if int(self._data.shape[0]) != max(new_used, 1):
            self._data.resize((max(new_used, 1),))
        if len(data_arr):
            self._data[start:new_used] = data_arr
        if int(self._offsets.shape[0]) != new_rows + 1:
            self._offsets.resize((new_rows + 1,))
        if values:
            self._offsets[pos + 1 : new_rows + 1] = start + np.cumsum(lengths)
        self._persisted_rows = new_rows
        self._bytes_used_cache = new_used

    # ------------------------------------------------------------------
    # Public write interface
    # ------------------------------------------------------------------

    def append(self, value: Any) -> None:
        """Append one string row (``None`` maps to the null sentinel)."""
        value = self._coerce(value)
        self._pending.append(value)
        self._pending_chars += len(value)
        self._flush_if_needed()

    def extend(self, values: Iterable[Any]) -> None:
        """Append many string rows.

        Processed in ``_FLUSH_ROWS``-sized chunks so the char-count flush
        bound is only checked once per chunk rather than once per row; an
        unusual batch of many multi-MB strings can therefore overshoot
        ``_FLUSH_CHARS`` by up to one chunk before a flush is triggered.
        """
        it = iter(values)
        while True:
            chunk = list(itertools.islice(it, _FLUSH_ROWS))
            if not chunk:
                break
            if all(type(v) is str for v in chunk):
                self._pending.extend(chunk)
                self._pending_chars += sum(map(len, chunk))
            else:
                for v in chunk:
                    v = self._coerce(v)
                    self._pending.append(v)
                    self._pending_chars += len(v)
            self._flush_if_needed()

    def flush(self) -> None:
        """Write pending rows to the backing offsets/data NDArrays."""
        if not self._pending:
            return
        values, self._pending = self._pending, []
        self._pending_chars = 0
        self._rewrite_from(self._persisted_rows, values)

    def set_all(self, values: Iterable[Any]) -> None:
        """Replace the whole column content in one bulk write.

        Writes through the existing backing offsets/data NDArrays, so a
        store-backed column stays persistent (unlike building a fresh
        in-memory ``Utf8Array``).  Used by ``sort_by(inplace=True)`` and
        ``compact()`` to rewrite a column in a new row order.
        """
        coerced = [self._coerce(v) for v in values]
        self._pending = []
        self._pending_chars = 0
        self._rewrite_from(0, coerced)

    # ------------------------------------------------------------------
    # Public read interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self._persisted_rows + len(self._pending)

    def __iter__(self) -> Iterator[str]:
        yield from self[:]

    def __getitem__(self, index: int | slice | list | tuple | np.ndarray):
        if isinstance(index, (int, np.integer)):
            n = len(self)
            index = int(index)
            if index < 0:
                index += n
            if not (0 <= index < n):
                raise IndexError("Utf8Array index out of range")
            if index >= self._persisted_rows:
                return self._pending[index - self._persisted_rows]
            return str(self._read_persisted_span(index, index + 1)[0])

        if isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            if step == 1:
                return self._read_span(start, stop)
            return self._get_many(np.arange(start, stop, step, dtype=np.int64))

        if isinstance(index, np.ndarray) and index.dtype == np.bool_:
            if len(index) != len(self):
                raise IndexError(f"Boolean mask length {len(index)} does not match array length {len(self)}")
            return self._get_many(np.flatnonzero(index))

        if isinstance(index, (list, tuple, np.ndarray)):
            return self._get_many(np.asarray(index, dtype=np.int64))

        raise TypeError(f"Utf8Array indices must be int, slice, or array; got {type(index)!r}")

    def __setitem__(self, index: int, value: Any) -> None:
        """Overwrite the value at *index*.

        Because row values have variable byte lengths, overwriting a persisted
        row rewrites the byte blob and offsets of all subsequent rows —
        an O(n - index) operation.
        """
        if not isinstance(index, (int, np.integer)):
            raise TypeError(f"Utf8Array assignment index must be int, got {type(index)!r}")
        value = self._coerce(value)
        n = len(self)
        index = int(index)
        if index < 0:
            index += n
        if not (0 <= index < n):
            raise IndexError("Utf8Array index out of range")
        if index >= self._persisted_rows:
            self._pending[index - self._persisted_rows] = value
            return
        # The tail bytes themselves are unchanged; only their position moves
        # by the length delta of the replaced value, so shift raw bytes and
        # add the delta to the tail offsets instead of decoding/re-encoding
        # every following row.
        encoded = value.encode("utf-8")
        bounds = np.asarray(self._offsets[index : index + 2], dtype=np.int64)
        old_start, old_end = int(bounds[0]), int(bounds[1])
        delta = len(encoded) - (old_end - old_start)
        old_used = self._bytes_used
        new_used = old_used + delta
        tail = np.asarray(self._data[old_end:old_used]).copy() if delta != 0 and old_used > old_end else None
        if delta > 0 and int(self._data.shape[0]) < new_used:
            self._data.resize((new_used,))
        if encoded:
            self._data[old_start : old_start + len(encoded)] = np.frombuffer(encoded, dtype=np.uint8)
        if tail is not None:
            self._data[old_end + delta : new_used] = tail
        if delta < 0 and int(self._data.shape[0]) != max(new_used, 1):
            self._data.resize((max(new_used, 1),))
        if delta != 0:
            n = self._persisted_rows
            self._offsets[index + 1 : n + 1] = (
                np.asarray(self._offsets[index + 1 : n + 1], dtype=np.int64) + delta
            )
            self._bytes_used_cache = new_used

    # ------------------------------------------------------------------
    # Properties mirroring the interface expected by CTable
    # ------------------------------------------------------------------

    @property
    def spec(self):
        return self._spec

    @property
    def dtype(self):
        """The ``StringDType`` used for materialized reads."""
        return self._dtype

    @property
    def offsets(self):
        """The underlying ``int64`` NDArray of row offsets (length ``n + 1``)."""
        return self._offsets

    @property
    def data(self):
        """The underlying ``uint8`` NDArray with the concatenated UTF-8 bytes."""
        return self._data

    @property
    def schunk(self):
        return self._offsets.schunk

    @property
    def urlpath(self) -> str | None:
        return getattr(self._offsets, "urlpath", None)

    @property
    def nbytes(self) -> int:
        return self._offsets.schunk.nbytes + self._data.schunk.nbytes

    @property
    def cbytes(self) -> int:
        return self._offsets.schunk.cbytes + self._data.schunk.cbytes

    @property
    def cratio(self) -> float:
        cb = self.cbytes
        if cb == 0:
            return float("inf")
        return self.nbytes / cb

    def factorizer(self) -> Utf8Factorizer:
        """Return a fresh incremental factorizer over this column's rows."""
        return Utf8Factorizer(self)

    def factorize_span(self, a: int, b: int) -> tuple[np.ndarray, np.ndarray]:
        """Factorize rows ``[a, b)`` without decoding them.

        Returns ``(codes, uniques)`` where ``uniques`` is a ``StringDType``
        array of the distinct values sorted ascending and ``codes`` (int64,
        length ``b - a``) maps each row to its value — the same contract as
        ``np.unique(values, return_inverse=True)``, but computed from the raw
        offsets/bytes buffers via :class:`Utf8Factorizer`: only the distinct
        values are ever decoded to ``str``.  Pending rows are flushed first.
        """
        fact = self.factorizer()
        codes = fact.codes_for_span(a, b)
        uniques = fact.uniques()
        order = np.argsort(uniques, kind="stable")
        rank = np.empty(len(order), dtype=np.int64)
        rank[order] = np.arange(len(order))
        return rank[codes], uniques[order]

    def equal_mask_span(self, value: str, a: int, b: int) -> np.ndarray:
        """Boolean mask for persisted rows ``[a, b)``: row bytes == *value*'s UTF-8 bytes.

        Compares raw bytes with no decode: rows are first filtered by byte
        length, then the surviving candidates are compared byte-by-byte with
        whole-array gathers (one comparison per byte position of *value*).
        """
        self.flush()
        enc = value.encode("utf-8")
        length = len(enc)
        target = np.frombuffer(enc, dtype=np.uint8)
        offs = np.asarray(self._offsets[a : b + 1], dtype=np.int64)
        rel = offs - int(offs[0])
        data = np.asarray(self._data[int(offs[0]) : int(offs[-1])])
        lengths = np.diff(rel)
        mask = lengths == length
        if length and mask.any():
            cand = np.flatnonzero(mask)
            idx = rel[cand]
            hit = np.ones(len(idx), dtype=np.bool_)
            for i in range(length):
                hit &= data[idx] == target[i]
                idx = idx + 1
            mask[cand[~hit]] = False
        return mask

    def order_masks_span(self, value: str, a: int, b: int) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(lt, gt)`` boolean masks comparing rows ``[a, b)`` to *value*.

        Byte-lexicographic comparison, which equals Unicode code-point order
        for valid UTF-8.  A row where neither mask is True is equal to
        *value*.  Rows are grouped by byte length (each row touched once);
        within a group, bytes are compared position by position against
        *value*'s bytes with whole-array gathers.
        """
        self.flush()
        target = np.frombuffer(value.encode("utf-8"), dtype=np.uint8)
        probe_len = len(target)
        n = b - a
        lt = np.zeros(n, dtype=np.bool_)
        gt = np.zeros(n, dtype=np.bool_)
        if n == 0:
            return lt, gt
        offs = np.asarray(self._offsets[a : b + 1], dtype=np.int64)
        rel = offs - int(offs[0])
        data = np.asarray(self._data[int(offs[0]) : int(offs[-1])])
        lengths = np.diff(rel)
        starts = rel[:-1]
        if int(lengths.max()) <= 65536:
            distinct_lengths = np.flatnonzero(np.bincount(lengths))
        else:
            distinct_lengths = np.unique(lengths)
        for length in distinct_lengths:
            length = int(length)
            rows = np.flatnonzero(lengths == length)
            m = min(length, probe_len)
            undecided = np.ones(len(rows), dtype=np.bool_)
            row_lt = np.zeros(len(rows), dtype=np.bool_)
            row_gt = np.zeros(len(rows), dtype=np.bool_)
            idx = starts[rows]
            for i in range(m):
                byte = data[idx]
                row_lt |= undecided & (byte < target[i])
                row_gt |= undecided & (byte > target[i])
                undecided &= byte == target[i]
                idx = idx + 1
            if length < probe_len:
                row_lt |= undecided  # row is a strict byte-prefix of value
            elif length > probe_len:
                row_gt |= undecided  # value is a strict byte-prefix of row
            lt[rows] = row_lt
            gt[rows] = row_gt
        return lt, gt

    def arrow_slice(self, pa, a: int, b: int, null_value: str | None = None):
        """Persisted rows ``[a, b)`` as a ``pyarrow.LargeStringArray``.

        The storage layout (int64 offsets + UTF-8 byte blob) is exactly
        Arrow's ``large_string`` layout, so the array is built directly from
        the raw buffers with no per-row decode or Python string objects.
        When *null_value* is given, rows equal to the sentinel become Arrow
        nulls (matched on raw bytes, still without decoding).
        """
        n = b - a
        offs = np.ascontiguousarray(self._offsets[a : b + 1], dtype=np.int64)
        start, end = int(offs[0]), int(offs[-1])
        rel = offs - start
        data = np.ascontiguousarray(self._data[start:end]) if end > start else np.empty(0, dtype=np.uint8)
        validity = None
        null_count = 0
        if null_value is not None and n > 0:
            mask = self.equal_mask_span(null_value, a, b)
            null_count = int(mask.sum())
            if null_count:
                validity = pa.array(~mask).buffers()[1]
        return pa.LargeStringArray.from_buffers(
            n, pa.py_buffer(rel), pa.py_buffer(data), validity, null_count
        )

    def copy(self, spec=None, **kwargs: Any) -> Utf8Array:
        """Return an in-memory copy."""
        if spec is None:
            spec = self._spec
        out = Utf8Array(spec)
        out.extend(self)
        out.flush()
        return out


class Utf8Factorizer:
    """Incremental factorizer over a :class:`Utf8Array`'s rows.

    Successive :meth:`codes_for_span` calls share a vocabulary: each row is
    hashed from its raw bytes and matched against the known values of its
    byte length (``np.searchsorted`` on sorted hashes plus a vectorized
    byte-equality verify), so only rows carrying *new* values pay a sort
    (via :func:`_factorize_byte_rows`).  Codes are global across spans, in
    first-appearance order; :meth:`uniques` returns the decoded vocabulary in
    that same order.  No row is ever decoded — only the distinct values.

    A hash collision between a new value and an already-known value of the
    same byte length is not resolved: the new value always gets a fresh
    code, even in the (unhandled) case where its decoded string equals an
    existing vocabulary entry's. No downstream consumer merges groups by
    decoded value, so such a collision would silently split what should be
    one group into two rather than being caught -- an accepted, unmitigated
    risk given how vanishingly rare it is with 64-bit hashes, not a
    guarantee that it is otherwise resolved.
    """

    def __init__(self, arr: Utf8Array) -> None:
        arr.flush()
        self._arr = arr
        self._values: list[str] = []
        self._empty_code: int | None = None
        # byte length -> (sorted hashes, codes aligned to them, (D, L) rep bytes)
        self._by_len: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def uniques(self) -> np.ndarray:
        """The vocabulary so far, decoded, in global code order."""
        return np.array(self._values, dtype=self._arr._dtype)

    def codes_for_span(self, a: int, b: int) -> np.ndarray:
        """Global codes (int64) for persisted rows ``[a, b)``."""
        n = b - a
        codes = np.empty(max(n, 0), dtype=np.int64)
        if n <= 0:
            return codes
        offs = np.asarray(self._arr._offsets[a : b + 1], dtype=np.int64)
        start, end = int(offs[0]), int(offs[-1])
        rel = offs - start
        data = np.asarray(self._arr._data[start:end]) if end > start else np.empty(0, dtype=np.uint8)
        lengths = np.diff(rel)
        if int(lengths.max()) <= 65536:
            distinct_lengths = np.flatnonzero(np.bincount(lengths))
        else:
            distinct_lengths = np.unique(lengths)
        for length in distinct_lengths:
            length = int(length)
            rows = np.flatnonzero(lengths == length)
            if length == 0:
                if self._empty_code is None:
                    self._empty_code = len(self._values)
                    self._values.append("")
                codes[rows] = self._empty_code
                continue
            # Column-wise gather: reusing one index vector is ~2x faster than
            # materializing the (k, L) int64 index matrix of a 2-D gather.
            mat = np.empty((len(rows), length), dtype=np.uint8)
            idx = rel[rows]
            for i in range(length):
                mat[:, i] = data[idx]
                idx += 1
            h = mat[:, 0].astype(np.uint64)
            for i in range(1, length):
                h = (h * _HASH_MIX) ^ mat[:, i]
            miss = np.ones(len(rows), dtype=np.bool_)
            entry = self._by_len.get(length)
            if entry is not None:
                sorted_h, code_at, rep_mat = entry
                pos = np.searchsorted(sorted_h, h)
                pos[pos == len(sorted_h)] = 0  # equality check below rejects
                hit = sorted_h[pos] == h
                if hit.any():
                    # One full-row compare beats gathering the hit rows out of
                    # mat first (hits are usually nearly all rows).
                    good = np.flatnonzero(hit & (mat == rep_mat[pos]).all(axis=1))
                    codes[rows[good]] = code_at[pos[good]]
                    miss[good] = False
            if miss.any():
                miss_mat = mat[miss]
                rep_rows, inverse = _factorize_byte_rows(miss_mat)
                base = len(self._values)
                self._values.extend(bytes(miss_mat[j]).decode("utf-8") for j in rep_rows)
                codes[rows[miss]] = base + inverse
                new_h = np.empty(len(rep_rows), dtype=np.uint64)
                new_rep = miss_mat[rep_rows]
                new_h[:] = new_rep[:, 0]
                for i in range(1, length):
                    new_h = (new_h * _HASH_MIX) ^ new_rep[:, i]
                new_codes = np.arange(base, base + len(rep_rows), dtype=np.int64)
                if entry is not None:
                    new_h = np.concatenate([entry[0], new_h])
                    new_codes = np.concatenate([entry[1], new_codes])
                    new_rep = np.concatenate([entry[2], new_rep])
                order = np.argsort(new_h, kind="stable")
                self._by_len[length] = (new_h[order], new_codes[order], new_rep[order])
        return codes
