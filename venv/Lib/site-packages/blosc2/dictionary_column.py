#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

"""Dictionary-encoded string column for CTable.

Physical layout
---------------
A dictionary column is stored as two components:

* **codes** — a fixed-width ``int32`` NDArray with one code per physical row
  slot.  The special code ``null_code`` (default ``-1``) marks null slots.
* **dict_store** — a variable-length string array (:class:`_ScalarVarLenArray`)
  holding unique category values in first-seen order.

An in-memory mapping ``_value_to_code: dict[str, int]`` is built lazily from
the persisted dict_store on open and kept in sync during writes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from blosc2.scalar_array import _ScalarVarLenArray
    from blosc2.schema import DictionarySpec


_NULL_INT32 = np.int32(-1)


class DictionaryColumn:
    """Row-wise dictionary-encoded string column wrapping codes + dict_store.

    This class is internal; obtain instances via
    ``storage.create_dictionary_column()`` or ``storage.open_dictionary_column()``.

    Parameters
    ----------
    spec:
        The :class:`~blosc2.schema.DictionarySpec` that describes this column.
    codes:
        A ``blosc2.NDArray`` of dtype ``int32`` with one slot per physical row.
    dict_store:
        A :class:`~blosc2.scalar_array._ScalarVarLenArray` holding unique
        category strings in insertion order (no nulls).
    """

    def __init__(self, spec: DictionarySpec, codes, dict_store: _ScalarVarLenArray) -> None:
        self._spec = spec
        self._codes = codes  # int32 NDArray (physical slot array)
        self._dict_store = dict_store  # _ScalarVarLenArray of vlstring (unique values)
        # Cache: str → int32 code.  Built lazily from dict_store on first access.
        self._value_to_code: dict[str, int] | None = None

    # ------------------------------------------------------------------
    # Cache management
    # ------------------------------------------------------------------

    def _ensure_cache(self) -> None:
        """Build the value→code mapping from the persisted dict_store."""
        if self._value_to_code is not None:
            return
        self._dict_store.flush()
        cache: dict[str, int] = {}
        for code, value in enumerate(self._dict_store):
            if value is not None:
                cache[value] = code
        self._value_to_code = cache

    def _invalidate_cache(self) -> None:
        self._value_to_code = None

    # ------------------------------------------------------------------
    # Encoding / decoding
    # ------------------------------------------------------------------

    def encode(self, value: str | None) -> int:
        """Encode *value* to an int32 code.  Appends new values to the dictionary."""
        if value is None:
            if not self._spec.nullable:
                raise ValueError(f"Dictionary column {self._spec!r} is not nullable; received None.")
            return self._spec.null_code
        if not isinstance(value, str):
            raise TypeError(f"Dictionary column expects str or None values, got {type(value).__name__!r}.")
        self._ensure_cache()
        assert self._value_to_code is not None
        code = self._value_to_code.get(value)
        if code is not None:
            return code
        # New category — append to dictionary.
        new_code = len(self._value_to_code)
        if new_code > np.iinfo(np.int32).max:
            raise OverflowError(
                "Dictionary column has exceeded the maximum number of unique values (2^31 - 1)."
            )
        self._dict_store.append(value)
        self._value_to_code[value] = new_code
        return new_code

    def decode(self, code: int) -> str | None:
        """Decode an int32 *code* to its string value, or ``None`` for null codes."""
        if code == self._spec.null_code:
            return None
        self._ensure_cache()
        return self._dict_store[int(code)]

    def decode_batch(self, codes) -> list[str | None]:
        """Decode an array of int32 *codes* to a list of strings (``None`` for null codes).

        Reads the whole dictionary store once (cardinality ``D``, not the ``N``
        rows) instead of indexing the backing var-length array per code, which
        decompresses a whole msgpack batch on every call.  For many codes this
        is dramatically cheaper than looping over :meth:`decode`.
        """
        codes = np.asarray(codes)
        self._dict_store.flush()
        all_strings = np.asarray(self._dict_store[:])  # D unique values, no nulls
        null_code = int(self._spec.null_code)
        result: list[str | None] = [None] * len(codes)
        non_null_idx = np.nonzero(codes != null_code)[0]
        if non_null_idx.size:
            picked = all_strings[codes[non_null_idx].astype(np.intp)].tolist()
            for pos, value in zip(non_null_idx.tolist(), picked, strict=True):
                result[pos] = value
        return result

    def encode_batch(self, values) -> np.ndarray:
        """Encode a sequence of str/None to a numpy ``int32`` array of codes."""
        result = np.empty(len(values), dtype=np.int32)
        for i, v in enumerate(values):
            result[i] = self.encode(v)
        return result

    def value_to_code(self, value: str) -> int:
        """Return the code for *value*.  Raises :exc:`KeyError` if absent."""
        self._ensure_cache()
        assert self._value_to_code is not None
        if value not in self._value_to_code:
            raise KeyError(value)
        return self._value_to_code[value]

    def code_to_value(self, code: int) -> str | None:
        """Return the category string for *code*."""
        return self.decode(code)

    # ------------------------------------------------------------------
    # Arrow-optimised batch import
    # ------------------------------------------------------------------

    def extend_from_arrow(self, pa, arrow_col, pos: int, m: int, *, ordered: bool = False) -> None:
        """Write *m* rows from an Arrow dictionary array into the codes NDArray at *pos*.

        Performs global dictionary unification: chunk-local codes are remapped
        to global codes.  ``ordered=True`` raises if chunk dictionary order
        differs from the established global order.
        """
        local_dict = arrow_col.dictionary.to_pylist()

        # Build a local-code → global-code lookup table.  The chunk-local
        # dictionary is small (one entry per distinct value in the batch), so
        # this Python loop is cheap; the per-row translation below is vectorised.
        lut = np.empty(len(local_dict), dtype=np.int32)
        for local_code, value in enumerate(local_dict):
            lut[local_code] = self._spec.null_code if value is None else self.encode(value)

        if ordered and len(local_dict) > 0:
            self._validate_ordered_chunk_dict(local_dict)

        # Translate Arrow indices to global int32 codes with a single numpy
        # gather (lut[indices]) instead of a per-row Python loop.
        indices = arrow_col.indices
        if indices.null_count:
            if not self._spec.nullable:
                raise ValueError("Dictionary column is not nullable but Arrow input contains nulls.")
            if len(lut) == 0:
                global_codes = np.full(m, self._spec.null_code, dtype=np.int32)
            else:
                null_mask = np.asarray(indices.is_null())
                local_codes = indices.fill_null(0).to_numpy(zero_copy_only=False)
                global_codes = lut[local_codes]
                global_codes[null_mask] = self._spec.null_code
        elif len(lut) == 0:
            # No local entries and no nulls means an empty batch.
            global_codes = np.empty(m, dtype=np.int32)
        else:
            local_codes = indices.to_numpy(zero_copy_only=False)
            global_codes = lut[local_codes]

        self._codes[pos : pos + m] = global_codes

    def _validate_ordered_chunk_dict(self, local_dict: list) -> None:
        """Raise if *local_dict* order differs from the existing global order."""
        self._ensure_cache()
        assert self._value_to_code is not None
        for local_code, value in enumerate(local_dict):
            if value is None:
                continue
            global_code = self._value_to_code.get(value)
            if global_code is not None and global_code != local_code:
                raise ValueError(
                    f"ordered=True dictionary column has inconsistent ordering across Arrow "
                    f"batches: value {value!r} has global code {global_code} but appears as "
                    f"local code {local_code} in this chunk."
                )

    # ------------------------------------------------------------------
    # Core interface: __len__, __getitem__, __setitem__
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        """Return the physical slot capacity (same as the codes NDArray length)."""
        return len(self._codes)

    def __getitem__(self, key) -> str | None | list:
        """Return decoded value(s) for the given index.

        - ``int``                → ``str | None``
        - ``slice``              → ``list``
        - ``numpy.ndarray``/``list`` → ``list``
        """
        if isinstance(key, (int, np.integer)):
            return self.decode(int(self._codes[int(key)]))
        if isinstance(key, slice):
            codes_arr = np.asarray(self._codes[key], dtype=np.int32)
            return [self.decode(int(c)) for c in codes_arr]
        if isinstance(key, (list, np.ndarray)):
            codes_arr = self._codes[key]
            if isinstance(codes_arr, np.ndarray):
                return [self.decode(int(c)) for c in codes_arr.ravel()]
            return [self.decode(int(codes_arr))]
        raise TypeError(f"DictionaryColumn indices must be int, slice, or array; got {type(key)!r}")

    def __setitem__(self, key, value) -> None:
        """Encode *value* (str/None or list thereof) and write the code(s)."""
        if isinstance(key, (int, np.integer)):
            self._codes[int(key)] = np.int32(self.encode(value))
        elif isinstance(key, slice):
            if isinstance(value, (list, tuple, np.ndarray)):
                self._codes[key] = self.encode_batch(list(value))
            else:
                # scalar broadcast
                code = np.int32(self.encode(value))
                self._codes[key] = code
        elif isinstance(key, (list, np.ndarray)):
            self._codes[key] = self.encode_batch(list(value))
        else:
            raise TypeError(f"DictionaryColumn indices must be int, slice, or array; got {type(key)!r}")

    def resize(self, shape: tuple) -> None:
        """Resize the underlying codes NDArray (delegates to the NDArray)."""
        self._codes.resize(shape)

    # ------------------------------------------------------------------
    # Flush / close
    # ------------------------------------------------------------------

    def flush(self) -> None:
        """Flush pending dict_store batches to the backend."""
        self._dict_store.flush()

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def codes(self):
        """The underlying ``int32`` NDArray of category codes."""
        return self._codes

    @property
    def dictionary(self) -> list[str]:
        """Return the list of unique dictionary values in insertion order."""
        self._dict_store.flush()
        return list(self._dict_store)

    @property
    def spec(self) -> DictionarySpec:
        return self._spec

    @property
    def dtype(self):
        """Always ``None`` — dictionary columns have no fixed NumPy dtype."""
        return None

    @property
    def urlpath(self) -> str | None:
        return getattr(self._codes, "urlpath", None)

    @property
    def nbytes(self) -> int:
        return self._codes.nbytes + self._dict_store.nbytes

    @property
    def cbytes(self) -> int:
        return self._codes.cbytes + self._dict_store.cbytes

    @property
    def cratio(self) -> float:
        cb = self.cbytes
        if cb == 0:
            return float("inf")
        return self.nbytes / cb
