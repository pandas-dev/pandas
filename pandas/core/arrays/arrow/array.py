from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    TypeVar,
)

import numpy as np

from pandas._typing import (
    TakeIndexer,
    npt,
)
from pandas.compat import (
    pa_version_under1p01,
    pa_version_under2p0,
    pa_version_under5p0,
    pa_version_under6p0,
)
from pandas.util._decorators import doc

from pandas.core.dtypes.common import (
    is_array_like,
    is_bool_dtype,
    is_integer,
    is_scalar,
)
from pandas.core.dtypes.missing import isna

from pandas.core.arrays.base import ExtensionArray
from pandas.core.indexers import (
    check_array_indexer,
    validate_indices,
)

if not pa_version_under1p01:
    import pyarrow as pa
    import pyarrow.compute as pc

    from pandas.core.arrays.arrow._arrow_utils import fallback_performancewarning

if TYPE_CHECKING:
    from pandas import Series

ArrowExtensionArrayT = TypeVar("ArrowExtensionArrayT", bound="ArrowExtensionArray")


class ArrowExtensionArray(ExtensionArray):
    """
    Base class for ExtensionArray backed by Arrow array.
    """

    _data: pa.ChunkedArray

    def __init__(self, values: pa.ChunkedArray) -> None:
        self._data = values

    def __arrow_array__(self, type=None):
        """Convert myself to a pyarrow Array or ChunkedArray."""
        return self._data

    def equals(self, other) -> bool:
        if not isinstance(other, ArrowExtensionArray):
            return False
        # I'm told that pyarrow makes __eq__ behave like pandas' equals;
        #  TODO: is this documented somewhere?
        return self._data == other._data

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        return self._data.nbytes

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self._data)

    def isna(self) -> npt.NDArray[np.bool_]:
        """
        Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        if pa_version_under2p0:
            return self._data.is_null().to_pandas().values
        else:
            return self._data.is_null().to_numpy()

    def copy(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Return a shallow copy of the array.

        Underlying ChunkedArray is immutable, so a deep copy is unnecessary.

        Returns
        -------
        type(self)
        """
        return type(self)(self._data)

    def dropna(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Return ArrowExtensionArray without NA values.

        Returns
        -------
        ArrowExtensionArray
        """
        if pa_version_under6p0:
            fallback_performancewarning(version="6")
            return super().dropna()
        else:
            return type(self)(pc.drop_null(self._data))

    @doc(ExtensionArray.factorize)
    def factorize(self, na_sentinel: int = -1) -> tuple[np.ndarray, ExtensionArray]:
        encoded = self._data.dictionary_encode()
        indices = pa.chunked_array(
            [c.indices for c in encoded.chunks], type=encoded.type.index_type
        ).to_pandas()
        if indices.dtype.kind == "f":
            indices[np.isnan(indices)] = na_sentinel
        indices = indices.astype(np.int64, copy=False)

        if encoded.num_chunks:
            uniques = type(self)(encoded.chunk(0).dictionary)
        else:
            uniques = type(self)(pa.array([], type=encoded.type.value_type))

        return indices.values, uniques

    def take(
        self,
        indices: TakeIndexer,
        allow_fill: bool = False,
        fill_value: Any = None,
    ):
        """
        Take elements from an array.

        Parameters
        ----------
        indices : sequence of int or one-dimensional np.ndarray of int
            Indices to be taken.
        allow_fill : bool, default False
            How to handle negative values in `indices`.

            * False: negative values in `indices` indicate positional indices
              from the right (the default). This is similar to
              :func:`numpy.take`.

            * True: negative values in `indices` indicate
              missing values. These values are set to `fill_value`. Any other
              other negative values raise a ``ValueError``.

        fill_value : any, optional
            Fill value to use for NA-indices when `allow_fill` is True.
            This may be ``None``, in which case the default NA value for
            the type, ``self.dtype.na_value``, is used.

            For many ExtensionArrays, there will be two representations of
            `fill_value`: a user-facing "boxed" scalar, and a low-level
            physical NA value. `fill_value` should be the user-facing version,
            and the implementation should handle translating that to the
            physical version for processing the take if necessary.

        Returns
        -------
        ExtensionArray

        Raises
        ------
        IndexError
            When the indices are out of bounds for the array.
        ValueError
            When `indices` contains negative values other than ``-1``
            and `allow_fill` is True.

        See Also
        --------
        numpy.take
        api.extensions.take

        Notes
        -----
        ExtensionArray.take is called by ``Series.__getitem__``, ``.loc``,
        ``iloc``, when `indices` is a sequence of values. Additionally,
        it's called by :meth:`Series.reindex`, or any other method
        that causes realignment, with a `fill_value`.
        """
        # TODO: Remove once we got rid of the (indices < 0) check
        if not is_array_like(indices):
            indices_array = np.asanyarray(indices)
        else:
            # error: Incompatible types in assignment (expression has type
            # "Sequence[int]", variable has type "ndarray")
            indices_array = indices  # type: ignore[assignment]

        if len(self._data) == 0 and (indices_array >= 0).any():
            raise IndexError("cannot do a non-empty take")
        if indices_array.size > 0 and indices_array.max() >= len(self._data):
            raise IndexError("out of bounds value in 'indices'.")

        if allow_fill:
            fill_mask = indices_array < 0
            if fill_mask.any():
                validate_indices(indices_array, len(self._data))
                # TODO(ARROW-9433): Treat negative indices as NULL
                indices_array = pa.array(indices_array, mask=fill_mask)
                result = self._data.take(indices_array)
                if isna(fill_value):
                    return type(self)(result)
                # TODO: ArrowNotImplementedError: Function fill_null has no
                # kernel matching input types (array[string], scalar[string])
                result = type(self)(result)
                result[fill_mask] = fill_value
                return result
                # return type(self)(pc.fill_null(result, pa.scalar(fill_value)))
            else:
                # Nothing to fill
                return type(self)(self._data.take(indices))
        else:  # allow_fill=False
            # TODO(ARROW-9432): Treat negative indices as indices from the right.
            if (indices_array < 0).any():
                # Don't modify in-place
                indices_array = np.copy(indices_array)
                indices_array[indices_array < 0] += len(self._data)
            return type(self)(self._data.take(indices_array))

    def unique(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Compute the ArrowExtensionArray of unique values.

        Returns
        -------
        ArrowExtensionArray
        """
        if pa_version_under2p0:
            fallback_performancewarning(version="2")
            return super().unique()
        else:
            return type(self)(pc.unique(self._data))

    def value_counts(self, dropna: bool = True) -> Series:
        """
        Return a Series containing counts of each unique value.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of missing values.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts
        """
        from pandas import (
            Index,
            Series,
        )

        vc = self._data.value_counts()

        values = vc.field(0)
        counts = vc.field(1)
        if dropna and self._data.null_count > 0:
            mask = values.is_valid()
            values = values.filter(mask)
            counts = counts.filter(mask)

        # No missing values so we can adhere to the interface and return a numpy array.
        counts = np.array(counts)

        index = Index(type(self)(values))

        return Series(counts, index=index).astype("Int64")

    @classmethod
    def _concat_same_type(
        cls: type[ArrowExtensionArrayT], to_concat
    ) -> ArrowExtensionArrayT:
        """
        Concatenate multiple ArrowExtensionArrays.

        Parameters
        ----------
        to_concat : sequence of ArrowExtensionArrays

        Returns
        -------
        ArrowExtensionArray
        """
        import pyarrow as pa

        chunks = [array for ea in to_concat for array in ea._data.iterchunks()]
        arr = pa.chunked_array(chunks)
        return cls(arr)

    def __setitem__(self, key: int | slice | np.ndarray, value: Any) -> None:
        """Set one or more values inplace.

        Parameters
        ----------
        key : int, ndarray, or slice
            When called from, e.g. ``Series.__setitem__``, ``key`` will be
            one of

            * scalar int
            * ndarray of integers.
            * boolean ndarray
            * slice object

        value : ExtensionDtype.type, Sequence[ExtensionDtype.type], or object
            value or values to be set of ``key``.

        Returns
        -------
        None
        """
        key = check_array_indexer(self, key)
        indices = self._indexing_key_to_indices(key)
        value = self._maybe_convert_setitem_value(value)

        argsort = np.argsort(indices)
        indices = indices[argsort]

        if is_scalar(value):
            value = np.broadcast_to(value, len(self))
        elif len(indices) != len(value):
            raise ValueError("Length of indexer and values mismatch")
        else:
            value = np.asarray(value)[argsort]

        self._data = self._set_via_chunk_iteration(indices=indices, value=value)

    def _indexing_key_to_indices(
        self, key: int | slice | np.ndarray
    ) -> npt.NDArray[np.intp]:
        """
        Convert indexing key for self into positional indices.

        Parameters
        ----------
        key : int | slice | np.ndarray

        Returns
        -------
        npt.NDArray[np.intp]
        """
        n = len(self)
        if isinstance(key, slice):
            indices = np.arange(n)[key]
        elif is_integer(key):
            indices = np.arange(n)[[key]]  # type: ignore[index]
        elif is_bool_dtype(key):
            key = np.asarray(key)
            if len(key) != n:
                raise ValueError("Length of indexer and values mismatch")
            indices = key.nonzero()[0]
        else:
            key = np.asarray(key)
            indices = np.arange(n)[key]
        return indices

    def _maybe_convert_setitem_value(self, value):
        """Maybe convert value to be pyarrow compatible."""
        raise NotImplementedError()

    def _set_via_chunk_iteration(
        self, indices: npt.NDArray[np.intp], value: npt.NDArray[Any]
    ) -> pa.ChunkedArray:
        """
        Loop through the array chunks and set the new values while
        leaving the chunking layout unchanged.

        Parameters
        ----------
        indices : npt.NDArray[np.intp]
            Position indices for the underlying ChunkedArray.

        value : ExtensionDtype.type, Sequence[ExtensionDtype.type], or object
            value or values to be set of ``key``.

        Notes
        -----
        Assumes that indices is sorted. Caller is responsible for sorting.
        """
        new_data = []
        stop = 0
        for chunk in self._data.iterchunks():
            start, stop = stop, stop + len(chunk)
            if len(indices) == 0 or stop <= indices[0]:
                new_data.append(chunk)
            else:
                n = int(np.searchsorted(indices, stop, side="left"))
                c_ind = indices[:n] - start
                indices = indices[n:]
                n = len(c_ind)
                c_value, value = value[:n], value[n:]
                new_data.append(self._replace_with_indices(chunk, c_ind, c_value))
        return pa.chunked_array(new_data)

    @classmethod
    def _replace_with_indices(
        cls,
        chunk: pa.Array,
        indices: npt.NDArray[np.intp],
        value: npt.NDArray[Any],
    ) -> pa.Array:
        """
        Replace items selected with a set of positional indices.

        Analogous to pyarrow.compute.replace_with_mask, except that replacement
        positions are identified via indices rather than a mask.

        Parameters
        ----------
        chunk : pa.Array
        indices : npt.NDArray[np.intp]
        value : npt.NDArray[Any]
            Replacement value(s).

        Returns
        -------
        pa.Array
        """
        n = len(indices)

        if n == 0:
            return chunk

        start, stop = indices[[0, -1]]

        if (stop - start) == (n - 1):
            # fast path for a contiguous set of indices
            arrays = [
                chunk[:start],
                pa.array(value, type=chunk.type),
                chunk[stop + 1 :],
            ]
            arrays = [arr for arr in arrays if len(arr)]
            if len(arrays) == 1:
                return arrays[0]
            return pa.concat_arrays(arrays)

        mask = np.zeros(len(chunk), dtype=np.bool_)
        mask[indices] = True

        if pa_version_under5p0:
            arr = chunk.to_numpy(zero_copy_only=False)
            arr[mask] = value
            return pa.array(arr, type=chunk.type)

        if isna(value).all():
            return pc.if_else(mask, None, chunk)

        return pc.replace_with_mask(chunk, mask, value)
