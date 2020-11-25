from __future__ import annotations

from distutils.version import LooseVersion
from typing import TYPE_CHECKING, Any, Sequence, Type, Union

import numpy as np

from pandas._libs import lib, missing as libmissing
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import register_extension_dtype
from pandas.core.dtypes.missing import isna

from pandas.api.types import (
    is_array_like,
    is_bool_dtype,
    is_integer,
    is_integer_dtype,
    is_scalar,
)
from pandas.core.arraylike import OpsMixin
from pandas.core.arrays.base import ExtensionArray
from pandas.core.indexers import check_array_indexer, validate_indices
from pandas.core.missing import get_fill_func

try:
    import pyarrow as pa
except ImportError:
    pa = None
else:
    # our min supported version of pyarrow, 0.15.1, does not have a compute
    # module
    try:
        import pyarrow.compute as pc
    except ImportError:
        pass
    else:
        ARROW_CMP_FUNCS = {
            "eq": pc.equal,
            "ne": pc.not_equal,
            "lt": pc.less,
            "gt": pc.greater,
            "le": pc.less_equal,
            "ge": pc.greater_equal,
        }


if TYPE_CHECKING:
    from pandas import Series


@register_extension_dtype
class ArrowStringDtype(ExtensionDtype):
    """
    Extension dtype for string data in a ``pyarrow.ChunkedArray``.

    .. versionadded:: 1.2.0

    .. warning::

       ArrowStringDtype is considered experimental. The implementation and
       parts of the API may change without warning.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Examples
    --------
    >>> from pandas.core.arrays.string_arrow import ArrowStringDtype
    >>> ArrowStringDtype()
    ArrowStringDtype
    """

    name = "arrow_string"

    #: StringDtype.na_value uses pandas.NA
    na_value = libmissing.NA

    @property
    def type(self) -> Type[str]:
        return str

    @classmethod
    def construct_array_type(cls) -> Type["ArrowStringArray"]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return ArrowStringArray

    def __hash__(self) -> int:
        return hash("ArrowStringDtype")

    def __repr__(self) -> str:
        return "ArrowStringDtype"

    def __from_arrow__(
        self, array: Union["pa.Array", "pa.ChunkedArray"]
    ) -> "ArrowStringArray":
        """
        Construct StringArray from pyarrow Array/ChunkedArray.
        """
        return ArrowStringArray(array)

    def __eq__(self, other) -> bool:
        """Check whether 'other' is equal to self.

        By default, 'other' is considered equal if
        * it's a string matching 'self.name'.
        * it's an instance of this type.

        Parameters
        ----------
        other : Any

        Returns
        -------
        bool
        """
        if isinstance(other, ArrowStringDtype):
            return True
        elif isinstance(other, str) and other == "arrow_string":
            return True
        else:
            return False


class ArrowStringArray(OpsMixin, ExtensionArray):
    """
    Extension array for string data in a ``pyarrow.ChunkedArray``.

    .. versionadded:: 1.2.0

    .. warning::

       ArrowStringArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : pyarrow.Array or pyarrow.ChunkedArray
        The array of data.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    array
        The recommended function for creating a ArrowStringArray.
    Series.str
        The string methods are available on Series backed by
        a ArrowStringArray.

    Notes
    -----
    ArrowStringArray returns a BooleanArray for comparison methods.

    Examples
    --------
    >>> pd.array(['This is', 'some text', None, 'data.'], dtype="arrow_string")
    <ArrowStringArray>
    ['This is', 'some text', <NA>, 'data.']
    Length: 4, dtype: arrow_string
    """

    _dtype = ArrowStringDtype()

    def __init__(self, values):
        self._chk_pyarrow_available()
        if isinstance(values, pa.Array):
            self._data = pa.chunked_array([values])
        elif isinstance(values, pa.ChunkedArray):
            self._data = values
        else:
            raise ValueError(f"Unsupported type '{type(values)}' for ArrowStringArray")

        if not pa.types.is_string(self._data.type):
            raise ValueError(
                "ArrowStringArray requires a PyArrow (chunked) array of string type"
            )

    @classmethod
    def _chk_pyarrow_available(cls) -> None:
        # TODO: maybe update import_optional_dependency to allow a minimum
        # version to be specified rather than use the global minimum
        if pa is None or LooseVersion(pa.__version__) < "1.0.0":
            msg = "pyarrow>=1.0.0 is required for PyArrow backed StringArray."
            raise ImportError(msg)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        cls._chk_pyarrow_available()
        # convert non-na-likes to str, and nan-likes to ArrowStringDtype.na_value
        scalars = lib.ensure_string_array(scalars, copy=False)
        return cls(pa.array(scalars, type=pa.string(), from_pandas=True))

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

    @property
    def dtype(self) -> ArrowStringDtype:
        """
        An instance of 'ArrowStringDtype'.
        """
        return self._dtype

    def __array__(self, dtype=None) -> np.ndarray:
        """Correctly construct numpy arrays when passed to `np.asarray()`."""
        return self.to_numpy(dtype=dtype)

    def __arrow_array__(self, type=None):
        """Convert myself to a pyarrow Array or ChunkedArray."""
        return self._data

    def to_numpy(
        self, dtype=None, copy: bool = False, na_value=lib.no_default
    ) -> np.ndarray:
        """
        Convert to a NumPy ndarray.
        """
        # TODO: copy argument is ignored

        if na_value is lib.no_default:
            na_value = self._dtype.na_value
        result = self._data.__array__(dtype=dtype)
        result[isna(result)] = na_value
        return result

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self._data)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls._from_sequence(values)

    @classmethod
    def _concat_same_type(cls, to_concat) -> ArrowStringArray:
        """
        Concatenate multiple ArrowStringArray.

        Parameters
        ----------
        to_concat : sequence of ArrowStringArray

        Returns
        -------
        ArrowStringArray
        """
        return cls(
            pa.chunked_array(
                [array for ea in to_concat for array in ea._data.iterchunks()]
            )
        )

    def __getitem__(self, item: Any) -> Any:
        """Select a subset of self.

        Parameters
        ----------
        item : int, slice, or ndarray
            * int: The position in 'self' to get.
            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None
            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : scalar or ExtensionArray

        Notes
        -----
        For scalar ``item``, return a scalar value suitable for the array's
        type. This should be an instance of ``self.dtype.type``.
        For slice ``key``, return an instance of ``ExtensionArray``, even
        if the slice is length 0 or 1.
        For a boolean mask, return an instance of ``ExtensionArray``, filtered
        to the values where ``item`` is True.
        """
        item = check_array_indexer(self, item)

        if isinstance(item, np.ndarray):
            if not len(item):
                return type(self)(pa.chunked_array([], type=pa.string()))
            elif is_integer_dtype(item.dtype):
                return self.take(item)
            elif is_bool_dtype(item.dtype):
                return type(self)(self._data.filter(item))
            else:
                raise IndexError(
                    "Only integers, slices and integer or "
                    "boolean arrays are valid indices."
                )

        # We are not an array indexer, so maybe e.g. a slice or integer
        # indexer. We dispatch to pyarrow.
        value = self._data[item]
        if isinstance(value, pa.ChunkedArray):
            return type(self)(value)
        else:
            return self._as_pandas_scalar(value)

    def _as_pandas_scalar(self, arrow_scalar: pa.Scalar):
        scalar = arrow_scalar.as_py()
        if scalar is None:
            return self._dtype.na_value
        else:
            return scalar

    def fillna(self, value=None, method=None, limit=None):
        """
        Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, array-like
            If a scalar value is passed it is used to fill all missing values.
            Alternatively, an array-like 'value' can be given. It's expected
            that the array-like have the same length as 'self'.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap.
        limit : int, default None
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        ExtensionArray
            With NA/NaN filled.
        """
        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError(
                    f"Length of 'value' does not match. Got ({len(value)}) "
                    f"expected {len(self)}"
                )
            value = value[mask]

        if mask.any():
            if method is not None:
                func = get_fill_func(method)
                new_values = func(self.to_numpy(object), limit=limit, mask=mask)
                new_values = self._from_sequence(new_values)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def _reduce(self, name, skipna=True, **kwargs):
        if name in ["min", "max"]:
            return getattr(self, name)(skipna=skipna)

        raise TypeError(f"Cannot perform reduction '{name}' with string dtype")

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        return self._data.nbytes

    def isna(self) -> np.ndarray:
        """
        Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        # TODO: Implement .to_numpy for ChunkedArray
        return self._data.is_null().to_pandas().values

    def copy(self) -> ArrowStringArray:
        """
        Return a shallow copy of the array.

        Returns
        -------
        ArrowStringArray
        """
        return type(self)(self._data)

    def _cmp_method(self, other, op):
        from pandas.arrays import BooleanArray

        pc_func = ARROW_CMP_FUNCS[op.__name__]
        if isinstance(other, ArrowStringArray):
            result = pc_func(self._data, other._data)
        elif isinstance(other, np.ndarray):
            result = pc_func(self._data, other)
        elif is_scalar(other):
            try:
                result = pc_func(self._data, pa.scalar(other))
            except (pa.lib.ArrowNotImplementedError, pa.lib.ArrowInvalid):
                mask = isna(self) | isna(other)
                valid = ~mask
                result = np.zeros(len(self), dtype="bool")
                result[valid] = op(np.array(self)[valid], other)
                return BooleanArray(result, mask)
        else:
            return NotImplemented

        # TODO(ARROW-9429): Add a .to_numpy() to ChunkedArray
        return BooleanArray._from_sequence(result.to_pandas().values)

    def __setitem__(self, key: Union[int, np.ndarray], value: Any) -> None:
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

        if is_integer(key):
            if not is_scalar(value):
                raise ValueError("Must pass scalars with scalar indexer")
            elif isna(value):
                value = None
            elif not isinstance(value, str):
                raise ValueError("Scalar must be NA or str")

            # Slice data and insert inbetween
            new_data = [
                *self._data[0:key].chunks,
                pa.array([value], type=pa.string()),
                *self._data[(key + 1) :].chunks,
            ]
            self._data = pa.chunked_array(new_data)
        else:
            # Convert to integer indices and iteratively assign.
            # TODO: Make a faster variant of this in Arrow upstream.
            #       This is probably extremely slow.

            # Convert all possible input key types to an array of integers
            if is_bool_dtype(key):
                # TODO(ARROW-9430): Directly support setitem(booleans)
                key_array = np.argwhere(key).flatten()
            elif isinstance(key, slice):
                key_array = np.array(range(len(self))[key])
            else:
                # TODO(ARROW-9431): Directly support setitem(integers)
                key_array = np.asanyarray(key)

            if is_scalar(value):
                value = np.broadcast_to(value, len(key_array))
            else:
                value = np.asarray(value)

            if len(key_array) != len(value):
                raise ValueError("Length of indexer and values mismatch")

            for k, v in zip(key_array, value):
                self[k] = v

    def take(
        self, indices: Sequence[int], allow_fill: bool = False, fill_value: Any = None
    ) -> "ExtensionArray":
        """
        Take elements from an array.

        Parameters
        ----------
        indices : sequence of int
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
            indices_array = indices

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
        from pandas import Index, Series

        vc = self._data.value_counts()

        # Index cannot hold ExtensionArrays yet
        index = Index(type(self)(vc.field(0)).astype(object))
        # No missings, so we can adhere to the interface and return a numpy array.
        counts = np.array(vc.field(1))

        if dropna and self._data.null_count > 0:
            raise NotImplementedError("yo")

        return Series(counts, index=index).astype("Int64")
