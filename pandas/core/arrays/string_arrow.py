from __future__ import annotations

from collections import abc
from distutils.version import LooseVersion
from typing import Any, Optional, Sequence, Tuple, Type, Union

import numpy as np

from pandas._libs import missing as libmissing
from pandas._typing import ArrayLike

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import register_extension_dtype

import pandas as pd
from pandas.api.types import (
    is_array_like,
    is_bool_dtype,
    is_integer,
    is_integer_dtype,
    is_scalar,
)
from pandas.core.arrays.base import ExtensionArray
from pandas.core.indexers import check_array_indexer

try:
    import pyarrow as pa
except ImportError:
    pa = None
else:
    try:
        import pyarrow.compute as pc
    except ImportError:
        pass


def _as_pandas_scalar(arrow_scalar: pa.Scalar) -> Optional[str]:
    scalar = arrow_scalar.as_py()
    if scalar is None:
        return libmissing.NA
    else:
        return scalar


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


class ArrowStringArray(ExtensionArray):
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

    def __init__(self, values):
        self._chk_pyarrow_available()
        if isinstance(values, pa.Array):
            self.data = pa.chunked_array([values])
        elif isinstance(values, pa.ChunkedArray):
            self.data = values
        else:
            raise ValueError(f"Unsupported type '{type(values)}' for ArrowStringArray")

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

        # TODO(ARROW-9407): Accept pd.NA in Arrow
        scalars_corrected = [None if pd.isna(x) else x for x in scalars]
        return cls(pa.array(scalars_corrected, type=pa.string()))

    @property
    def dtype(self) -> ArrowStringDtype:
        """
        An instance of 'ArrowStringDtype'.
        """
        return ArrowStringDtype()

    def __array__(self, *args, **kwargs) -> "np.ndarray":
        """Correctly construct numpy arrays when passed to `np.asarray()`."""
        return self.data.__array__(*args, **kwargs)

    def __arrow_array__(self, type=None):
        """Convert myself to a pyarrow Array or ChunkedArray."""
        return self.data

    @property
    def size(self) -> int:
        """
        Return the number of elements in this array.

        Returns
        -------
        size : int
        """
        return len(self.data)

    @property
    def shape(self) -> Tuple[int]:
        """Return the shape of the data."""
        # This may be patched by pandas to support pseudo-2D operations.
        return (len(self.data),)

    @property
    def ndim(self) -> int:
        """Return the number of dimensions of the underlying data."""
        return 1

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self.data)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

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

        if isinstance(item, abc.Iterable):
            if not is_array_like(item):
                item = np.array(item)
            if len(item) == 0:
                return type(self)(pa.chunked_array([], type=pa.string()))
            elif is_integer_dtype(item):
                return self.take(item)
            elif is_bool_dtype(item):
                return type(self)(self.data.filter(item))
            else:
                raise IndexError(
                    "Only integers, slices and integer or "
                    "boolean arrays are valid indices."
                )
        elif is_integer(item):
            if item < 0:
                item += len(self)
            if item >= len(self):
                raise IndexError("index out of bounds")

        value = self.data[item]
        if isinstance(value, pa.ChunkedArray):
            return type(self)(value)
        else:
            return _as_pandas_scalar(value)

    def fillna(self, value=None, method=None, limit=None):
        raise NotImplementedError("fillna")

    def _reduce(self, name, skipna=True, **kwargs):
        if name in ["min", "max"]:
            return getattr(self, name)(skipna=skipna)

        raise TypeError(f"Cannot perform reduction '{name}' with string dtype")

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        return self.data.nbytes

    def isna(self) -> np.ndarray:
        """
        Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        # TODO: Implement .to_numpy for ChunkedArray
        return self.data.is_null().to_pandas().values

    def copy(self) -> ExtensionArray:
        """
        Return a copy of the array.

        Parameters
        ----------
        deep : bool, default False
            Also copy the underlying data backing this array.

        Returns
        -------
        ExtensionArray
        """
        return type(self)(self.data)

    def __eq__(self, other: Any) -> ArrayLike:
        """
        Return for `self == other` (element-wise equality).
        """
        if isinstance(other, (pd.Series, pd.DataFrame, pd.Index)):
            return NotImplemented
        if isinstance(other, ArrowStringArray):
            result = pc.equal(self.data, other.data)
        elif is_scalar(other):
            result = pc.equal(self.data, pa.scalar(other))
        else:
            raise NotImplementedError("Neither scalar nor ArrowStringArray")

        # TODO(ARROW-9429): Add a .to_numpy() to ChunkedArray
        return pd.array(result.to_pandas().values)

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
            if not pd.api.types.is_scalar(value):
                raise ValueError("Must pass scalars with scalar indexer")
            elif pd.isna(value):
                value = None
            elif not isinstance(value, str):
                raise ValueError("Scalar must be NA or str")

            # Slice data and insert inbetween
            new_data = [
                *self.data[0:key].chunks,
                pa.array([value], type=pa.string()),
                *self.data[(key + 1) :].chunks,
            ]
            self.data = pa.chunked_array(new_data)
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

            if pd.api.types.is_scalar(value):
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

        if len(self.data) == 0 and (indices_array >= 0).any():
            raise IndexError("cannot do a non-empty take")
        if indices_array.max() >= len(self.data):
            raise IndexError("out of bounds value in 'indices'.")

        if allow_fill:
            if (indices_array < 0).any():
                # TODO(ARROW-9433): Treat negative indices as NULL
                indices_array = pa.array(indices_array, mask=indices_array < 0)
                result = self.data.take(indices_array)
                if pd.isna(fill_value):
                    return type(self)(result)
                return type(self)(pc.fill_null(result, pa.scalar(fill_value)))
            else:
                # Nothing to fill
                return type(self)(self.data.take(indices))
        else:  # allow_fill=False
            # TODO(ARROW-9432): Treat negative indices as indices from the right.
            if (indices_array < 0).any():
                # Don't modify in-place
                indices_array = np.copy(indices_array)
                indices_array[indices_array < 0] += len(self.data)
            return type(self)(self.data.take(indices_array))
