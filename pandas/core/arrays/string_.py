import operator
from typing import TYPE_CHECKING, Type, Union

import numpy as np

from pandas._libs import lib, missing as libmissing

from pandas.core.dtypes.base import ExtensionDtype, register_extension_dtype
from pandas.core.dtypes.common import pandas_dtype
from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries
from pandas.core.dtypes.inference import is_array_like

from pandas import compat
from pandas.core import ops
from pandas.core.arrays import IntegerArray, PandasArray
from pandas.core.arrays.integer import _IntegerDtype
from pandas.core.construction import extract_array
from pandas.core.indexers import check_array_indexer
from pandas.core.missing import isna

if TYPE_CHECKING:
    import pyarrow  # noqa: F401


@register_extension_dtype
class StringDtype(ExtensionDtype):
    """
    Extension dtype for string data.

    .. versionadded:: 1.0.0

    .. warning::

       StringDtype is considered experimental. The implementation and
       parts of the API may change without warning.

       In particular, StringDtype.na_value may change to no longer be
       ``numpy.nan``.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Examples
    --------
    >>> pd.StringDtype()
    StringDtype
    """

    name = "string"

    #: StringDtype.na_value uses pandas.NA
    na_value = libmissing.NA

    @property
    def type(self) -> Type[str]:
        return str

    @classmethod
    def construct_array_type(cls) -> Type["StringArray"]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return StringArray

    def __repr__(self) -> str:
        return "StringDtype"

    def __from_arrow__(
        self, array: Union["pyarrow.Array", "pyarrow.ChunkedArray"]
    ) -> "StringArray":
        """
        Construct StringArray from pyarrow Array/ChunkedArray.
        """
        import pyarrow  # noqa: F811

        if isinstance(array, pyarrow.Array):
            chunks = [array]
        else:
            # pyarrow.ChunkedArray
            chunks = array.chunks

        results = []
        for arr in chunks:
            # using _from_sequence to ensure None is converted to NA
            str_arr = StringArray._from_sequence(np.array(arr))
            results.append(str_arr)

        return StringArray._concat_same_type(results)


class StringArray(PandasArray):
    """
    Extension array for string data.

    .. versionadded:: 1.0.0

    .. warning::

       StringArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : array-like
        The array of data.

        .. warning::

           Currently, this expects an object-dtype ndarray
           where the elements are Python strings or :attr:`pandas.NA`.
           This may change without warning in the future. Use
           :meth:`pandas.array` with ``dtype="string"`` for a stable way of
           creating a `StringArray` from any sequence.

    copy : bool, default False
        Whether to copy the array of data.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    array
        The recommended function for creating a StringArray.
    Series.str
        The string methods are available on Series backed by
        a StringArray.

    Notes
    -----
    StringArray returns a BooleanArray for comparison methods.

    Examples
    --------
    >>> pd.array(['This is', 'some text', None, 'data.'], dtype="string")
    <StringArray>
    ['This is', 'some text', <NA>, 'data.']
    Length: 4, dtype: string

    Unlike arrays instantiated with ``dtype="object"``, ``StringArray``
    will convert the values to strings.

    >>> pd.array(['1', 1], dtype="object")
    <PandasArray>
    ['1', 1]
    Length: 2, dtype: object
    >>> pd.array(['1', 1], dtype="string")
    <StringArray>
    ['1', '1']
    Length: 2, dtype: string

    However, instantiating StringArrays directly with non-strings will raise an error.

    For comparison methods, `StringArray` returns a :class:`pandas.BooleanArray`:

    >>> pd.array(["a", None, "c"], dtype="string") == "a"
    <BooleanArray>
    [True, <NA>, False]
    Length: 3, dtype: boolean
    """

    # undo the PandasArray hack
    _typ = "extension"

    def __init__(self, values, copy=False):
        values = extract_array(values)

        super().__init__(values, copy=copy)
        self._dtype = StringDtype()
        if not isinstance(values, type(self)):
            self._validate()

    def _validate(self):
        """Validate that we only store NA or strings."""
        if len(self._ndarray) and not lib.is_string_array(self._ndarray, skipna=True):
            raise ValueError("StringArray requires a sequence of strings or pandas.NA")
        if self._ndarray.dtype != "object":
            raise ValueError(
                "StringArray requires a sequence of strings or pandas.NA. Got "
                f"'{self._ndarray.dtype}' dtype instead."
            )

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        if dtype:
            assert dtype == "string"

        result = np.asarray(scalars, dtype="object")

        # convert non-na-likes to str, and nan-likes to StringDtype.na_value
        result = lib.ensure_string_array(
            result, na_value=StringDtype.na_value, copy=copy
        )

        return cls(result)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

    def __arrow_array__(self, type=None):
        """
        Convert myself into a pyarrow Array.
        """
        import pyarrow as pa

        if type is None:
            type = pa.string()

        values = self._ndarray.copy()
        values[self.isna()] = None
        return pa.array(values, type=type, from_pandas=True)

    def _values_for_factorize(self):
        arr = self._ndarray.copy()
        mask = self.isna()
        arr[mask] = -1
        return arr, -1

    def __setitem__(self, key, value):
        value = extract_array(value, extract_numpy=True)
        if isinstance(value, type(self)):
            # extract_array doesn't extract PandasArray subclasses
            value = value._ndarray

        key = check_array_indexer(self, key)
        scalar_key = lib.is_scalar(key)
        scalar_value = lib.is_scalar(value)
        if scalar_key and not scalar_value:
            raise ValueError("setting an array element with a sequence.")

        # validate new items
        if scalar_value:
            if isna(value):
                value = StringDtype.na_value
            elif not isinstance(value, str):
                raise ValueError(
                    f"Cannot set non-string value '{value}' into a StringArray."
                )
        else:
            if not is_array_like(value):
                value = np.asarray(value, dtype=object)
            if len(value) and not lib.is_string_array(value, skipna=True):
                raise ValueError("Must provide strings.")

        super().__setitem__(key, value)

    def fillna(self, value=None, method=None, limit=None):
        # TODO: validate dtype
        return super().fillna(value, method, limit)

    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if isinstance(dtype, StringDtype):
            if copy:
                return self.copy()
            return self
        elif isinstance(dtype, _IntegerDtype):
            arr = self._ndarray.copy()
            mask = self.isna()
            arr[mask] = 0
            values = arr.astype(dtype.numpy_dtype)
            return IntegerArray(values, mask, copy=False)

        return super().astype(dtype, copy)

    def _reduce(self, name: str, skipna: bool = True, **kwargs):
        if name in ["min", "max"]:
            return getattr(self, name)(skipna=skipna)

        raise TypeError(f"Cannot perform reduction '{name}' with string dtype")

    def value_counts(self, dropna=False):
        from pandas import value_counts

        return value_counts(self._ndarray, dropna=dropna).astype("Int64")

    def memory_usage(self, deep=False):
        result = self._ndarray.nbytes
        if deep:
            return result + lib.memory_usage_of_objects(self._ndarray)
        return result

    # Override parent because we have different return types.
    @classmethod
    def _create_arithmetic_method(cls, op):
        # Note: this handles both arithmetic and comparison methods.
        def method(self, other):
            from pandas.arrays import BooleanArray

            assert op.__name__ in ops.ARITHMETIC_BINOPS | ops.COMPARISON_BINOPS

            if isinstance(other, (ABCIndexClass, ABCSeries, ABCDataFrame)):
                return NotImplemented

            elif isinstance(other, cls):
                other = other._ndarray

            mask = isna(self) | isna(other)
            valid = ~mask

            if not lib.is_scalar(other):
                if len(other) != len(self):
                    # prevent improper broadcasting when other is 2D
                    raise ValueError(
                        f"Lengths of operands do not match: {len(self)} != {len(other)}"
                    )

                other = np.asarray(other)
                other = other[valid]

            if op.__name__ in ops.ARITHMETIC_BINOPS:
                result = np.empty_like(self._ndarray, dtype="object")
                result[mask] = StringDtype.na_value
                result[valid] = op(self._ndarray[valid], other)
                return StringArray(result)
            else:
                # logical
                result = np.zeros(len(self._ndarray), dtype="bool")
                result[valid] = op(self._ndarray[valid], other)
                return BooleanArray(result, mask)

        return compat.set_function_name(method, f"__{op.__name__}__", cls)

    @classmethod
    def _add_arithmetic_ops(cls):
        cls.__add__ = cls._create_arithmetic_method(operator.add)
        cls.__radd__ = cls._create_arithmetic_method(ops.radd)

        cls.__mul__ = cls._create_arithmetic_method(operator.mul)
        cls.__rmul__ = cls._create_arithmetic_method(ops.rmul)

    _create_comparison_method = _create_arithmetic_method


StringArray._add_arithmetic_ops()
StringArray._add_comparison_ops()
