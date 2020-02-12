import numbers
from typing import Any, Tuple, Type
import warnings

import numpy as np

from pandas._libs import lib, missing as libmissing
from pandas.compat import set_function_name
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.cast import astype_nansafe
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import register_extension_dtype
from pandas.core.dtypes.missing import isna

from pandas.core import nanops, ops
from pandas.core.indexers import check_array_indexer
from pandas.core.ops import invalid_comparison
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.tools.numeric import to_numeric

from .masked import BaseMaskedArray


class _IntegerDtype(ExtensionDtype):
    """
    An ExtensionDtype to hold a single size & kind of integer dtype.

    These specific implementations are subclasses of the non-public
    _IntegerDtype. For example we have Int8Dtype to represent signed int 8s.

    The attributes name & type are set when these subclasses are created.
    """

    name: str
    base = None
    type: Type
    na_value = libmissing.NA

    def __repr__(self) -> str:
        sign = "U" if self.is_unsigned_integer else ""
        return f"{sign}Int{8 * self.itemsize}Dtype()"

    @cache_readonly
    def is_signed_integer(self):
        return self.kind == "i"

    @cache_readonly
    def is_unsigned_integer(self):
        return self.kind == "u"

    @property
    def _is_numeric(self):
        return True

    @cache_readonly
    def numpy_dtype(self):
        """ Return an instance of our numpy dtype """
        return np.dtype(self.type)

    @cache_readonly
    def kind(self):
        return self.numpy_dtype.kind

    @cache_readonly
    def itemsize(self):
        """ Return the number of bytes in this dtype """
        return self.numpy_dtype.itemsize

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return IntegerArray

    def __from_arrow__(self, array):
        """Construct IntegerArray from passed pyarrow Array/ChunkedArray"""
        import pyarrow
        from pandas.core.arrays._arrow_utils import pyarrow_array_to_numpy_and_mask

        pyarrow_type = pyarrow.from_numpy_dtype(self.type)
        if not array.type.equals(pyarrow_type):
            array = array.cast(pyarrow_type)

        if isinstance(array, pyarrow.Array):
            chunks = [array]
        else:
            # pyarrow.ChunkedArray
            chunks = array.chunks

        results = []
        for arr in chunks:
            data, mask = pyarrow_array_to_numpy_and_mask(arr, dtype=self.type)
            int_arr = IntegerArray(data.copy(), ~mask, copy=False)
            results.append(int_arr)

        return IntegerArray._concat_same_type(results)


def integer_array(values, dtype=None, copy=False):
    """
    Infer and return an integer array of the values.

    Parameters
    ----------
    values : 1D list-like
    dtype : dtype, optional
        dtype to coerce
    copy : bool, default False

    Returns
    -------
    IntegerArray

    Raises
    ------
    TypeError if incompatible types
    """
    values, mask = coerce_to_array(values, dtype=dtype, copy=copy)
    return IntegerArray(values, mask)


def safe_cast(values, dtype, copy):
    """
    Safely cast the values to the dtype if they
    are equivalent, meaning floats must be equivalent to the
    ints.

    """

    try:
        return values.astype(dtype, casting="safe", copy=copy)
    except TypeError:

        casted = values.astype(dtype, copy=copy)
        if (casted == values).all():
            return casted

        raise TypeError(
            f"cannot safely cast non-equivalent {values.dtype} to {np.dtype(dtype)}"
        )


def coerce_to_array(values, dtype, mask=None, copy=False):
    """
    Coerce the input values array to numpy arrays with a mask

    Parameters
    ----------
    values : 1D list-like
    dtype : integer dtype
    mask : bool 1D array, optional
    copy : bool, default False
        if True, copy the input

    Returns
    -------
    tuple of (values, mask)
    """
    # if values is integer numpy array, preserve it's dtype
    if dtype is None and hasattr(values, "dtype"):
        if is_integer_dtype(values.dtype):
            dtype = values.dtype

    if dtype is not None:
        if isinstance(dtype, str) and (
            dtype.startswith("Int") or dtype.startswith("UInt")
        ):
            # Avoid DeprecationWarning from NumPy about np.dtype("Int64")
            # https://github.com/numpy/numpy/pull/7476
            dtype = dtype.lower()

        if not issubclass(type(dtype), _IntegerDtype):
            try:
                dtype = _dtypes[str(np.dtype(dtype))]
            except KeyError:
                raise ValueError(f"invalid dtype specified {dtype}")

    if isinstance(values, IntegerArray):
        values, mask = values._data, values._mask
        if dtype is not None:
            values = values.astype(dtype.numpy_dtype, copy=False)

        if copy:
            values = values.copy()
            mask = mask.copy()
        return values, mask

    values = np.array(values, copy=copy)
    if is_object_dtype(values):
        inferred_type = lib.infer_dtype(values, skipna=True)
        if inferred_type == "empty":
            values = np.empty(len(values))
            values.fill(np.nan)
        elif inferred_type not in [
            "floating",
            "integer",
            "mixed-integer",
            "integer-na",
            "mixed-integer-float",
        ]:
            raise TypeError(f"{values.dtype} cannot be converted to an IntegerDtype")

    elif is_bool_dtype(values) and is_integer_dtype(dtype):
        values = np.array(values, dtype=int, copy=copy)

    elif not (is_integer_dtype(values) or is_float_dtype(values)):
        raise TypeError(f"{values.dtype} cannot be converted to an IntegerDtype")

    if mask is None:
        mask = isna(values)
    else:
        assert len(mask) == len(values)

    if not values.ndim == 1:
        raise TypeError("values must be a 1D list-like")
    if not mask.ndim == 1:
        raise TypeError("mask must be a 1D list-like")

    # infer dtype if needed
    if dtype is None:
        dtype = np.dtype("int64")
    else:
        dtype = dtype.type

    # if we are float, let's make sure that we can
    # safely cast

    # we copy as need to coerce here
    if mask.any():
        values = values.copy()
        values[mask] = 1
        values = safe_cast(values, dtype, copy=False)
    else:
        values = safe_cast(values, dtype, copy=False)

    return values, mask


class IntegerArray(BaseMaskedArray):
    """
    Array of integer (optional missing) values.

    .. versionadded:: 0.24.0

    .. versionchanged:: 1.0.0

       Now uses :attr:`pandas.NA` as the missing value rather
       than :attr:`numpy.nan`.

    .. warning::

       IntegerArray is currently experimental, and its API or internal
       implementation may change without warning.

    We represent an IntegerArray with 2 numpy arrays:

    - data: contains a numpy integer array of the appropriate dtype
    - mask: a boolean array holding a mask on the data, True is missing

    To construct an IntegerArray from generic array-like input, use
    :func:`pandas.array` with one of the integer dtypes (see examples).

    See :ref:`integer_na` for more.

    Parameters
    ----------
    values : numpy.ndarray
        A 1-d integer-dtype array.
    mask : numpy.ndarray
        A 1-d boolean-dtype array indicating missing values.
    copy : bool, default False
        Whether to copy the `values` and `mask`.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Returns
    -------
    IntegerArray

    Examples
    --------
    Create an IntegerArray with :func:`pandas.array`.

    >>> int_array = pd.array([1, None, 3], dtype=pd.Int32Dtype())
    >>> int_array
    <IntegerArray>
    [1, <NA>, 3]
    Length: 3, dtype: Int32

    String aliases for the dtypes are also available. They are capitalized.

    >>> pd.array([1, None, 3], dtype='Int32')
    <IntegerArray>
    [1, <NA>, 3]
    Length: 3, dtype: Int32

    >>> pd.array([1, None, 3], dtype='UInt16')
    <IntegerArray>
    [1, <NA>, 3]
    Length: 3, dtype: UInt16
    """

    # The value used to fill '_data' to avoid upcasting
    _internal_fill_value = 1

    @cache_readonly
    def dtype(self):
        return _dtypes[str(self._data.dtype)]

    def __init__(self, values, mask, copy=False):
        if not (isinstance(values, np.ndarray) and is_integer_dtype(values.dtype)):
            raise TypeError(
                "values should be integer numpy array. Use "
                "the 'integer_array' function instead"
            )
        if not (isinstance(mask, np.ndarray) and is_bool_dtype(mask.dtype)):
            raise TypeError(
                "mask should be boolean numpy array. Use "
                "the 'integer_array' function instead"
            )

        if copy:
            values = values.copy()
            mask = mask.copy()

        self._data = values
        self._mask = mask

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return integer_array(scalars, dtype=dtype, copy=copy)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        scalars = to_numeric(strings, errors="raise")
        return cls._from_sequence(scalars, dtype, copy)

    @classmethod
    def _from_factorized(cls, values, original):
        return integer_array(values, dtype=original.dtype)

    _HANDLED_TYPES = (np.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # For IntegerArray inputs, we apply the ufunc to ._data
        # and mask the result.
        if method == "reduce":
            # Not clear how to handle missing values in reductions. Raise.
            raise NotImplementedError("The 'reduce' method is not supported.")
        out = kwargs.get("out", ())

        for x in inputs + out:
            if not isinstance(x, self._HANDLED_TYPES + (IntegerArray,)):
                return NotImplemented

        # for binary ops, use our custom dunder methods
        result = ops.maybe_dispatch_ufunc_to_dunder_op(
            self, ufunc, method, *inputs, **kwargs
        )
        if result is not NotImplemented:
            return result

        mask = np.zeros(len(self), dtype=bool)
        inputs2 = []
        for x in inputs:
            if isinstance(x, IntegerArray):
                mask |= x._mask
                inputs2.append(x._data)
            else:
                inputs2.append(x)

        def reconstruct(x):
            # we don't worry about scalar `x` here, since we
            # raise for reduce up above.

            if is_integer_dtype(x.dtype):
                m = mask.copy()
                return IntegerArray(x, m)
            else:
                x[mask] = np.nan
            return x

        result = getattr(ufunc, method)(*inputs2, **kwargs)
        if isinstance(result, tuple):
            tuple(reconstruct(x) for x in result)
        else:
            return reconstruct(result)

    def __setitem__(self, key, value):
        _is_scalar = is_scalar(value)
        if _is_scalar:
            value = [value]
        value, mask = coerce_to_array(value, dtype=self.dtype)

        if _is_scalar:
            value = value[0]
            mask = mask[0]

        key = check_array_indexer(self, key)
        self._data[key] = value
        self._mask[key] = mask

    def astype(self, dtype, copy=True):
        """
        Cast to a NumPy array or IntegerArray with 'dtype'.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        copy : bool, default True
            Whether to copy the data, even if not necessary. If False,
            a copy is made only if the old dtype does not match the
            new dtype.

        Returns
        -------
        array : ndarray or IntegerArray
            NumPy ndarray or IntergerArray with 'dtype' for its dtype.

        Raises
        ------
        TypeError
            if incompatible type with an IntegerDtype, equivalent of same_kind
            casting
        """
        from pandas.core.arrays.boolean import BooleanArray, BooleanDtype

        dtype = pandas_dtype(dtype)

        # if we are astyping to an existing IntegerDtype we can fastpath
        if isinstance(dtype, _IntegerDtype):
            result = self._data.astype(dtype.numpy_dtype, copy=False)
            return type(self)(result, mask=self._mask, copy=False)
        elif isinstance(dtype, BooleanDtype):
            result = self._data.astype("bool", copy=False)
            return BooleanArray(result, mask=self._mask, copy=False)

        # coerce
        if is_float_dtype(dtype):
            # In astype, we consider dtype=float to also mean na_value=np.nan
            kwargs = dict(na_value=np.nan)
        else:
            kwargs = {}

        data = self.to_numpy(dtype=dtype, **kwargs)
        return astype_nansafe(data, dtype, copy=False)

    @property
    def _ndarray_values(self) -> np.ndarray:
        """Internal pandas method for lossy conversion to a NumPy ndarray.

        This method is not part of the pandas interface.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.
        """
        return self._data

    def _values_for_factorize(self) -> Tuple[np.ndarray, Any]:
        # TODO: https://github.com/pandas-dev/pandas/issues/30037
        # use masked algorithms, rather than object-dtype / np.nan.
        return self.to_numpy(na_value=np.nan), np.nan

    def _values_for_argsort(self) -> np.ndarray:
        """Return values for sorting.

        Returns
        -------
        ndarray
            The transformed values should maintain the ordering between values
            within the array.

        See Also
        --------
        ExtensionArray.argsort
        """
        data = self._data.copy()
        data[self._mask] = data.min() - 1
        return data

    @classmethod
    def _create_comparison_method(cls, op):
        op_name = op.__name__

        @unpack_zerodim_and_defer(op.__name__)
        def cmp_method(self, other):
            from pandas.arrays import BooleanArray

            mask = None

            if isinstance(other, (BooleanArray, IntegerArray)):
                other, mask = other._data, other._mask

            elif is_list_like(other):
                other = np.asarray(other)
                if other.ndim > 1:
                    raise NotImplementedError(
                        "can only perform ops with 1-d structures"
                    )
                if len(self) != len(other):
                    raise ValueError("Lengths must match to compare")

            if other is libmissing.NA:
                # numpy does not handle pd.NA well as "other" scalar (it returns
                # a scalar False instead of an array)
                # This may be fixed by NA.__array_ufunc__. Revisit this check
                # once that's implemented.
                result = np.zeros(self._data.shape, dtype="bool")
                mask = np.ones(self._data.shape, dtype="bool")
            else:
                with warnings.catch_warnings():
                    # numpy may show a FutureWarning:
                    #     elementwise comparison failed; returning scalar instead,
                    #     but in the future will perform elementwise comparison
                    # before returning NotImplemented. We fall back to the correct
                    # behavior today, so that should be fine to ignore.
                    warnings.filterwarnings("ignore", "elementwise", FutureWarning)
                    with np.errstate(all="ignore"):
                        method = getattr(self._data, f"__{op_name}__")
                        result = method(other)

                    if result is NotImplemented:
                        result = invalid_comparison(self._data, other, op)

            # nans propagate
            if mask is None:
                mask = self._mask.copy()
            else:
                mask = self._mask | mask

            return BooleanArray(result, mask)

        name = f"__{op.__name__}__"
        return set_function_name(cmp_method, name, cls)

    def _reduce(self, name, skipna=True, **kwargs):
        data = self._data
        mask = self._mask

        # coerce to a nan-aware float if needed
        # (we explicitly use NaN within reductions)
        if self._hasna:
            data = self.to_numpy("float64", na_value=np.nan)

        op = getattr(nanops, "nan" + name)
        result = op(data, axis=0, skipna=skipna, mask=mask, **kwargs)

        if np.isnan(result):
            return libmissing.NA

        # if we have a boolean op, don't coerce
        if name in ["any", "all"]:
            pass

        # if we have a preservable numeric op,
        # provide coercion back to an integer type if possible
        elif name in ["sum", "min", "max", "prod"]:
            int_result = int(result)
            if int_result == result:
                result = int_result

        return result

    def _maybe_mask_result(self, result, mask, other, op_name):
        """
        Parameters
        ----------
        result : array-like
        mask : array-like bool
        other : scalar or array-like
        op_name : str
        """

        # if we have a float operand we are by-definition
        # a float result
        # or our op is a divide
        if (is_float_dtype(other) or is_float(other)) or (
            op_name in ["rtruediv", "truediv"]
        ):
            result[mask] = np.nan
            return result

        return type(self)(result, mask, copy=False)

    @classmethod
    def _create_arithmetic_method(cls, op):
        op_name = op.__name__

        @unpack_zerodim_and_defer(op.__name__)
        def integer_arithmetic_method(self, other):

            omask = None

            if getattr(other, "ndim", 0) > 1:
                raise NotImplementedError("can only perform ops with 1-d structures")

            if isinstance(other, IntegerArray):
                other, omask = other._data, other._mask

            elif is_list_like(other):
                other = np.asarray(other)
                if other.ndim > 1:
                    raise NotImplementedError(
                        "can only perform ops with 1-d structures"
                    )
                if len(self) != len(other):
                    raise ValueError("Lengths must match")
                if not (is_float_dtype(other) or is_integer_dtype(other)):
                    raise TypeError("can only perform ops with numeric values")

            else:
                if not (is_float(other) or is_integer(other) or other is libmissing.NA):
                    raise TypeError("can only perform ops with numeric values")

            if omask is None:
                mask = self._mask.copy()
                if other is libmissing.NA:
                    mask |= True
            else:
                mask = self._mask | omask

            if op_name == "pow":
                # 1 ** x is 1.
                mask = np.where((self._data == 1) & ~self._mask, False, mask)
                # x ** 0 is 1.
                if omask is not None:
                    mask = np.where((other == 0) & ~omask, False, mask)
                elif other is not libmissing.NA:
                    mask = np.where(other == 0, False, mask)

            elif op_name == "rpow":
                # 1 ** x is 1.
                if omask is not None:
                    mask = np.where((other == 1) & ~omask, False, mask)
                elif other is not libmissing.NA:
                    mask = np.where(other == 1, False, mask)
                # x ** 0 is 1.
                mask = np.where((self._data == 0) & ~self._mask, False, mask)

            if other is libmissing.NA:
                result = np.ones_like(self._data)
            else:
                with np.errstate(all="ignore"):
                    result = op(self._data, other)

            # divmod returns a tuple
            if op_name == "divmod":
                div, mod = result
                return (
                    self._maybe_mask_result(div, mask, other, "floordiv"),
                    self._maybe_mask_result(mod, mask, other, "mod"),
                )

            return self._maybe_mask_result(result, mask, other, op_name)

        name = f"__{op.__name__}__"
        return set_function_name(integer_arithmetic_method, name, cls)


IntegerArray._add_arithmetic_ops()
IntegerArray._add_comparison_ops()


_dtype_docstring = """
An ExtensionDtype for {dtype} integer data.

.. versionchanged:: 1.0.0

   Now uses :attr:`pandas.NA` as its missing value,
   rather than :attr:`numpy.nan`.

Attributes
----------
None

Methods
-------
None
"""

# create the Dtype
Int8Dtype = register_extension_dtype(
    type(
        "Int8Dtype",
        (_IntegerDtype,),
        {
            "type": np.int8,
            "name": "Int8",
            "__doc__": _dtype_docstring.format(dtype="int8"),
        },
    )
)

Int16Dtype = register_extension_dtype(
    type(
        "Int16Dtype",
        (_IntegerDtype,),
        {
            "type": np.int16,
            "name": "Int16",
            "__doc__": _dtype_docstring.format(dtype="int16"),
        },
    )
)

Int32Dtype = register_extension_dtype(
    type(
        "Int32Dtype",
        (_IntegerDtype,),
        {
            "type": np.int32,
            "name": "Int32",
            "__doc__": _dtype_docstring.format(dtype="int32"),
        },
    )
)

Int64Dtype = register_extension_dtype(
    type(
        "Int64Dtype",
        (_IntegerDtype,),
        {
            "type": np.int64,
            "name": "Int64",
            "__doc__": _dtype_docstring.format(dtype="int64"),
        },
    )
)

UInt8Dtype = register_extension_dtype(
    type(
        "UInt8Dtype",
        (_IntegerDtype,),
        {
            "type": np.uint8,
            "name": "UInt8",
            "__doc__": _dtype_docstring.format(dtype="uint8"),
        },
    )
)

UInt16Dtype = register_extension_dtype(
    type(
        "UInt16Dtype",
        (_IntegerDtype,),
        {
            "type": np.uint16,
            "name": "UInt16",
            "__doc__": _dtype_docstring.format(dtype="uint16"),
        },
    )
)

UInt32Dtype = register_extension_dtype(
    type(
        "UInt32Dtype",
        (_IntegerDtype,),
        {
            "type": np.uint32,
            "name": "UInt32",
            "__doc__": _dtype_docstring.format(dtype="uint32"),
        },
    )
)

UInt64Dtype = register_extension_dtype(
    type(
        "UInt64Dtype",
        (_IntegerDtype,),
        {
            "type": np.uint64,
            "name": "UInt64",
            "__doc__": _dtype_docstring.format(dtype="uint64"),
        },
    )
)

_dtypes = {
    "int8": Int8Dtype(),
    "int16": Int16Dtype(),
    "int32": Int32Dtype(),
    "int64": Int64Dtype(),
    "uint8": UInt8Dtype(),
    "uint16": UInt16Dtype(),
    "uint32": UInt32Dtype(),
    "uint64": UInt64Dtype(),
}
