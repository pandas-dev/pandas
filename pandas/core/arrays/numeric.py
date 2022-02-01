from __future__ import annotations

import numbers
from typing import (
    TYPE_CHECKING,
    TypeVar,
)

import numpy as np

from pandas._libs import (
    lib,
    missing as libmissing,
)
from pandas._typing import (
    Dtype,
    DtypeObj,
)
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_string_dtype,
    pandas_dtype,
)

from pandas.core import ops
from pandas.core.arrays.base import ExtensionArray
from pandas.core.arrays.masked import (
    BaseMaskedArray,
    BaseMaskedDtype,
)
from pandas.core.construction import ensure_wrapped_if_datetimelike

if TYPE_CHECKING:
    import pyarrow


T = TypeVar("T", bound="NumericArray")


class NumericDtype(BaseMaskedDtype):
    _default_np_dtype: np.dtype

    def __from_arrow__(
        self, array: pyarrow.Array | pyarrow.ChunkedArray
    ) -> BaseMaskedArray:
        """
        Construct IntegerArray/FloatingArray from pyarrow Array/ChunkedArray.
        """
        import pyarrow

        from pandas.core.arrays._arrow_utils import pyarrow_array_to_numpy_and_mask

        array_class = self.construct_array_type()

        pyarrow_type = pyarrow.from_numpy_dtype(self.type)
        if not array.type.equals(pyarrow_type):
            # test_from_arrow_type_error raise for string, but allow
            #  through itemsize conversion GH#31896
            rt_dtype = pandas_dtype(array.type.to_pandas_dtype())
            if rt_dtype.kind not in ["i", "u", "f"]:
                # Could allow "c" or potentially disallow float<->int conversion,
                #  but at the moment we specifically test that uint<->int works
                raise TypeError(
                    f"Expected array of {self} type, got {array.type} instead"
                )

            array = array.cast(pyarrow_type)

        if isinstance(array, pyarrow.Array):
            chunks = [array]
        else:
            # pyarrow.ChunkedArray
            chunks = array.chunks

        results = []
        for arr in chunks:
            data, mask = pyarrow_array_to_numpy_and_mask(arr, dtype=self.type)
            num_arr = array_class(data.copy(), ~mask, copy=False)
            results.append(num_arr)

        if not results:
            return array_class(
                np.array([], dtype=self.numpy_dtype), np.array([], dtype=np.bool_)
            )
        elif len(results) == 1:
            # avoid additional copy in _concat_same_type
            return results[0]
        else:
            return array_class._concat_same_type(results)

    @classmethod
    def _standardize_dtype(cls, dtype) -> NumericDtype:
        """
        Convert a string representation or a numpy dtype to NumericDtype.
        """
        raise AbstractMethodError(cls)

    @classmethod
    def _safe_cast(cls, values: np.ndarray, dtype: np.dtype, copy: bool) -> np.ndarray:
        """
        Safely cast the values to the given dtype.

        "safe" in this context means the casting is lossless.
        """
        raise AbstractMethodError(cls)


def _coerce_to_data_and_mask(values, mask, dtype, copy, dtype_cls, default_dtype):
    if default_dtype.kind == "f":
        checker = is_float_dtype
    else:
        checker = is_integer_dtype

    inferred_type = None

    if dtype is None and hasattr(values, "dtype"):
        if checker(values.dtype):
            dtype = values.dtype

    if dtype is not None:
        dtype = dtype_cls._standardize_dtype(dtype)

    cls = dtype_cls.construct_array_type()
    if isinstance(values, cls):
        values, mask = values._data, values._mask
        if dtype is not None:
            values = values.astype(dtype.numpy_dtype, copy=False)

        if copy:
            values = values.copy()
            mask = mask.copy()
        return values, mask, dtype, inferred_type

    values = np.array(values, copy=copy)
    inferred_type = None
    if is_object_dtype(values.dtype) or is_string_dtype(values.dtype):
        inferred_type = lib.infer_dtype(values, skipna=True)
        if inferred_type == "empty":
            pass
        elif inferred_type == "boolean":
            name = dtype_cls.__name__.strip("_")
            raise TypeError(f"{values.dtype} cannot be converted to {name}")

    elif is_bool_dtype(values) and checker(dtype):
        values = np.array(values, dtype=default_dtype, copy=copy)

    elif not (is_integer_dtype(values) or is_float_dtype(values)):
        name = dtype_cls.__name__.strip("_")
        raise TypeError(f"{values.dtype} cannot be converted to {name}")

    if values.ndim != 1:
        raise TypeError("values must be a 1D list-like")

    if mask is None:
        mask = libmissing.is_numeric_na(values)
    else:
        assert len(mask) == len(values)

    if mask.ndim != 1:
        raise TypeError("mask must be a 1D list-like")

    # infer dtype if needed
    if dtype is None:
        dtype = default_dtype
    else:
        dtype = dtype.type

    # we copy as need to coerce here
    if mask.any():
        values = values.copy()
        values[mask] = cls._internal_fill_value
    if inferred_type in ("string", "unicode"):
        # casts from str are always safe since they raise
        # a ValueError if the str cannot be parsed into a float
        values = values.astype(dtype, copy=copy)
    else:
        values = dtype_cls._safe_cast(values, dtype, copy=False)

    return values, mask, dtype, inferred_type


class NumericArray(BaseMaskedArray):
    """
    Base class for IntegerArray and FloatingArray.
    """

    _dtype_cls: type[NumericDtype]

    @classmethod
    def _coerce_to_array(
        cls, value, *, dtype: DtypeObj, copy: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        dtype_cls = cls._dtype_cls
        default_dtype = dtype_cls._default_np_dtype
        mask = None
        values, mask, _, _ = _coerce_to_data_and_mask(
            value, mask, dtype, copy, dtype_cls, default_dtype
        )
        return values, mask

    @classmethod
    def _from_sequence_of_strings(
        cls: type[T], strings, *, dtype: Dtype | None = None, copy: bool = False
    ) -> T:
        from pandas.core.tools.numeric import to_numeric

        scalars = to_numeric(strings, errors="raise")
        return cls._from_sequence(scalars, dtype=dtype, copy=copy)

    def _arith_method(self, other, op):
        op_name = op.__name__
        omask = None

        if isinstance(other, BaseMaskedArray):
            other, omask = other._data, other._mask

        elif is_list_like(other):
            if not isinstance(other, ExtensionArray):
                other = np.asarray(other)
            if other.ndim > 1:
                raise NotImplementedError("can only perform ops with 1-d structures")

        # We wrap the non-masked arithmetic logic used for numpy dtypes
        #  in Series/Index arithmetic ops.
        other = ops.maybe_prepare_scalar_for_op(other, (len(self),))
        pd_op = ops.get_array_op(op)
        other = ensure_wrapped_if_datetimelike(other)

        mask = self._propagate_mask(omask, other)

        if other is libmissing.NA:
            result = np.ones_like(self._data)
            if "truediv" in op_name and self.dtype.kind != "f":
                # The actual data here doesn't matter since the mask
                #  will be all-True, but since this is division, we want
                #  to end up with floating dtype.
                result = result.astype(np.float64)
        else:
            # Make sure we do this before the "pow" mask checks
            #  to get an expected exception message on shape mismatch.
            if self.dtype.kind in ["i", "u"] and op_name in ["floordiv", "mod"]:
                # ATM we don't match the behavior of non-masked types with
                #  respect to floordiv-by-zero
                pd_op = op

            with np.errstate(all="ignore"):
                result = pd_op(self._data, other)

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

        return self._maybe_mask_result(result, mask, other, op_name)

    _HANDLED_TYPES = (np.ndarray, numbers.Number)

    def __neg__(self):
        return type(self)(-self._data, self._mask.copy())

    def __pos__(self):
        return self.copy()

    def __abs__(self):
        return type(self)(abs(self._data), self._mask.copy())

    def round(self: T, decimals: int = 0, *args, **kwargs) -> T:
        """
        Round each value in the array a to the given number of decimals.

        Parameters
        ----------
        decimals : int, default 0
            Number of decimal places to round to. If decimals is negative,
            it specifies the number of positions to the left of the decimal point.
        *args, **kwargs
            Additional arguments and keywords have no effect but might be
            accepted for compatibility with NumPy.

        Returns
        -------
        NumericArray
            Rounded values of the NumericArray.

        See Also
        --------
        numpy.around : Round values of an np.array.
        DataFrame.round : Round values of a DataFrame.
        Series.round : Round values of a Series.
        """
        nv.validate_round(args, kwargs)
        values = np.round(self._data, decimals=decimals, **kwargs)
        return type(self)(values, self._mask.copy())
