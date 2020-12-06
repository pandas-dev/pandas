import datetime

import numpy as np

from pandas._libs import Timedelta, missing as libmissing
from pandas.errors import AbstractMethodError

from pandas.core.dtypes.common import (
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
)

from .masked import BaseMaskedArray


class NumericArray(BaseMaskedArray):
    """
    Base class for IntegerArray and FloatingArray.
    """

    def _maybe_mask_result(self, result, mask, other, op_name: str):
        raise AbstractMethodError(self)

    def _arith_method(self, other, op):
        op_name = op.__name__
        omask = None

        if getattr(other, "ndim", 0) > 1:
            raise NotImplementedError("can only perform ops with 1-d structures")

        if isinstance(other, NumericArray):
            other, omask = other._data, other._mask

        elif is_list_like(other):
            other = np.asarray(other)
            if other.ndim > 1:
                raise NotImplementedError("can only perform ops with 1-d structures")
            if len(self) != len(other):
                raise ValueError("Lengths must match")
            if not (is_float_dtype(other) or is_integer_dtype(other)):
                raise TypeError("can only perform ops with numeric values")

        elif isinstance(other, (datetime.timedelta, np.timedelta64)):
            other = Timedelta(other)

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
