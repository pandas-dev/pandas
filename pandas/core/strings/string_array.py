from typing import Type, Union

import numpy as np

from pandas._libs import lib, missing as libmissing

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype,
)

from pandas.core.missing import isna
from pandas.core.strings.object_array import ObjectArrayMethods


class StringArrayMethods(ObjectArrayMethods):
    _default_na_value = libmissing.NA

    def _map(self, f, na_value=None, dtype=None):
        from pandas.arrays import BooleanArray, IntegerArray, StringArray
        from pandas.core.arrays.string_ import StringDtype

        if dtype is None:
            dtype = StringDtype()
        if na_value is None:
            na_value = self._default_na_value

        arr = self._array
        mask = isna(arr)

        arr = np.asarray(arr)

        if is_integer_dtype(dtype) or is_bool_dtype(dtype):
            constructor: Union[Type[IntegerArray], Type[BooleanArray]]
            if is_integer_dtype(dtype):
                constructor = IntegerArray
            else:
                constructor = BooleanArray

            na_value_is_na = isna(na_value)
            if na_value_is_na:
                na_value = 1
            result = lib.map_infer_mask(
                arr,
                f,
                mask.view("uint8"),
                convert=False,
                na_value=na_value,
                dtype=np.dtype(dtype),
            )

            if not na_value_is_na:
                mask[:] = False

            return constructor(result, mask)

        elif is_string_dtype(dtype) and not is_object_dtype(dtype):
            # i.e. StringDtype
            result = lib.map_infer_mask(
                arr, f, mask.view("uint8"), convert=False, na_value=na_value
            )
            return StringArray(result)
        else:
            # This is when the result type is object. We reach this when
            # -> We know the result type is truly object (e.g. .encode returns bytes
            #    or .findall returns a list).
            # -> We don't know the result type. E.g. `.get` can return anything.
            return lib.map_infer_mask(arr, f, mask.view("uint8"))
