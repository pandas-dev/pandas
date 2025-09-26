from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np

from pandas._config import is_nan_na

from pandas._libs import lib
from pandas._libs.missing import NA
from pandas.errors import LossySetitemError

from pandas.core.dtypes.cast import np_can_hold_element
from pandas.core.dtypes.common import is_numeric_dtype

if TYPE_CHECKING:
    from pandas._typing import (
        npt,
    )

    from pandas.core.arrays.base import ExtensionArray


def to_numpy_dtype_inference(
    arr: ExtensionArray,
    dtype: npt.DTypeLike | None,
    na_value,
    hasna: bool,
) -> tuple[np.dtype | None, Any]:
    result_dtype: np.dtype | None
    inferred_numeric_dtype = False
    if dtype is None and is_numeric_dtype(arr.dtype):
        inferred_numeric_dtype = True
        if hasna:
            if arr.dtype.kind == "b":
                result_dtype = np.dtype(np.object_)
            else:
                if arr.dtype.kind in "iu":
                    result_dtype = np.dtype(np.float64)
                else:
                    result_dtype = arr.dtype.numpy_dtype  # type: ignore[attr-defined]
                if na_value is lib.no_default:
                    if not is_nan_na():
                        na_value = NA
                        dtype = np.dtype(object)
                    else:
                        na_value = np.nan
        else:
            result_dtype = arr.dtype.numpy_dtype  # type: ignore[attr-defined]
    elif dtype is not None:
        result_dtype = np.dtype(dtype)
    else:
        result_dtype = None

    if na_value is lib.no_default:
        if result_dtype is None or not hasna:
            na_value = arr.dtype.na_value
        elif result_dtype.kind == "f":
            na_value = np.nan
        elif result_dtype.kind == "M":
            na_value = np.datetime64("nat")
        elif result_dtype.kind == "m":
            na_value = np.timedelta64("nat")
        else:
            na_value = arr.dtype.na_value

    if inferred_numeric_dtype and hasna:
        try:
            np_can_hold_element(result_dtype, na_value)  # type: ignore[arg-type]
        except LossySetitemError:
            result_dtype = np.dtype(np.object_)
    return result_dtype, na_value
