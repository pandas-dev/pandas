from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)
import warnings

import numpy as np

from pandas._config import is_nan_na

from pandas._libs import lib
from pandas._libs.missing import NA
from pandas._libs.tslibs import is_unitless
from pandas.errors import (
    LossySetitemError,
    Pandas4Warning,
)
from pandas.util._exceptions import find_stack_level

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
            unit = np.datetime_data(result_dtype)[0]  # type: ignore[arg-type]
            na_value = np.datetime64("NaT", unit)  # type: ignore[call-overload]
        elif result_dtype.kind == "m":
            unit = np.datetime_data(result_dtype)[0]  # type: ignore[arg-type]
            na_value = np.timedelta64("NaT", unit)  # type: ignore[call-overload]
        else:
            na_value = arr.dtype.na_value

    if inferred_numeric_dtype and hasna:
        try:
            np_can_hold_element(result_dtype, na_value)  # type: ignore[arg-type]
        except LossySetitemError:
            result_dtype = np.dtype(np.object_)
    return result_dtype, na_value


def warn_unitless_to_numpy_dtype(dtype: npt.DTypeLike | None) -> bool:
    """
    Warn if ``dtype`` is a unit-less datetime64 or timedelta64 dtype.

    Returns whether it warned.
    """
    if dtype is None:
        return False
    try:
        dtype = np.dtype(dtype)
    except TypeError:
        # not a numpy dtype, so not unit-less; any error is raised by the caller
        return False
    if dtype.kind in "mM" and is_unitless(dtype):
        # GH#59772
        warnings.warn(
            f"Passing a unit-less dtype '{dtype.name}' to to_numpy is deprecated "
            "and will raise in a future version. "
            f"Pass e.g. '{dtype.name}[ns]' instead.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
        return True
    return False
