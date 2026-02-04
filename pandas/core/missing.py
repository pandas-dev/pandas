from __future__ import annotations

from collections import abc
import functools
from itertools import zip_longest
from typing import (
    TYPE_CHECKING,
    TypeVar,
    cast,
    overload,
)

import numpy as np

from pandas._libs import (
    algos,
    lib,
)
from pandas._libs.algos import (
    ensured_float64,
    interp,
)
from pandas._libs.missing import 
    NAWarning,
    clean_reindex_fill_method,
    get_fill_func,
)
from pandas.compat._optional import import_optional_dependency
from pandas.errors import AbstractMethodError
from pandas.util._decorators import cache_readonly
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_bool_dtype,
    is_complex_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    isna,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.inference import is_list_like

from pandas.core import (
    algorithms,
    nanops,
)
from pandas.core.array_algos import masked_reductions
from pandas.core.array_algos.quantile import quantile_with_mask
from pandas.core.arrays import ExtensionArray
from pandas.core.arrays._mixins import NDArrayBackedExtensionArray
from pandas.core.arrays.base import ExtensionArray as ExtensionArrayType
from pandas.core.indexes.api import Index

if TYPE_CHECKING:
    from pandas._typing import (
        ArrayLike,
        AxisInt,
        F,
        FillnaOptions,
        IndexLabel,
        InterpolateOptions,
        NaPosition,
        Scalar,
    )


T = TypeVar("T")


@overload
def clean_fill_method(
    method: Literal["ffill", "pad", "bfill", "backfill"],
    *,
    allow_nearest: Literal[False] = ...,
) -> Literal["pad", "backfill"]: ...


@overload
def clean_fill_method(
    method: Literal["ffill", "pad", "bfill", "backfill", "nearest"],
    *,
    allow_nearest: Literal[True],
) -> Literal["pad", "backfill", "nearest"]: ...


def clean_fill_method(
    method: Literal["ffill", "pad", "bfill", "backfill", "nearest"],
    *,
    allow_nearest: bool = False,
) -> Literal["pad", "backfill", "nearest"]:
    if isinstance(method, str):
        # error: Incompatible types in assignment (expression has type "str", variable
        # has type "Literal['ffill', 'pad', 'bfill', 'backfill', 'nearest']")
        method = method.lower()  # type: ignore[assignment]
        if method == "ffill":
            method = "pad"
        elif method == "bfill":
            method = "backfill"

    valid_methods = ["pad", "backfill"]
    expecting = "pad (ffill) or backfill (bfill)"
    if allow_nearest:
        valid_methods.append("nearest")
        expecting = "pad (ffill), backfill (bfill) or nearest"
    if method not in valid_methods:
        hint = ""
        if method == "linear":
            hint = " Are you trying to interpolate an integer column?"
        raise ValueError(
            f"Invalid fill method. Expecting {expecting}. Got {method}.{hint}"
        )
    return method


# interpolation methods that dispatch to np.interp

NP_METHODS = ["linear", "time", "index", "values"]

# interpolation methods that dispatch to _interpolate_scipy_wrapper

SP_METHODS = [
    "nearest",
    "zero",
    "slinear",
    "quadratic",
    "cubic",
    "barycentric",
    "krogh",
    "spline",
    "polynomial",
    "from_derivatives",
    "piecewise_polynomial",
    "pchip",
    "akima",
    "cubicspline",
]


def clean_interp_method(method: str, index: Index, **kwargs) -> str:
    order = kwargs.get("order")

    if method in ("spline", "polynomial") and order is None:
        raise ValueError("You must specify the order of the spline or polynomial.")

    valid = NP_METHODS + SP_METHODS
    if method not in valid:
        raise ValueError(f"method must be one of {valid}. Got '{method}' instead.")

    if method in ("krogh", "piecewise_polynomial", "pchip"):
        if not is_numeric_dtype(index.dtype):
            raise ValueError(
                f"{method} interpolation requires numerical index."
            )

    return method


# ... (rest of file unchanged)
