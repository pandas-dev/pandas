"""
This file is here as an example, this code will live in the Numba and
Bodo libraries.
"""
from __future__ import annotations
from collections.abc import Callable
from typing import TYPE_CHECKING, Literal, Any

import pandas as pd
import bodo

if TYPE_CHECKING:
    from pandas._typing import Axis, AggFuncType

def __pandas_udf__(
    jit_decorator: Callable,
    obj: pd.Series | pd.DataFrame,
    method: Literal["apply", "map"],
    func: AggFuncType,
    axis: Axis,
    raw: bool,
    result_type: Literal["expand", "reduce", "broadcast"] | None,
    args: tuple,
    kwargs: dict[str, Any],
    by_row: Literal[False, "compat"],
):

    if isinstance(obj, pd.DataFrame) and method == "apply":
        if result_type is not None:
            raise NotImplementedError(
                "engine='bodo' not supported when result_type is not None"
            )

        if raw:
            raise NotImplementedError(
                "engine='bodo' not supported when raw=True"
            )
        if isinstance(func, str) and axis != 1:
            raise NotImplementedError(
                "engine='bodo' only supports axis=1 when func is the name of a "
                "user-defined function"
            )
        if args or kwargs:
            raise NotImplementedError(
                "engine='bodo' not supported when args or kwargs are specified"
            )
        @jit_decorator
        def jit_func(df, func, axis):
            return df.apply(func, axis=axis)

        return jit_func(obj, func, axis)
    else:
        raise NotImplementedError(
            f"engine='bodo' not supported for {obj.__class__.__name__}.{method}"
        )

bodo.jit.__pandas_udf__ = __pandas_udf__
