"""
This file is here as an example, this code will live in the Numba and
Bodo libraries.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

import bodo
import numpy as np

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import (
        AggFuncType,
        Axis,
    )


class BodoExecutionEngine(pd.api.executors.BaseExecutionEngine):
    @staticmethod
    def map(
        data: pd.Series | pd.DataFrame | np.ndarray,
        func: AggFuncType,
        args: tuple,
        kwargs: dict[str, Any],
        decorator: Callable,
        skip_na: bool,
    ):
        raise NotImplementedError("engine='bodo' not supported for map")

    @staticmethod
    def apply(
        data: pd.Series | pd.DataFrame | np.ndarray,
        func: AggFuncType,
        args: tuple,
        kwargs: dict[str, Any],
        decorator: Callable,
        axis: Axis,
    ):
        if isinstance(data, pd.Series):
            raise NotImplementedError("engine='bodo' not supported for Series.apply")

        if isinstance(data, np.ndarray):
            raise NotImplementedError("engine='bodo' not supported when raw=True")

        if args or kwargs:
            raise NotImplementedError(
                "engine='bodo' not supported when args or kwargs are specified"
            )

        if isinstance(func, str) and axis != 1:
            raise NotImplementedError(
                "engine='bodo' only supports axis=1 when func is the name of a "
                "user-defined function"
            )

        def jit_func(df, func, axis):
            return df.apply(func, axis=axis)

        jit_func = decorator(jit_func)

        return jit_func(data, func, axis)


bodo.jit.__pandas_udf__ = BodoExecutionEngine
