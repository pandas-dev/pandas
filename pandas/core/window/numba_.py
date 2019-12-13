from typing import Callable, Dict, Optional, Tuple

import numpy as np

from pandas.compat._optional import import_optional_dependency


def _generate_numba_apply_func(
    args: Tuple, kwargs: Dict, func: Callable, engine_kwargs: Optional[Dict]
):
    """
    Generate a numba jitted apply function specified by values from engine_kwargs.

    1. jit the user's function
    2. Return a rolling apply function with the jitted function inline

    Configurations specified in engine_kwargs apply to both the user's
    function _AND_ the rolling apply function.
    """
    numba = import_optional_dependency("numba")

    if engine_kwargs is None:
        engine_kwargs = {"nopython": True, "nogil": False, "parallel": False}

    nopython = engine_kwargs.get("nopython", True)
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)

    if kwargs and nopython:
        raise ValueError(
            "numba does not support kwargs with nopython=True: "
            "https://github.com/numba/numba/issues/2916"
        )

    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    def make_rolling_apply(func):
        @numba.generated_jit(nopython=nopython)
        def numba_func(window, *_args):
            if getattr(np, func.__name__, False) is func:

                def impl(window, *_args):
                    return func(window, *_args)

                return impl
            else:
                jf = numba.jit(func, nopython=nopython)

                def impl(window, *_args):
                    return jf(window, *_args)

                return impl

        @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
        def roll_apply(
            values: np.ndarray,
            begin: np.ndarray,
            end: np.ndarray,
            minimum_periods: int,
        ):
            result = np.empty(len(begin))
            for i in loop_range(len(result)):
                start = begin[i]
                stop = end[i]
                window = values[start:stop]
                count_nan = np.sum(np.isnan(window))
                if len(window) - count_nan >= minimum_periods:
                    result[i] = numba_func(window, *args)
                else:
                    result[i] = np.nan
            return result

        return roll_apply

    return make_rolling_apply(func)
