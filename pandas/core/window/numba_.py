from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

from pandas._typing import Scalar
from pandas.compat._optional import import_optional_dependency

from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    get_jit_arguments,
    jit_user_function,
)


def generate_numba_apply_func(
    args: Tuple,
    kwargs: Dict[str, Any],
    func: Callable[..., Scalar],
    engine_kwargs: Optional[Dict[str, bool]],
):
    """
    Generate a numba jitted apply function specified by values from engine_kwargs.

    1. jit the user's function
    2. Return a rolling apply function with the jitted function inline

    Configurations specified in engine_kwargs apply to both the user's
    function _AND_ the rolling apply function.

    Parameters
    ----------
    args : tuple
        *args to be passed into the function
    kwargs : dict
        **kwargs to be passed into the function
    func : function
        function to be applied to each window and will be JITed
    engine_kwargs : dict
        dictionary of arguments to be passed into numba.jit

    Returns
    -------
    Numba function
    """
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs, kwargs)

    cache_key = (func, "rolling_apply")
    if cache_key in NUMBA_FUNC_CACHE:
        return NUMBA_FUNC_CACHE[cache_key]

    numba_func = jit_user_function(func, nopython, nogil, parallel)
    numba = import_optional_dependency("numba")
    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def roll_apply(
        values: np.ndarray, begin: np.ndarray, end: np.ndarray, minimum_periods: int
    ) -> np.ndarray:
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


def generate_numba_groupby_ewma_func(
    engine_kwargs: Optional[Dict[str, bool]],
    com: float,
    adjust: bool,
    ignore_na: bool,
):
    """
    Generate a numba jitted groupby ewma function specified by values
    from engine_kwargs.

    Parameters
    ----------
    engine_kwargs : dict
        dictionary of arguments to be passed into numba.jit
    com : float
    adjust : bool
    ignore_na : bool

    Returns
    -------
    Numba function
    """
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs)

    cache_key = (lambda x: x, "groupby_ewma")
    if cache_key in NUMBA_FUNC_CACHE:
        return NUMBA_FUNC_CACHE[cache_key]

    numba = import_optional_dependency("numba")
    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def groupby_ewma(
        values: np.ndarray,
        begin: np.ndarray,
        end: np.ndarray,
        minimum_periods: int,
    ) -> np.ndarray:
        result = np.empty(len(values))
        alpha = 1.0 / (1.0 + com)
        for i in loop_range(len(begin)):
            start = begin[i]
            stop = end[i]
            window = values[start:stop]
            sub_result = np.empty(len(window))

            old_wt_factor = 1.0 - alpha
            new_wt = 1.0 if adjust else alpha

            weighted_avg = window[0]
            nobs = int(not np.isnan(weighted_avg))
            sub_result[0] = weighted_avg if nobs >= minimum_periods else np.nan
            old_wt = 1.0

            for j in range(1, len(window)):
                cur = window[j]
                is_observation = not np.isnan(cur)
                nobs += is_observation
                if not np.isnan(weighted_avg):

                    if is_observation or not ignore_na:

                        old_wt *= old_wt_factor
                        if is_observation:

                            # avoid numerical errors on constant series
                            if weighted_avg != cur:
                                weighted_avg = (
                                    (old_wt * weighted_avg) + (new_wt * cur)
                                ) / (old_wt + new_wt)
                            if adjust:
                                old_wt += new_wt
                            else:
                                old_wt = 1.0
                elif is_observation:
                    weighted_avg = cur

                sub_result[j] = weighted_avg if nobs >= minimum_periods else np.nan

            result[start:stop] = sub_result

        return result

    return groupby_ewma
