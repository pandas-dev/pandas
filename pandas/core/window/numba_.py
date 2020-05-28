from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

from pandas._typing import Scalar
from pandas.compat._optional import import_optional_dependency

from pandas.core.util.numba_ import (
    check_kwargs_and_nopython,
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
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs)

    check_kwargs_and_nopython(kwargs, nopython)

    numba_func = jit_user_function(func, nopython, nogil, parallel)

    numba = import_optional_dependency("numba")

    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def roll_apply(
        values: np.ndarray, begin: np.ndarray, end: np.ndarray, minimum_periods: int,
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
