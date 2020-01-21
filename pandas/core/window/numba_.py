import types
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

from pandas._typing import Scalar
from pandas.compat._optional import import_optional_dependency


def make_rolling_apply(
    func: Callable[..., Scalar],
    args: Tuple,
    nogil: bool,
    parallel: bool,
    nopython: bool,
):
    """
    Creates a JITted rolling apply function with a JITted version of
    the user's function.

    Parameters
    ----------
    func : function
        function to be applied to each window and will be JITed
    args : tuple
        *args to be passed into the function
    nogil : bool
        nogil parameter from engine_kwargs for numba.jit
    parallel : bool
        parallel parameter from engine_kwargs for numba.jit
    nopython : bool
        nopython parameter from engine_kwargs for numba.jit

    Returns
    -------
    Numba function
    """
    numba = import_optional_dependency("numba")

    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    if isinstance(func, numba.targets.registry.CPUDispatcher):
        # Don't jit a user passed jitted function
        numba_func = func
    else:

        @numba.generated_jit(nopython=nopython, nogil=nogil, parallel=parallel)
        def numba_func(window, *_args):
            if getattr(np, func.__name__, False) is func or isinstance(
                func, types.BuiltinFunctionType
            ):
                jf = func
            else:
                jf = numba.jit(func, nopython=nopython, nogil=nogil)

            def impl(window, *_args):
                return jf(window, *_args)

            return impl

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

    if engine_kwargs is None:
        engine_kwargs = {}

    nopython = engine_kwargs.get("nopython", True)
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)

    if kwargs and nopython:
        raise ValueError(
            "numba does not support kwargs with nopython=True: "
            "https://github.com/numba/numba/issues/2916"
        )

    return make_rolling_apply(func, args, nogil, parallel, nopython)
