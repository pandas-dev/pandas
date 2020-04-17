"""Common utilities for Numba operations"""
import inspect
import types
from typing import Callable, Dict, Optional, Tuple

import numpy as np

from pandas._typing import FrameOrSeries
from pandas.compat._optional import import_optional_dependency


def check_kwargs_and_nopython(
    kwargs: Optional[Dict] = None, nopython: Optional[bool] = None
) -> None:
    """
    Validate that **kwargs and nopython=True was passed
    https://github.com/numba/numba/issues/2916

    Parameters
    ----------
    kwargs : dict, default None
        user passed keyword arguments to pass into the JITed function
    nopython : bool, default None
        nopython parameter

    Returns
    -------
    None

    Raises
    ------
    ValueError
    """
    if kwargs and nopython:
        raise ValueError(
            "numba does not support kwargs with nopython=True: "
            "https://github.com/numba/numba/issues/2916"
        )


def get_jit_arguments(
    engine_kwargs: Optional[Dict[str, bool]] = None
) -> Tuple[bool, bool, bool]:
    """
    Return arguments to pass to numba.JIT, falling back on pandas default JIT settings.

    Parameters
    ----------
    engine_kwargs : dict, default None
        user passed keyword arguments for numba.JIT

    Returns
    -------
    (bool, bool, bool)
        nopython, nogil, parallel
    """
    if engine_kwargs is None:
        engine_kwargs = {}

    nopython = engine_kwargs.get("nopython", True)
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)
    return nopython, nogil, parallel


def jit_user_function(
    func: Callable, nopython: bool, nogil: bool, parallel: bool
) -> Callable:
    """
    JIT the user's function given the configurable arguments.

    Parameters
    ----------
    func : function
        user defined function
    nopython : bool
        nopython parameter for numba.JIT
    nogil : bool
        nogil parameter for numba.JIT
    parallel : bool
        parallel parameter for numba.JIT

    Returns
    -------
    function
        Numba JITed function
    """
    numba = import_optional_dependency("numba")

    if isinstance(func, numba.targets.registry.CPUDispatcher):
        # Don't jit a user passed jitted function
        numba_func = func
    else:

        @numba.generated_jit(nopython=nopython, nogil=nogil, parallel=parallel)
        def numba_func(data, *_args):
            if getattr(np, func.__name__, False) is func or isinstance(
                func, types.BuiltinFunctionType
            ):
                jf = func
            else:
                jf = numba.jit(func, nopython=nopython, nogil=nogil)

            def impl(data, *_args):
                return jf(data, *_args)

            return impl

    return numba_func


def split_for_numba(arg: FrameOrSeries) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split pandas object into its components as numpy arrays for numba functions.

    Parameters
    ----------
    arg : Series or DataFrame

    Returns
    -------
    (ndarray, ndarray)
        values, index
    """
    return arg.to_numpy(), arg.index.to_numpy()


def validate_udf(func: Callable) -> None:
    """
    Validate user defined function for ops when using Numba.

    The first signature arguments should include:

    def f(values, index, ...):
        ...

    Parameters
    ----------
    func : function, default False
        user defined function

    Returns
    -------
    None
    """
    udf_signature = list(inspect.signature(func).parameters.keys())
    expected_args = ["values", "index"]
    min_number_args = len(expected_args)
    if (
        len(udf_signature) < min_number_args
        or udf_signature[:min_number_args] != expected_args
    ):
        raise ValueError(
            f"The first {min_number_args} arguments to {func.__name__} must be "
            f"{expected_args}"
        )
