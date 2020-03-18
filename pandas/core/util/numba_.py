"""Common utilities for Numba operations"""
import types
from typing import Callable, Dict, Optional

import numpy as np

from pandas.compat._optional import import_optional_dependency


def check_kwargs_and_nopython(
    kwargs: Optional[Dict] = None, nopython: Optional[bool] = None
):
    if kwargs and nopython:
        raise ValueError(
            "numba does not support kwargs with nopython=True: "
            "https://github.com/numba/numba/issues/2916"
        )


def get_jit_arguments(engine_kwargs: Optional[Dict[str, bool]] = None):
    """
    Return arguments to pass to numba.JIT, falling back on pandas default JIT settings.
    """
    if engine_kwargs is None:
        engine_kwargs = {}

    nopython = engine_kwargs.get("nopython", True)
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)
    return nopython, nogil, parallel


def jit_user_function(func: Callable, nopython: bool, nogil: bool, parallel: bool):
    """
    JIT the user's function given the configurable arguments.
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
