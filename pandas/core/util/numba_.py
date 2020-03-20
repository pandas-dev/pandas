"""Common utilities for Numba operations"""
import inspect
import types
from typing import Callable, Dict, Optional

import numpy as np

from pandas.compat._optional import import_optional_dependency
from pandas._typing import FrameOrSeries


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


def split_for_numba(arg: FrameOrSeries):
    """
    Split pandas object into its components as numpy arrays for numba functions.
    """
    if getattr(arg, "columns", None) is not None:
        columns_as_array = arg.columns.to_numpy()
    else:
        columns_as_array = None
    return arg.to_numpy(), arg.index.to_numpy(), columns_as_array


def validate_udf(func: Callable, include_columns: bool = False):
    """
    Validate user defined function for ops when using Numba.

    For routines that pass Series objects, the first signature arguments should include:

    def f(values, index, ...):
        ...

    For routines that pass DataFrame objects, the first signature arguments should
    include:

    def f(values, index, columns, ...):
        ...
    """
    udf_signature = list(inspect.signature(func).parameters.keys())
    expected_args = ["values", "index"]
    if include_columns:
        expected_args.append("columns")
    min_number_args = len(expected_args)
    if (
        len(udf_signature) < min_number_args
        or udf_signature[:min_number_args] != expected_args
    ):
        raise ValueError(
            f"The first {min_number_args} arguments to {func.__name__} must be "
            f"{expected_args}"
        )
