"""Common utilities for Numba operations"""

from __future__ import annotations

import inspect
import types
from typing import TYPE_CHECKING

import numpy as np

from pandas.compat._optional import import_optional_dependency
from pandas.errors import NumbaUtilError

if TYPE_CHECKING:
    from collections.abc import Callable

GLOBAL_USE_NUMBA: bool = False


def maybe_use_numba(engine: str | None) -> bool:
    """Signal whether to use numba routines."""
    return engine == "numba" or (engine is None and GLOBAL_USE_NUMBA)


def set_use_numba(enable: bool = False) -> None:
    global GLOBAL_USE_NUMBA
    if enable:
        import_optional_dependency("numba")
    GLOBAL_USE_NUMBA = enable


def get_jit_arguments(engine_kwargs: dict[str, bool] | None = None) -> dict[str, bool]:
    """
    Return arguments to pass to numba.JIT, falling back on pandas default JIT settings.

    Parameters
    ----------
    engine_kwargs : dict, default None
        user passed keyword arguments for numba.JIT

    Returns
    -------
    dict[str, bool]
        nopython, nogil, parallel

    Raises
    ------
    NumbaUtilError
    """
    if engine_kwargs is None:
        engine_kwargs = {}

    nopython = engine_kwargs.get("nopython", True)
    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)
    return {"nopython": nopython, "nogil": nogil, "parallel": parallel}


def jit_user_function(func: Callable) -> Callable:
    """
    If user function is not jitted already, mark the user's function
    as jitable.

    Parameters
    ----------
    func : function
        user defined function

    Returns
    -------
    function
        Numba JITed function, or function marked as JITable by numba
    """
    if TYPE_CHECKING:
        import numba
    else:
        numba = import_optional_dependency("numba")

    if numba.extending.is_jitted(func):
        # Don't jit a user passed jitted function
        numba_func = func
    elif getattr(np, func.__name__, False) is func or isinstance(
        func, types.BuiltinFunctionType
    ):
        # Not necessary to jit builtins or np functions
        # This will mess up register_jitable
        numba_func = func
    else:
        numba_func = numba.extending.register_jitable(func)

    return numba_func


_sentinel = object()


def prepare_function_arguments(
    func: Callable, args: tuple, kwargs: dict, *, num_required_args: int
) -> tuple[tuple, dict]:
    """
    Prepare arguments for jitted function. As numba functions do not support kwargs,
    we try to move kwargs into args if possible.

    Parameters
    ----------
    func : function
        User defined function
    args : tuple
        User input positional arguments
    kwargs : dict
        User input keyword arguments
    num_required_args : int
        The number of leading positional arguments we will pass to udf.
        These are not supplied by the user.
        e.g. for groupby we require "values", "index" as the first two arguments:
        `numba_func(group, group_index, *args)`, in this case num_required_args=2.
        See :func:`pandas.core.groupby.numba_.generate_numba_agg_func`

    Returns
    -------
    tuple[tuple, dict]
        args, kwargs

    """
    if not kwargs:
        return args, kwargs

    # the udf should have this pattern: def udf(arg1, arg2, ..., *args, **kwargs):...
    signature = inspect.signature(func)
    arguments = signature.bind(*[_sentinel] * num_required_args, *args, **kwargs)
    arguments.apply_defaults()
    # Ref: https://peps.python.org/pep-0362/
    # Arguments which could be passed as part of either *args or **kwargs
    # will be included only in the BoundArguments.args attribute.
    args = arguments.args
    kwargs = arguments.kwargs

    if kwargs:
        # Note: in case numba supports keyword-only arguments in
        # a future version, we should remove this check. But this
        # seems unlikely to happen soon.

        raise NumbaUtilError(
            "numba does not support keyword-only arguments"
            "https://github.com/numba/numba/issues/2916, "
            "https://github.com/numba/numba/issues/6846"
        )

    args = args[num_required_args:]
    return args, kwargs
