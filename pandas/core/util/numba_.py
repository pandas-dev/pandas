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
        nogil, parallel

    Raises
    ------
    NumbaUtilError
    """
    if engine_kwargs is None:
        engine_kwargs = {}

    nogil = engine_kwargs.get("nogil", False)
    parallel = engine_kwargs.get("parallel", False)
    return {"nogil": nogil, "parallel": parallel}


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


_JIT_APPLY_CACHE: dict[int, Callable] = {}
_GLOBAL_LOOP: Callable | None = None


def get_global_loop(numba) -> Callable:
    global _GLOBAL_LOOP
    if _GLOBAL_LOOP is None:
        @numba.njit(parallel=True)
        def global_numba_apply_loop(arr, jitted_func):
            n = len(arr)
            first_val = jitted_func(arr[0])
            res = np.empty(n, dtype=type(first_val))
            res[0] = first_val
            for i in numba.prange(1, n):
                res[i] = jitted_func(arr[i])
            return res
        _GLOBAL_LOOP = global_numba_apply_loop
    return _GLOBAL_LOOP


def maybe_run_numba_apply(
    series,
    func,
    engine: str | None = None,
    engine_kwargs: dict[str, Any] | None = None,
) -> np.ndarray | None:
    """
    Attempt to JIT compile the user function using Numba and run it over the Series.

    Parameters
    ----------
    series : Series
        The Pandas Series object.
    func : Callable
        The User Defined Function (UDF).
    engine : str, optional
        The execution engine (e.g. 'numba' or 'python').
    engine_kwargs : dict, optional
        Keyword arguments passed to numba.njit.

    Returns
    -------
    np.ndarray or None
        The result array if compilation and execution succeed, otherwise None.
    """
    if not maybe_use_numba(engine):
        return None

    # Size threshold only applies to automatic JIT (engine is None and GLOBAL_USE_NUMBA is True)
    if engine is None and len(series) < 50_000:
        return None

    if not np.issubdtype(series.dtype, np.number) and not np.issubdtype(series.dtype, np.datetime64) and not np.issubdtype(series.dtype, np.timedelta64):
        if engine == "numba":
            raise ValueError(f"Numba engine only supports numeric/datetime dtypes, got {series.dtype}")
        return None

    try:
        import numba
    except ImportError as err:
        if engine == "numba":
            raise ImportError("Numba is not installed. Please install numba to use engine='numba'") from err
        return None

    func_id = id(func)
    if func_id in _JIT_APPLY_CACHE:
        jitted_udf = _JIT_APPLY_CACHE[func_id]
    else:
        try:
            nopython = True
            nogil = True
            parallel = False
            if engine_kwargs:
                nopython = engine_kwargs.get("nopython", True)
                nogil = engine_kwargs.get("nogil", True)
                parallel = engine_kwargs.get("parallel", False)
            jitted_udf = numba.njit(func, nopython=nopython, nogil=nogil, parallel=parallel)
            _JIT_APPLY_CACHE[func_id] = jitted_udf
        except Exception as err:
            if engine == "numba":
                raise NumbaUtilError(f"Failed to JIT compile function: {err}") from err
            return None

    try:
        loop = get_global_loop(numba)
        return loop(series.to_numpy(), jitted_udf)
    except Exception as err:
        if engine == "numba":
            raise NumbaUtilError(f"Failed to execute Numba JIT loop: {err}") from err
        return None




