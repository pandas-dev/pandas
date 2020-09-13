"""Common utilities for Numba operations with groupby ops"""
import inspect
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

from pandas._typing import FrameOrSeries, Scalar
from pandas.compat._optional import import_optional_dependency

from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    NumbaUtilError,
    check_kwargs_and_nopython,
    get_jit_arguments,
    jit_user_function,
)


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
    Validate user defined function for ops when using Numba with groupby ops.

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

    Raises
    ------
    NumbaUtilError
    """
    udf_signature = list(inspect.signature(func).parameters.keys())
    expected_args = ["values", "index"]
    min_number_args = len(expected_args)
    if (
        len(udf_signature) < min_number_args
        or udf_signature[:min_number_args] != expected_args
    ):
        raise NumbaUtilError(
            f"The first {min_number_args} arguments to {func.__name__} must be "
            f"{expected_args}"
        )


def generate_numba_func(
    func: Callable,
    engine_kwargs: Optional[Dict[str, bool]],
    kwargs: dict,
    cache_key_str: str,
) -> Tuple[Callable, Tuple[Callable, str]]:
    """
    Return a JITed function and cache key for the NUMBA_FUNC_CACHE

    This _may_ be specific to groupby (as it's only used there currently).

    Parameters
    ----------
    func : function
        user defined function
    engine_kwargs : dict or None
        numba.jit arguments
    kwargs : dict
        kwargs for func
    cache_key_str : str
        string representing the second part of the cache key tuple

    Returns
    -------
    (JITed function, cache key)

    Raises
    ------
    NumbaUtilError
    """
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs)
    check_kwargs_and_nopython(kwargs, nopython)
    validate_udf(func)
    cache_key = (func, cache_key_str)
    numba_func = NUMBA_FUNC_CACHE.get(
        cache_key, jit_user_function(func, nopython, nogil, parallel)
    )
    return numba_func, cache_key


def generate_numba_agg_func(
    args: Tuple,
    kwargs: Dict[str, Any],
    func: Callable[..., Scalar],
    engine_kwargs: Optional[Dict[str, bool]],
) -> Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int], np.ndarray]:
    """
    Generate a numba jitted agg function specified by values from engine_kwargs.

    1. jit the user's function
    2. Return a groupby agg function with the jitted function inline

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

    validate_udf(func)

    numba_func = jit_user_function(func, nopython, nogil, parallel)

    numba = import_optional_dependency("numba")

    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def group_agg(
        values: np.ndarray,
        index: np.ndarray,
        begin: np.ndarray,
        end: np.ndarray,
        num_groups: int,
        num_columns: int,
    ) -> np.ndarray:
        result = np.empty((num_groups, num_columns))
        for i in loop_range(num_groups):
            group_index = index[begin[i] : end[i]]
            for j in loop_range(num_columns):
                group = values[begin[i] : end[i], j]
                result[i, j] = numba_func(group, group_index, *args)
        return result

    return group_agg


def generate_numba_transform_func(
    args: Tuple,
    kwargs: Dict[str, Any],
    func: Callable[..., Scalar],
    engine_kwargs: Optional[Dict[str, bool]],
) -> Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int], np.ndarray]:
    """
    Generate a numba jitted transform function specified by values from engine_kwargs.

    1. jit the user's function
    2. Return a groupby agg function with the jitted function inline

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

    validate_udf(func)

    numba_func = jit_user_function(func, nopython, nogil, parallel)

    numba = import_optional_dependency("numba")

    if parallel:
        loop_range = numba.prange
    else:
        loop_range = range

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def group_transform(
        values: np.ndarray,
        index: np.ndarray,
        begin: np.ndarray,
        end: np.ndarray,
        num_groups: int,
        num_columns: int,
    ) -> np.ndarray:
        result = np.empty((len(values), num_columns))
        for i in loop_range(num_groups):
            group_index = index[begin[i] : end[i]]
            for j in loop_range(num_columns):
                group = values[begin[i] : end[i], j]
                result[begin[i] : end[i], j] = numba_func(group, group_index, *args)
        return result

    return group_transform
