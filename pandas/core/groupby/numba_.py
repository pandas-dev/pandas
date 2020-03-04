import inspect
import types

import numpy as np

from pandas.compat._optional import import_optional_dependency


class InvalidApply(Exception):
    pass


def execute_groupby_function(splitter, f):
    """Mimics apply_frame_axis0 which is the Cython equivalent of this function."""
    results = []
    for _, group in splitter:
        # TODO: what about series names/dataframe columns
        index = group.index
        values_as_array = group.to_numpy()
        index_as_array = index.to_numpy()
        try:
            # TODO: support *args, **kwargs here
            group_result = f(values_as_array, index_as_array)
        except Exception:
            # We can't be more specific without knowing something about `f`
            # Like we do in Cython
            raise InvalidApply("Let this error raise above us")
        # Reconstruct the pandas object (expected downstream)
        # This construction will fail is there is mutation,
        # but we're banning it with numba?
        group_result = group._constructor(group_result, index=index)
        results.append(group_result)

    return results


def validate_apply_function_signature(func):
    """
    Validate that the apply function's first 2 arguments are 'values' and 'index'.

    func : function
        function to be applied to each group and will be JITed
    """
    apply_function_signature = list(inspect.signature(func).parameters.keys())[:2]
    if apply_function_signature != ["values", "index"]:
        raise ValueError(
            "The apply function's first 2 arguments must be 'values' and 'index'"
        )


def make_groupby_apply(
    func, args, nogil, parallel, nopython,
):
    """
    Creates a JITted groupby apply function with a JITted version of
    the user's function.

    Parameters
    ----------
    func : function
        function to be applied to each group and will be JITed
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

    if isinstance(func, numba.targets.registry.CPUDispatcher):
        # Don't jit a user passed jitted function
        numba_func = func
    else:

        @numba.generated_jit(nopython=nopython, nogil=nogil, parallel=parallel)
        def numba_func(group, *_args):
            if getattr(np, func.__name__, False) is func or isinstance(
                func, types.BuiltinFunctionType
            ):
                jf = func
            else:
                jf = numba.jit(func, nopython=nopython, nogil=nogil)

            def impl(group, *_args):
                return jf(group, *_args)

            return impl

    return numba_func


def generate_numba_apply_func(
    args, kwargs, func, engine_kwargs,
):
    """
    Generate a numba jitted apply function specified by values from engine_kwargs.

    1. jit the user's function

    Configurations specified in engine_kwargs apply to both the user's
    function _AND_ the rolling apply function.

    Parameters
    ----------
    args : tuple
        *args to be passed into the function
    kwargs : dict
        **kwargs to be passed into the function
    func : function
        function to be applied to each group and will be JITed
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

    return make_groupby_apply(func, args, nogil, parallel, nopython)
