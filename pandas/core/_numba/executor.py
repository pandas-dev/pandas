from __future__ import annotations

import functools
from typing import (
    TYPE_CHECKING,
    Callable,
)

if TYPE_CHECKING:
    from pandas._typing import Scalar

import numpy as np

from pandas.compat._optional import import_optional_dependency


@functools.lru_cache(maxsize=None)
def make_looper(func, result_dtype, nopython, nogil, parallel):
    if TYPE_CHECKING:
        import numba
    else:
        numba = import_optional_dependency("numba")

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def column_looper(
        values: np.ndarray,
        start: np.ndarray,
        end: np.ndarray,
        min_periods: int,
        *args,
    ):
        result = np.empty((len(start), values.shape[0]), dtype=result_dtype)
        for i in numba.prange(values.shape[0]):
            result[:, i] = func(values[i], start, end, min_periods, *args)
        return result.T

    return column_looper


default_dtype_mapping = {
    np.dtype("int8"): np.int64,
    np.dtype("int16"): np.int64,
    np.dtype("int32"): np.int64,
    np.dtype("int64"): np.int64,
    np.dtype("uint8"): np.uint64,
    np.dtype("uint16"): np.uint64,
    np.dtype("uint32"): np.uint64,
    np.dtype("uint64"): np.uint64,
    np.dtype("float32"): np.float64,
    np.dtype("float64"): np.float64,
}


def generate_shared_aggregator(
    func: Callable[..., Scalar],
    dtype_mapping: dict[np.dtype, np.dtype],
    nopython: bool,
    nogil: bool,
    parallel: bool,
):
    """
    Generate a Numba function that loops over the columns 2D object and applies
    a 1D numba kernel over each column.

    Parameters
    ----------
    func : function
        aggregation function to be applied to each column
    dtype_mapping: dict or None
        If not None, maps a dtype to a result dtype.
        Otherwise, will fall back to default mapping.
        all integers -> int64
        all floats -> float64
    nopython : bool
        nopython to be passed into numba.jit
    nogil : bool
        nogil to be passed into numba.jit
    parallel : bool
        parallel to be passed into numba.jit

    Returns
    -------
    Numba function
    """

    # A wrapper around the looper function,
    # to dispatch based on dtype since numba is unable to do that in nopython mode
    def looper_wrapper(values, start, end, min_periods, **kwargs):
        result_dtype = dtype_mapping[values.dtype]
        column_looper = make_looper(func, result_dtype, nopython, nogil, parallel)
        # Need to unpack kwargs since numba only supports *args
        return column_looper(values, start, end, min_periods, *kwargs.values())

    return looper_wrapper
