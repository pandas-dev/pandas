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
        result = np.empty((values.shape[0], len(start)), dtype=result_dtype)
        na_positions = {}
        for i in numba.prange(values.shape[0]):
            output, na_pos = func(
                values[i], result_dtype, start, end, min_periods, *args
            )
            result[i] = output
            if len(na_pos) > 0:
                na_positions[i] = np.array(na_pos)
        return result, na_positions

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
    np.dtype("complex64"): np.complex64,
    np.dtype("complex128"): np.complex128,
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

    # It also post-processes the values by inserting nans where number of observations
    # is less than min_periods
    # Cannot do this in numba nopython mode
    # (you'll run into type-unification error when you cast int -> float)
    def looper_wrapper(values, start, end, min_periods, **kwargs):
        result_dtype = dtype_mapping[values.dtype]
        column_looper = make_looper(func, result_dtype, nopython, nogil, parallel)
        # Need to unpack kwargs since numba only supports *args
        result, na_positions = column_looper(
            values, start, end, min_periods, *kwargs.values()
        )
        if result.dtype.kind == "i":
            # Look if na_positions is not empty
            # If so, convert the whole block
            # This is OK since int dtype cannot hold nan,
            # so if min_periods not satisfied for 1 col, it is not satisfied for
            # all columns at that index
            for na_pos in na_positions.values():
                if len(na_pos) > 0:
                    result = result.astype("float64")
                    break
        # TODO: Optimize this
        for i, na_pos in na_positions.items():
            if len(na_pos) > 0:
                result[i, na_pos] = np.nan
        return result

    return looper_wrapper
