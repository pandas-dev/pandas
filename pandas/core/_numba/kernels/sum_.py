"""
Numba 1D sum kernels that can be shared by
* Dataframe / Series
* groupby
* rolling / expanding

Mirrors pandas/_libs/window/aggregation.pyx
"""
from __future__ import annotations

import numba
import numpy as np

from pandas.core._numba.kernels.shared import is_monotonic_increasing


@numba.jit(nopython=True, nogil=True, parallel=False)
def add_sum(
    val: float,
    nobs: int,
    sum_x: float,
    compensation: float,
    num_consecutive_same_value: int,
    prev_value: float,
) -> tuple[int, float, float, int, float]:
    if not np.isnan(val):
        nobs += 1
        y = val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t

        if val == prev_value:
            num_consecutive_same_value += 1
        else:
            num_consecutive_same_value = 1
        prev_value = val

    return nobs, sum_x, compensation, num_consecutive_same_value, prev_value


@numba.jit(nopython=True, nogil=True, parallel=False)
def remove_sum(
    val: float, nobs: int, sum_x: float, compensation: float
) -> tuple[int, float, float]:
    if not np.isnan(val):
        nobs -= 1
        y = -val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t
    return nobs, sum_x, compensation


def sliding_sum(values, start, end, min_periods):
    # Need a dummy function to overload
    pass


@numba.extending.overload(sliding_sum)
def sliding_sum_wrapper(
    values: np.ndarray, start: np.ndarray, end: np.ndarray, min_periods: int
) -> np.ndarray:
    dtype = values.dtype

    if isinstance(dtype, numba.types.Integer):
        na_val = 0

        @numba.extending.register_jitable(inline="always")
        def process_na_func(x, na_pos):
            return na_pos.append(x)

    else:
        na_val = np.nan

        @numba.extending.register_jitable(inline="always")
        def process_na_func(x, na_pos):
            pass

    @numba.extending.register_jitable
    def sliding_sum(
        values: np.ndarray,
        start: np.ndarray,
        end: np.ndarray,
        min_periods: int,
    ) -> np.ndarray:
        N = len(start)
        nobs = 0
        sum_x = 0
        compensation_add = 0
        compensation_remove = 0
        # Stores positions of the na_values
        # Trick to force list to be inferred as int list
        # https://numba.pydata.org/numba-doc/latest/user/troubleshoot.html#my-code-has-an-untyped-list-problem
        na_pos = [0 for i in range(0)]

        is_monotonic_increasing_bounds = is_monotonic_increasing(
            start
        ) and is_monotonic_increasing(end)

        # output = np.empty(N, dtype=dtype)
        output = np.empty(N, dtype=np.float64)

        for i in range(N):
            s = start[i]
            e = end[i]
            if i == 0 or not is_monotonic_increasing_bounds:
                prev_value = values[s]
                num_consecutive_same_value = 0

                for j in range(s, e):
                    val = values[j]
                    (
                        nobs,
                        sum_x,
                        compensation_add,
                        num_consecutive_same_value,
                        prev_value,
                    ) = add_sum(
                        val,
                        nobs,
                        sum_x,
                        compensation_add,
                        num_consecutive_same_value,
                        prev_value,
                    )
            else:
                for j in range(start[i - 1], s):
                    val = values[j]
                    nobs, sum_x, compensation_remove = remove_sum(
                        val, nobs, sum_x, compensation_remove
                    )

                for j in range(end[i - 1], e):
                    val = values[j]
                    (
                        nobs,
                        sum_x,
                        compensation_add,
                        num_consecutive_same_value,
                        prev_value,
                    ) = add_sum(
                        val,
                        nobs,
                        sum_x,
                        compensation_add,
                        num_consecutive_same_value,
                        prev_value,
                    )

            if nobs == 0 == min_periods:
                result = 0
            elif nobs >= min_periods:
                if num_consecutive_same_value >= nobs:
                    result = prev_value * nobs
                else:
                    result = sum_x
            else:
                result = na_val
                na_pos = process_na_func(i, na_pos)

            output[i] = result

            if not is_monotonic_increasing_bounds:
                nobs = 0
                sum_x = 0
                compensation_remove = 0

        return output, na_pos

    if isinstance(dtype, numba.types.Integer):

        def sliding_sum_with_nan(
            values: np.ndarray, start: np.ndarray, end: np.ndarray, min_periods: int
        ):
            output, na_pos = sliding_sum(values, start, end, min_periods)
            output = output.astype("float64")
            # Fancy indexing w/ list  doesn't yet work for numba
            # xref https://github.com/numba/numba/issues/8616
            # output[na_pos] = np.nan
            for pos in na_pos:
                output[pos] = np.nan
            return output

        return sliding_sum_with_nan
    else:

        def postprocess_sliding_sum(
            values: np.ndarray, start: np.ndarray, end: np.ndarray, min_periods: int
        ):
            return sliding_sum(values, start, end, min_periods)[0]

        return postprocess_sliding_sum
