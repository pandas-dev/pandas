"""
Numba 1D min/max kernels that can be shared by
* Dataframe / Series
* groupby
* rolling / expanding

Mirrors pandas/_libs/window/aggregation.pyx
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

import numba
import numpy as np

if TYPE_CHECKING:
    from pandas._typing import npt


@numba.njit(nogil=True, parallel=False)
def bisect_left(a: list[Any], x: Any, lo: int = 0, hi: int = -1) -> int:
    """Same as https://docs.python.org/3/library/bisect.html; not in numba yet!"""
    if hi == -1:
        hi = len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@numba.jit(nopython=True, nogil=True, parallel=False)
def sliding_min_max(
    values: np.ndarray,
    result_dtype: np.dtype,
    start: np.ndarray,
    end: np.ndarray,
    min_periods: int,
    is_max: bool,
) -> tuple[np.ndarray, list[int]]:
    # Basic idea of the algorithm: https://stackoverflow.com/a/12239580
    # It was generalized to work with an arbitrary list of any window size and position
    # by adding the Dominators stack.

    N = len(start)
    na_pos = []
    output = np.empty(N, dtype=result_dtype)

    def cmp(a: Any, b: Any, is_max: bool) -> bool:
        if is_max:
            return a >= b
        else:
            return a <= b

    # Indices of bounded extrema in `values`. `candidates[i]` is always increasing.
    # `values[candidates[i]]` is decreasing for max and increasing for min.
    candidates: list[int] = []  # this is a queue
    # Indices of largest windows that "cover" preceding windows.
    dominators: list[int] = []  # this is a stack

    if min_periods < 1:
        min_periods = 1

    if N > 2:
        i_next = N - 1  # equivalent to i_next = i+1 inside the loop
        for i in range(N - 2, -1, -1):
            next_dominates = start[i_next] < start[i]
            if next_dominates and (
                not dominators or start[dominators[-1]] > start[i_next]
            ):
                dominators.append(i_next)
            i_next = i

    # NaN tracking to guarantee min_periods
    valid_start = -min_periods

    last_end = 0
    last_start = -1

    for i in range(N):
        this_start = start[i].item()
        this_end = end[i].item()

        if dominators and dominators[-1] == i:
            dominators.pop()

        if not (
            this_end > last_end or (this_end == last_end and this_start >= last_start)
        ):
            raise ValueError(
                "Start/End ordering requirement is violated at index " + str(i)
            )

        stash_start = (
            this_start if not dominators else min(this_start, start[dominators[-1]])
        )
        while candidates and candidates[0] < stash_start:
            candidates.pop(0)

        for k in range(last_end, this_end):
            if not np.isnan(values[k]):
                valid_start += 1
                while valid_start >= 0 and np.isnan(values[valid_start]):
                    valid_start += 1
                while candidates and cmp(values[k], values[candidates[-1]], is_max):
                    candidates.pop()  # Q.pop_back()
                candidates.append(k)  # Q.push_back(k)

        if not candidates or (this_start > valid_start):
            if values.dtype.kind != "i":
                output[i] = np.nan
            else:
                na_pos.append(i)
        elif candidates[0] >= this_start:
            # ^^ This is here to avoid costly bisection for fixed window sizes.
            output[i] = values[candidates[0]]
        else:
            q_idx = bisect_left(candidates, this_start, lo=1)
            output[i] = values[candidates[q_idx]]
        last_end = this_end
        last_start = this_start

    return output, na_pos


@numba.jit(nopython=True, nogil=True, parallel=False)
def grouped_min_max(
    values: np.ndarray,
    result_dtype: np.dtype,
    labels: npt.NDArray[np.intp],
    ngroups: int,
    min_periods: int,
    is_max: bool,
    skipna: bool = True,
) -> tuple[np.ndarray, list[int]]:
    N = len(labels)
    nobs = np.zeros(ngroups, dtype=np.int64)
    na_pos = []
    output = np.empty(ngroups, dtype=result_dtype)

    for i in range(N):
        lab = labels[i]
        val = values[i]
        if lab < 0 or (not skipna and nobs[lab] >= 1 and np.isnan(output[lab])):
            continue

        if values.dtype.kind == "i" or not np.isnan(val):
            nobs[lab] += 1
        else:
            if not skipna:
                # If skipna is False and we encounter a NaN,
                # both min and max of the group will be NaN
                output[lab] = np.nan
            continue

        if nobs[lab] == 1:
            # First element in group, set output equal to this
            output[lab] = val
            continue

        if is_max:
            if val > output[lab]:
                output[lab] = val
        elif val < output[lab]:
            output[lab] = val

    # Set labels that don't satisfy min_periods as np.nan
    for lab, count in enumerate(nobs):
        if count < min_periods:
            na_pos.append(lab)

    return output, na_pos
