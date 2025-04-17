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
    if hi == -1:
        hi = len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


# *** Notations ***
# N: size of the values[] array.
# NN: last accessed element (1-based) in the values[] array, that is max_<i>(end[<i>]).
#    In most cases NN==N (unless you are using a custom window indexer).
# M: number of min/max "jobs", that is, size of start[] and end[] arrays.
#    In pandas' common usage, M==N==NN, but it does not have to!
# k: maximum window size, that is max<i>(end[<i>] - start[<i>])
#
# *** Complexity ***
# - O(max(NN,M)) for constant window sizes.
# - O(max(NN,M)*log(k)) for arbitrary window sizes.
#
# *** Assumptions ***
# The min/max "jobs" have to be ordered in the lexiographical (end[i], start[i]) order.
# In the regular pandas' case with constant window size these array ARE properly
# sorted by construction.
# In other words, for each i = 0..N-2, this condition must hold:
# - (end[i+1] > end[i]) OR
# - (end[i+1] == end[i] AND start[i+1] >= start[i])
#
# To debug this with your favorite Python debugger:
# - Comment out the "numba.jit" line right below this comment above the function def.
# - Find and comment out a similar line in column_looper() defined in
#   generate_apply_looper() in executor.py.
# - Place a breakpoint in this function. Your Python debugger will stop there!
# - (Debugging was tested with VSCode on WSL.)
@numba.jit(nopython=True, nogil=True, parallel=False)
def sliding_min_max(
    values: np.ndarray,
    result_dtype: np.dtype,
    start: np.ndarray,
    end: np.ndarray,
    min_periods: int,
    is_max: bool,
) -> tuple[np.ndarray, list[int]]:
    # numba-only init part
    N = len(start)
    na_pos = []
    output = np.empty(N, dtype=result_dtype)

    def cmp(a: Any, b: Any, is_max: bool) -> bool:
        if is_max:
            return a >= b
        else:
            return a <= b

    # All comments below are for the case of a maximum, that is, is_max = True.
    # I will call this Q a "stash": preliminary calculations minimally necessary to
    # finish the job. Q will always be in ascending order regardless of min/max:
    # these are indices into the "values" array. The values[Q[i]], however, will
    # be in non-descending order for max, and non-ascending for min.
    # Think of this deque as indices of maximums found so far for varying window
    # positions. That is, there is only one maximum in the source array, but it may
    # not fit each window. So there are many "secondary maximums", and each of them
    # may potentially fit multiple windows (unless, of course you get a very special
    # case of an arary of strictly descending values and constant rolling window size).
    # In that case Q will be the longest, so here is an idea for the worst case
    # scenario testing.

    # We have to pretend, given that Numba has neither queue nor stack.
    Q: list = []  # this is a queue
    Dominators: list = []  # this is a stack
    # END-OF numba-only init part

    # Basic idea of the algorithm: https://stackoverflow.com/a/12239580
    # It was generalized to work with an arbitrary list of any window size and position

    # Zero is apparently passed here as a default.
    # It is important to have this value precise for calculations.
    if min_periods < 1:
        min_periods = 1

    # We will say that node j "dominates" node i if j comes after i, yet requires a
    # deeper deque Q at the time node i is processed in order to be able to finish
    # the job for node j. This precisely means the following two conditions:
    # - j > i
    # - start[j] < start[i].
    # We keep track of such nodes in the Dominators queue.
    # In addition, if it so happens that two nodes j1 and j2 dominate another node,
    # and j2 > j1, yet start[j2] <= start[j1], then we only need to keep track of j2.
    # (To understand why this works, note that the algorithm requires that
    # the "end" array is sorted in non-descending order, among other things.)
    if N > 2:
        i_next = N - 1  # equivalent to i_next = i+1 inside the loop
        for i in range(N - 2, -1, -1):
            next_dominates = start[i_next] < start[i]
            if next_dominates and (
                not Dominators or start[Dominators[-1]] > start[i_next]
            ):
                # Both ">" and ">=" would have been (logically) equivalent, but we are
                # shooting for the shortest size of the Dominators list, hence the
                # usage of ">"
                Dominators.append(i_next)
            i_next = i

    valid_start = -min_periods

    # Having these relieves us from having "if i>0" on each iteration for special
    # handling
    last_end = 0
    last_start = -1

    for i in range(N):
        this_start = start[i].item()
        this_end = end[i].item()

        if Dominators and Dominators[-1] == i:
            Dominators.pop()

        # TODO: Arguably there are benefits to having this consistency check before
        # this function is even called (e.g. in rolling.py).
        # Given the current implementation, it will be rather tricky at the moment
        # to have this check in rolling.py. Additionally, this is only required for
        # min/max, and may only ever be violated if user-defined window indexer is
        # used. Thus this is the best spot for it, given the circumstances.
        if not (
            this_end > last_end or (this_end == last_end and this_start >= last_start)
        ):
            raise ValueError(
                "Start/End ordering requirement is violated at index " + str(i)
            )

        # This is the least restrictive starting index that will take care of current
        # item (i) and all remaining items
        stash_start = (
            this_start if not Dominators else min(this_start, start[Dominators[-1]])
        )
        # Discard entries outside of the "needed" window. Do it first as to keep the
        # stash small.
        while Q and Q[0] < stash_start:
            Q.pop(0)

        for k in range(last_end, this_end):
            if not np.isnan(values[k]):
                valid_start += 1
                while valid_start >= 0 and np.isnan(values[valid_start]):
                    valid_start += 1
                while Q and cmp(values[k], values[Q[-1]], is_max):
                    Q.pop()  # Q.pop_back()
                Q.append(k)  # Q.push_back(k)

        if not Q or (this_start > valid_start):
            # The "not Q" condition means we have not seen anything but NaNs yet in
            # values[:this_end-1]. The "this_start > valid_start" condition means we
            # have not accumulated enough (min_periods or more) non-NaN values.
            if values.dtype.kind != "i":
                output[i] = np.nan
            else:
                na_pos.append(i)
        elif Q[0] >= this_start:
            # This is the only read-from-the-stash scenario that ever happens when
            # window size is constant across the set. This is also likely 99+% of
            # all use cases, thus we want to make sure we do not go into bisection
            # as to incur neither the *log(k) penalty nor the function call penalty
            # for this very common case. If we are here, then our stash is as "deep"
            # as what the current node ("job") requires. Thus take the front item.
            output[i] = values[Q[0]]
        else:
            # In this case our stash is bigger than what is necessary to compute this
            # node's output due to a wider search window at (one of) the nodes that
            # follow. We have to locate our value in the middle of the stash.
            # Since our stash is sorted, we can use binary search:
            # here we need to output the item closest to the front (idx=0) of the
            # stash that fits our window bounds. Item 0 has been looked at (and
            # discarded) by now, so lo=1
            q_idx = bisect_left(Q, this_start, lo=1)
            output[i] = values[Q[q_idx]]
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
        else:
            if val < output[lab]:
                output[lab] = val

    # Set labels that don't satisfy min_periods as np.nan
    for lab, count in enumerate(nobs):
        if count < min_periods:
            na_pos.append(lab)

    return output, na_pos
