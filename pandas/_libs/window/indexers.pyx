# cython: boundscheck=False, wraparound=False, cdivision=True

import numpy as np
from numpy cimport ndarray, int64_t

# Cython routines for window indexers


def calculate_variable_window_bounds(
    int64_t num_values,
    int64_t window_size,
    object min_periods,  # unused but here to match get_window_bounds signature
    object center,  # unused but here to match get_window_bounds signature
    object closed,
    const int64_t[:] index
):
    """
    Calculate window boundaries for rolling windows from a time offset.

    Parameters
    ----------
    num_values : int64
        total number of values

    window_size : int64
        window size calculated from the offset

    min_periods : object
        ignored, exists for compatibility

    center : object
        ignored, exists for compatibility

    closed : str
        string of side of the window that should be closed

    index : ndarray[int64]
        time series index to roll over

    Returns
    -------
    (ndarray[int64], ndarray[int64])
    """
    cdef:
        bint left_closed = False
        bint right_closed = False
        int index_growth_sign = 1
        ndarray[int64_t, ndim=1] start, end
        int64_t start_bound, end_bound
        Py_ssize_t i, j

    # if windows is variable, default is 'right', otherwise default is 'both'
    if closed is None:
        closed = 'right' if index is not None else 'both'

    if closed in ['right', 'both']:
        right_closed = True

    if closed in ['left', 'both']:
        left_closed = True

    if index[num_values - 1] < index[0]:
        index_growth_sign = -1

    start = np.empty(num_values, dtype='int64')
    start.fill(-1)
    end = np.empty(num_values, dtype='int64')
    end.fill(-1)

    start[0] = 0

    # right endpoint is closed
    if right_closed:
        end[0] = 1
    # right endpoint is open
    else:
        end[0] = 0

    with nogil:

        # start is start of slice interval (including)
        # end is end of slice interval (not including)
        for i in range(1, num_values):
            end_bound = index[i]
            start_bound = index[i] - index_growth_sign * window_size

            # left endpoint is closed
            if left_closed:
                start_bound -= 1

            # advance the start bound until we are
            # within the constraint
            start[i] = i
            for j in range(start[i - 1], i):
                if (index[j] - start_bound) * index_growth_sign > 0:
                    start[i] = j
                    break

            # end bound is previous end
            # or current index
            if (index[end[i - 1]] - end_bound) * index_growth_sign <= 0:
                end[i] = i + 1
            else:
                end[i] = end[i - 1]

            # right endpoint is open
            if not right_closed:
                end[i] -= 1
    return start, end
