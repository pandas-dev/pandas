# cython: boundscheck=False, wraparound=False, cdivision=True

import numpy as np

from numpy cimport (
    int64_t,
    ndarray,
)

# Cython routines for window indexers


def calculate_variable_window_bounds(
    int64_t num_values,
    int64_t window_size,
    object min_periods,  # unused but here to match get_window_bounds signature
    int64_t step_size,
    bint center,
    str closed,
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

    step_size : int64
        (minimum) step size of the moving window

    center : bint
        center the rolling window on the current observation

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
        ndarray[int64_t, ndim=1] start, end
        int64_t start_bound, end_bound, step_bound, index_growth_sign = 1
        Py_ssize_t i, j, last_valid

    # default is 'right'
    if closed is None:
        closed = 'right'

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
    if center:
        end_bound = index[0] + index_growth_sign * window_size / 2
        for j in range(0, num_values):
            if (index[j] - end_bound) * index_growth_sign < 0:
                end[0] = j + 1
            elif (index[j] - end_bound) * index_growth_sign == 0 and right_closed:
                end[0] = j + 1
            elif (index[j] - end_bound) * index_growth_sign >= 0:
                end[0] = j
                break

    step_bound = index[0] + index_growth_sign * step_size
    last_valid = 0

    with nogil:

        # start is start of slice interval (including)
        # end is end of slice interval (not including)
        for i in range(1, num_values):

            if index[i] - index_growth_sign * step_bound < 0:
                start[i] = end[i-1]
                end[i] = end[i-1]
                continue

            else:
                step_bound = index[i] + step_size

            if center:
                end_bound = index[i] + index_growth_sign * window_size / 2
                start_bound = index[i] - index_growth_sign * window_size / 2
            else:
                end_bound = index[i]
                start_bound = index[i] - index_growth_sign * window_size

            # left endpoint is closed
            if left_closed:
                start_bound -= 1 * index_growth_sign

            # advance the start bound until we are
            # within the constraint
            start[i] = i
            for j in range(start[last_valid], i):
                if (index[j] - start_bound) * index_growth_sign > 0:
                    start[i] = j
                    break

            # for centered window advance the end bound until we are
            # outside the constraint
            end[i] = num_values
            for j in range(max(end[last_valid], i), num_values):
                if (index[j] - end_bound) * index_growth_sign == 0 and right_closed:
                    end[i] = j + 1
                elif (index[j] - end_bound) * index_growth_sign >= 0:
                    end[i] = j
                    break

            last_valid = i

    return start, end
