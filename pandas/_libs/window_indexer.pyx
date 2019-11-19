# cython: boundscheck=False, wraparound=False, cdivision=True

import numpy as np
from numpy cimport ndarray, int64_t

# ----------------------------------------------------------------------
# The indexer objects for rolling
# These define start/end indexers to compute offsets


class MockFixedWindowIndexer:
    """

    We are just checking parameters of the indexer,
    and returning a consistent API with fixed/variable
    indexers.

    Parameters
    ----------
    values: ndarray
        values data array
    win: int64_t
        window size
    index: object
        index of the values
    closed: string
        closed behavior
    """
    def __init__(self, ndarray values, int64_t win, object closed, object index=None):

        self.start = np.empty(0, dtype='int64')
        self.end = np.empty(0, dtype='int64')

    def get_window_bounds(self):
        return self.start, self.end


class FixedWindowIndexer:
    """
    create a fixed length window indexer object
    that has start & end, that point to offsets in
    the index object; these are defined based on the win
    arguments

    Parameters
    ----------
    values: ndarray
        values data array
    win: int64_t
        window size
    index: object
        index of the values
    closed: string
        closed behavior
    """
    def __init__(self, ndarray values, int64_t win, object closed, object index=None):
        cdef:
            ndarray start_s, start_e, end_s, end_e
            int64_t N = len(values)

        start_s = np.zeros(win, dtype='int64')
        start_e = np.arange(win, N, dtype='int64') - win + 1
        self.start = np.concatenate([start_s, start_e])[:N]

        end_s = np.arange(win, dtype='int64') + 1
        end_e = start_e + win
        self.end = np.concatenate([end_s, end_e])[:N]

    def get_window_bounds(self):
        return self.start, self.end

class VariableWindowIndexer:
    """
    create a variable length window indexer object
    that has start & end, that point to offsets in
    the index object; these are defined based on the win
    arguments

    Parameters
    ----------
    values: ndarray
        values data array
    win: int64_t
        window size
    index: ndarray
        index of the values
    closed: string
        closed behavior
    """
    def __init__(self, ndarray values, int64_t win, object closed, ndarray index):
        cdef:
            bint left_closed = False
            bint right_closed = False
            int64_t N = len(index)

        # if windows is variable, default is 'right', otherwise default is 'both'
        if closed is None:
            closed = 'right' if index is not None else 'both'

        if closed in ['right', 'both']:
            right_closed = True

        if closed in ['left', 'both']:
            left_closed = True

        self.start, self.end = self.build(index, win, left_closed, right_closed, N)

        # TODO: Maybe will need to use this?
        # max window size
        #self.win = (self.end - self.start).max()

    def build(self, const int64_t[:] index, int64_t win, bint left_closed,
              bint right_closed, int64_t N):

        cdef:
            ndarray[int64_t] start, end
            int64_t start_bound, end_bound
            Py_ssize_t i, j

        start = np.empty(N, dtype='int64')
        start.fill(-1)
        end = np.empty(N, dtype='int64')
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
            for i in range(1, N):
                end_bound = index[i]
                start_bound = index[i] - win

                # left endpoint is closed
                if left_closed:
                    start_bound -= 1

                # advance the start bound until we are
                # within the constraint
                start[i] = i
                for j in range(start[i - 1], i):
                    if index[j] > start_bound:
                        start[i] = j
                        break

                # end bound is previous end
                # or current index
                if index[end[i - 1]] <= end_bound:
                    end[i] = i + 1
                else:
                    end[i] = end[i - 1]

                # right endpoint is open
                if not right_closed:
                    end[i] -= 1
        return start, end

    def get_window_bounds(self):
        return self.start, self.end