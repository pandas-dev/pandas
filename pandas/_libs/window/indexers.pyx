# cython: boundscheck=False, wraparound=False, cdivision=True

from typing import Optional, Tuple

import numpy as np
from numpy cimport ndarray, int64_t

# ----------------------------------------------------------------------
# The indexer objects for rolling
# These define start/end indexers to compute offsets


class BaseIndexer:
    """Base class for window bounds calculations"""

    def __init__(
        self,
        **kwargs,
    ):
        """
        Parameters
        ----------
        **kwargs :
            keyword argument that will be available when get_window_bounds is called
        """
        self.__dict__.update(kwargs)

    def get_window_bounds(
        self,
        num_values: int = 0,
        window_size: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
        win_type: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the bounds of a window.

        Parameters
        ----------
        num_values : int, default 0
            number of values that will be aggregated over
        window_size : int, default 0
            the number of rows in a window
        min_periods : int, default None
            min_periods passed from the top level rolling API
        center : bool, default None
            center passed from the top level rolling API
        closed : str, default None
            closed passed from the top level rolling API
        win_type : str, default None
            win_type passed from the top level rolling API

        Returns
        -------
        A tuple of ndarray[int64]s, indicating the boundaries of each
        window
        """
        raise NotImplementedError


class FixedWindowIndexer(BaseIndexer):
    """Creates window boundaries that are of fixed length."""

    def get_window_bounds(self,
        num_values: int = 0,
        window_size: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
        win_type: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the fixed bounds of a window.

        Parameters
        ----------
        num_values : int, default 0
            number of values that will be aggregated over
        window_size : int, default 0
            the number of rows in a window
        min_periods : int, default None
            min_periods passed from the top level rolling API
        center : bool, default None
            center passed from the top level rolling API
        closed : str, default None
            closed passed from the top level rolling API
        win_type : str, default None
            win_type passed from the top level rolling API

        Returns
        -------
        A tuple of ndarray[int64]s, indicating the boundaries of each
        window
        """
        cdef:
            ndarray[int64_t, ndim=1] start, start_s, start_e, end, end_s, end_e

        start_s = np.zeros(window_size, dtype='int64')
        start_e = np.arange(window_size, num_values, dtype='int64') - window_size + 1
        start = np.concatenate([start_s, start_e])[:num_values]

        end_s = np.arange(window_size, dtype='int64') + 1
        end_e = start_e + window_size
        end = np.concatenate([end_s, end_e])[:num_values]
        return start, end


class VariableWindowIndexer(BaseIndexer):
    """Creates window boundaries that are of variable length, namely for time series."""

    @staticmethod
    def _get_window_bound(
        int64_t num_values,
        int64_t window_size,
        object min_periods,
        object center,
        object closed,
        object win_type,
        const int64_t[:] index
    ):
        cdef:
            bint left_closed = False
            bint right_closed = False
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
                start_bound = index[i] - window_size

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

    def get_window_bounds(self,
        num_values: int = 0,
        window_size: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
        win_type: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the variable bounds of a window.

        Parameters
        ----------
        num_values : int, default 0
            number of values that will be aggregated over
        window_size : int, default 0
            the number of rows in a window
        min_periods : int, default None
            min_periods passed from the top level rolling API
        center : bool, default None
            center passed from the top level rolling API
        closed : str, default None
            closed passed from the top level rolling API
        win_type : str, default None
            win_type passed from the top level rolling API

        Returns
        -------
        A tuple of ndarray[int64]s, indicating the boundaries of each
        window
        """
        # We do this since cython doesn't like accessing class attributes in nogil
        return self._get_window_bound(
            num_values, window_size, min_periods, center, closed, win_type, self.index
        )
