"""Indexer objects for computing start/end window bounds for rolling operations"""
from datetime import timedelta
from typing import Dict, Optional, Tuple, Type, Union

import numpy as np

from pandas._libs.window.indexers import calculate_variable_window_bounds
from pandas.util._decorators import Appender

from pandas.tseries.offsets import Nano

get_window_bounds_doc = """
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


class BaseIndexer:
    """Base class for window bounds calculations."""

    def __init__(
        self, index_array: Optional[np.ndarray] = None, window_size: int = 0, **kwargs,
    ):
        """
        Parameters
        ----------
        **kwargs :
            keyword arguments that will be available when get_window_bounds is called
        """
        self.index_array = index_array
        self.window_size = window_size
        # Set user defined kwargs as attributes that can be used in get_window_bounds
        for key, value in kwargs.items():
            setattr(self, key, value)

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        raise NotImplementedError


class FixedWindowIndexer(BaseIndexer):
    """Creates window boundaries that are of fixed length."""

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        start_s = np.zeros(self.window_size, dtype="int64")
        start_e = (
            np.arange(self.window_size, num_values, dtype="int64")
            - self.window_size
            + 1
        )
        start = np.concatenate([start_s, start_e])[:num_values]

        end_s = np.arange(self.window_size, dtype="int64") + 1
        end_e = start_e + self.window_size
        end = np.concatenate([end_s, end_e])[:num_values]
        return start, end


class VariableWindowIndexer(BaseIndexer):
    """Creates window boundaries that are of variable length, namely for time series."""

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        return calculate_variable_window_bounds(
            num_values, self.window_size, min_periods, center, closed, self.index_array,
        )


class VariableOffsetWindowIndexer(BaseIndexer):
    """Calculate window boundaries based on a non-fixed offset such as a BusinessDay"""

    def __init__(
        self,
        index_array: Optional[np.ndarray] = None,
        window_size: int = 0,
        index=None,
        offset=None,
        **kwargs,
    ):
        super().__init__(index_array, window_size, **kwargs)
        self.index = index
        self.offset = offset

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        # if windows is variable, default is 'right', otherwise default is 'both'
        if closed is None:
            closed = "right" if self.index is not None else "both"

        right_closed = closed in ["right", "both"]
        left_closed = closed in ["left", "both"]

        if self.index[num_values - 1] < self.index[0]:
            index_growth_sign = -1
        else:
            index_growth_sign = 1

        start = np.empty(num_values, dtype="int64")
        start.fill(-1)
        end = np.empty(num_values, dtype="int64")
        end.fill(-1)

        start[0] = 0

        # right endpoint is closed
        if right_closed:
            end[0] = 1
        # right endpoint is open
        else:
            end[0] = 0

        # start is start of slice interval (including)
        # end is end of slice interval (not including)
        for i in range(1, num_values):
            end_bound = self.index[i]
            start_bound = self.index[i] - index_growth_sign * self.offset

            # left endpoint is closed
            if left_closed:
                start_bound -= Nano(1)

            # advance the start bound until we are
            # within the constraint
            start[i] = i
            for j in range(start[i - 1], i):
                if (self.index[j] - start_bound) * index_growth_sign > timedelta(0):
                    start[i] = j
                    break

            # end bound is previous end
            # or current index
            if (self.index[end[i - 1]] - end_bound) * index_growth_sign <= timedelta(0):
                end[i] = i + 1
            else:
                end[i] = end[i - 1]

            # right endpoint is open
            if not right_closed:
                end[i] -= 1

        return start, end


class ExpandingIndexer(BaseIndexer):
    """Calculate expanding window bounds, mimicking df.expanding()"""

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        return (
            np.zeros(num_values, dtype=np.int64),
            np.arange(1, num_values + 1, dtype=np.int64),
        )


class FixedForwardWindowIndexer(BaseIndexer):
    """
    Creates window boundaries for fixed-length windows that include the
    current row.

    Examples
    --------
    >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
    >>> df
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    >>> indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=2)
    >>> df.rolling(window=indexer, min_periods=1).sum()
         B
    0  1.0
    1  3.0
    2  2.0
    3  4.0
    4  4.0
    """

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:

        if center:
            raise ValueError("Forward-looking windows can't have center=True")
        if closed is not None:
            raise ValueError(
                "Forward-looking windows don't support setting the closed argument"
            )

        start = np.arange(num_values, dtype="int64")
        end_s = start[: -self.window_size] + self.window_size
        end_e = np.full(self.window_size, num_values, dtype="int64")
        end = np.concatenate([end_s, end_e])

        return start, end


class GroupbyRollingIndexer(BaseIndexer):
    """Calculate bounds to compute groupby rolling, mimicking df.groupby().rolling()"""

    def __init__(
        self,
        index_array: Optional[np.ndarray],
        window_size: int,
        groupby_indicies: Dict,
        rolling_indexer: Union[Type[FixedWindowIndexer], Type[VariableWindowIndexer]],
        **kwargs,
    ):
        """
        Parameters
        ----------
        **kwargs :
            keyword arguments that will be available when get_window_bounds is called
        """
        self.groupby_indicies = groupby_indicies
        self.rolling_indexer = rolling_indexer
        super().__init__(index_array, window_size, **kwargs)

    @Appender(get_window_bounds_doc)
    def get_window_bounds(
        self,
        num_values: int = 0,
        min_periods: Optional[int] = None,
        center: Optional[bool] = None,
        closed: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        # 1) For each group, get the indices that belong to the group
        # 2) Use the indices to calculate the start & end bounds of the window
        # 3) Append the window bounds in group order
        start_arrays = []
        end_arrays = []
        window_indicies_start = 0
        for key, indicies in self.groupby_indicies.items():
            if self.index_array is not None:
                index_array = self.index_array.take(indicies)
            else:
                index_array = self.index_array
            indexer = self.rolling_indexer(
                index_array=index_array, window_size=self.window_size,
            )
            start, end = indexer.get_window_bounds(
                len(indicies), min_periods, center, closed
            )
            start = start.astype(np.int64)
            end = end.astype(np.int64)
            # Cannot use groupby_indicies as they might not be monotonic with the object
            # we're rolling over
            window_indicies = np.arange(
                window_indicies_start, window_indicies_start + len(indicies),
            )
            window_indicies_start += len(indicies)
            # Extend as we'll be slicing window like [start, end)
            window_indicies = np.append(
                window_indicies, [window_indicies[-1] + 1]
            ).astype(np.int64)
            start_arrays.append(window_indicies.take(start))
            end_arrays.append(window_indicies.take(end))
        start = np.concatenate(start_arrays)
        end = np.concatenate(end_arrays)
        return start, end
