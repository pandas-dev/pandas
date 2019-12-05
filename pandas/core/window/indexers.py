"""Indexer objects for computing start/end window bounds for rolling operations"""
from typing import Optional, Tuple

import numpy as np

from pandas._libs.window.indexers import calculate_variable_window_bounds
from pandas.util._decorators import Appender

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
    """Base class for window bounds calculations"""

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
