import numpy as np

from pandas import (
    Timedelta,
    Timestamp,
)

VALID_CLOSED: frozenset[str]


class IntervalMixin:
    closed: str

    @property
    def closed_left(self) -> bool: ...

    @property
    def closed_right(self) -> bool: ...

    @property
    def open_left(self) -> bool: ...

    @property
    def open_right(self) -> bool: ...

    @property
    def mid(self): ...

    @property
    def length(self): ...

    @property
    def is_empty(self): ...

    def _check_closed_matches(self, other, name: str = ...) -> None: ...


class Interval(IntervalMixin):
    left: int | float | Timestamp | Timedelta
    right: int | float | Timestamp | Timedelta

    def __init__(self, left, right, closed: str = ...): ...

    def __contains__(self, key) -> bool: ...
    def __str__(self) -> str: ...
    def __add__(self, y): ...
    def __sub__(self, y): ...
    def __mul__(self, y): ...
    def __truediv__(self, y): ...
    def __floordiv__(self, y): ...

    def overlaps(self, other: Interval) -> bool: ...


def intervals_to_interval_bounds(
    intervals: np.ndarray,
    validate_closed: bool = ...,
) -> tuple[np.ndarray, np.ndarray, str]: ...


class IntervalTree(IntervalMixin):
    def __init__(self, left, right, closed=..., leaf_size=...): ...

    @property
    def left_sorter(self) -> np.ndarray: ...  # np.ndarray[np.intp]

    @property
    def right_sorter(self) -> np.ndarray: ...  # np.ndarray[np.intp]

    @property
    def is_overlapping(self) -> bool: ...

    @property
    def is_monotonic_increasing(self) -> bool: ...

    def get_indexer(
        self,
        target: np.ndarray,  # scalar_t[:]
    ) -> np.ndarray: ...  #  np.ndarray[np.intp]

    def get_indexer_non_unique(
        self,
        target: np.ndarray,  # scalar_t[:]
    ) -> tuple[
        np.ndarray,  # np.ndarray[np.intp]
        np.ndarray,  # np.ndarray[np.intp]
    ]: ...

    def clear_mapping(self) -> None: ...
