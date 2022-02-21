from datetime import timedelta
import numbers
from typing import (
    Any,
    Generic,
    TypeVar,
    overload,
)

import numpy as np
import numpy.typing as npt

from pandas._libs import (
    Timedelta,
    Timestamp,
)
from pandas._typing import IntervalBound

_OrderableMixinT = TypeVar(
    "_OrderableMixinT", int, float, Timestamp, Timedelta, npt.NDArray[np.generic]
)
_OrderableT = TypeVar("_OrderableT", int, float, Timestamp, Timedelta)

# note: mypy doesn't support overloading properties
# based on github.com/microsoft/python-type-stubs/pull/167
class _LengthProperty:
    @overload
    def __get__(self, instance: IntervalMixin[Timestamp], owner: Any) -> Timedelta: ...
    @overload
    def __get__(
        self, instance: IntervalMixin[_OrderableMixinT], owner: Any
    ) -> _OrderableMixinT: ...

class IntervalMixin(Generic[_OrderableMixinT]):
    @property
    def closed_left(self) -> bool: ...
    @property
    def closed_right(self) -> bool: ...
    @property
    def open_left(self) -> bool: ...
    @property
    def open_right(self) -> bool: ...
    @property
    def mid(self) -> _OrderableT: ...
    length: _LengthProperty
    @property
    def is_empty(self) -> bool: ...
    def _check_closed_matches(self, other: IntervalMixin, name: str = ...) -> None: ...

class Interval(IntervalMixin[_OrderableT]):
    def __init__(
        self,
        left: _OrderableT,
        right: _OrderableT,
        closed: IntervalBound = ...,
    ) -> None: ...
    @property
    def closed(self) -> str: ...
    @property
    def left(self) -> _OrderableT: ...
    @property
    def right(self) -> _OrderableT: ...
    def __str__(self) -> str: ...
    # TODO: could return Interval with different type
    def __add__(
        self, y: numbers.Number | np.timedelta64 | timedelta
    ) -> Interval[_OrderableT]: ...
    def __radd__(
        self, y: numbers.Number | np.timedelta64 | timedelta
    ) -> Interval[_OrderableT]: ...
    def __sub__(
        self, y: numbers.Number | np.timedelta64 | timedelta
    ) -> Interval[_OrderableT]: ...
    def __mul__(self, y: numbers.Number) -> Interval[_OrderableT]: ...
    def __rmul__(self, y: numbers.Number) -> Interval[_OrderableT]: ...
    def __truediv__(self, y: numbers.Number) -> Interval[_OrderableT]: ...
    def __floordiv__(self, y: numbers.Number) -> Interval[_OrderableT]: ...
    def __hash__(self) -> int: ...
    def __contains__(self: Interval[_OrderableT], key: _OrderableT) -> bool: ...
    def overlaps(self, other: Interval[_OrderableT]) -> bool: ...

VALID_CLOSED: frozenset[str]

# takes npt.NDArray[Interval[_OrderableT]] and returns arrays of type
# _OrderableT but _Orderable is not a valid dtype
def intervals_to_interval_bounds(
    intervals: npt.NDArray[np.object_], validate_closed: bool = ...
) -> tuple[np.ndarray, np.ndarray, str]: ...

# from pandas/_libs/intervaltree.pxi.in
_GenericT = TypeVar("_GenericT", bound=np.generic)

# error: Value of type variable "_OrderableMixinT" of "IntervalMixin"
# cannot be "ndarray"
class IntervalTree(
    Generic[_GenericT],
    IntervalMixin[npt.NDArray[_GenericT]],  # type: ignore[type-var]
):
    _na_count: int
    def __init__(
        self,
        left: npt.NDArray[_GenericT],
        right: npt.NDArray[_GenericT],
        closed: IntervalBound = ...,
        leaf_size: int = ...,
    ) -> None: ...
    @property
    def left_sorter(self) -> npt.NDArray[_GenericT]: ...
    @property
    def right_sorter(self) -> npt.NDArray[_GenericT]: ...
    @property
    def is_overlapping(self) -> bool: ...
    @property
    def is_monotonic_increasing(self) -> bool: ...
    def get_indexer(self, target: np.ndarray) -> npt.NDArray[np.intp]: ...
    def get_indexer_non_unique(
        self, target: np.ndarray
    ) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp]]: ...
    def __repr__(self) -> str: ...
    def clear_mapping(self) -> None: ...
