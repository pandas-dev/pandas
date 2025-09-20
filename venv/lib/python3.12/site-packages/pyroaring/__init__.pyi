import array
from typing import overload, TypedDict
from collections.abc import Iterable, Iterator

from typing_extensions import Self

__version__: str
__croaring_version__: str


class _Statistics(TypedDict):
    # Type as observed in the returned values.
    # Warning: This type does not exist at runtime.

    n_containers: int
    n_array_containers: int
    n_run_containers: int
    n_bitset_containers: int
    n_values_array_containers: int
    n_values_run_containers: int
    n_values_bitset_containers: int
    n_bytes_array_containers: int
    n_bytes_run_containers: int
    n_bytes_bitset_containers: int
    max_value: int
    min_value: int
    sum_value: int
    cardinality: int


class AbstractBitMap:
    def __init__(self, values: Iterable[int] | None = None, copy_on_write: bool = False, optimize: bool = True) -> None:
        ...

    @property
    def copy_on_write(self) -> bool:
        ...

    def run_optimize(self) -> bool:
        ...

    def shrink_to_fit(self) -> int:
        ...

    def __contains__(self, value: int) -> bool:
        ...

    def __bool__(self) -> bool:
        ...

    def __len__(self) -> int:
        ...

    def __lt__(self, other: AbstractBitMap) -> bool:
        ...

    def __le__(self, other: AbstractBitMap) -> bool:
        ...

    def __eq__(self, other: object) -> bool:
        ...

    def __ne__(self, other: object) -> bool:
        ...

    def __gt__(self, other: AbstractBitMap) -> bool:
        ...

    def __ge__(self, other: AbstractBitMap) -> bool:
        ...

    def contains_range(self, range_start: int, range_end: int) -> bool:
        ...

    def range_cardinality(self, range_start: int, range_end: int) -> int:
        ...

    def iter_equal_or_larger(self, val: int) -> Iterator[int]:
        ...

    def __iter__(self) -> Iterator[int]:
        ...

    def flip(self, start: int, end: int) -> Self:
        ...

    def shift(self, offset: int) -> Self:
        ...

    def copy(self) -> Self:
        ...

    def isdisjoint(self, other: AbstractBitMap) -> bool:
        ...

    def issubset(self, other: AbstractBitMap) -> bool:
        ...

    def issuperset(self, other: AbstractBitMap) -> bool:
        ...

    # Note: `difference` and others are sort-of set up like they're meant to be
    # static methods (accepting _only_ `*bitmaps` in the underlying Cython
    # code), however at runtime they require at least one argument and return an
    # instance of the same type as that value -- like instance methods. Typing
    # them as instances methods ensures that mypy matches this behaviour (other
    # type checkers untested), even when used statically as their docstrings
    # suggest.

    def difference(self, *bitmaps: AbstractBitMap) -> Self:
        ...

    def symmetric_difference(self, other: AbstractBitMap) -> Self:
        ...

    def union(self, *bitmaps: AbstractBitMap) -> Self:
        ...

    def intersection(self, *bitmaps: AbstractBitMap) -> Self:
        ...

    def __or__(self, other: AbstractBitMap) -> Self:
        ...

    def __and__(self, other: AbstractBitMap) -> Self:
        ...

    def __xor__(self, other: AbstractBitMap) -> Self:
        ...

    def __sub__(self, other: AbstractBitMap) -> Self:
        ...

    def union_cardinality(self, other: AbstractBitMap) -> int:
        ...

    def intersection_cardinality(self, other: AbstractBitMap) -> int:
        ...

    def difference_cardinality(self, other: AbstractBitMap) -> int:
        ...

    def symmetric_difference_cardinality(self, other: AbstractBitMap) -> int:
        ...

    def intersect(self, other: AbstractBitMap) -> bool:
        ...

    def jaccard_index(self, other: AbstractBitMap) -> float:
        ...

    def get_statistics(self) -> _Statistics:
        ...

    def min(self) -> int:
        ...

    def max(self) -> int:
        ...

    def rank(self, value: int) -> int:
        ...

    def next_set_bit(self, value: int) -> int:
        ...

    @overload
    def __getitem__(self, value: int) -> int:
        ...

    @overload
    def __getitem__(self, value: slice) -> Self:
        ...

    def serialize(self) -> bytes:
        ...

    @classmethod
    def deserialize(cls, buff: bytes) -> Self:
        ...

    def __getstate__(self) -> bytes:
        ...

    def __setstate__(self, state: bytes) -> Self:
        ...

    def __sizeof__(self) -> int:
        ...

    def to_array(self) -> array.array[int]:
        ...


class FrozenBitMap(AbstractBitMap):
    def __hash__(self) -> int:
        ...


class BitMap(AbstractBitMap):
    def add(self, value: int) -> None:
        ...

    def add_checked(self, value: int) -> None:
        ...

    def update(self, *all_values: Iterable[int]) -> None:
        ...

    def discard(self, value: int) -> None:
        ...

    def remove(self, value: int) -> None:
        ...

    def __ior__(self, other: AbstractBitMap) -> Self:
        ...

    def __iand__(self, other: AbstractBitMap) -> Self:
        ...

    def __ixor__(self, other: AbstractBitMap) -> Self:
        ...

    def __isub__(self, other: AbstractBitMap) -> Self:
        ...

    def intersection_update(self, *all_values: Iterable[int]) -> None:
        ...

    def difference_update(self, *others: AbstractBitMap) -> None:
        ...

    def symmetric_difference_update(self, other: AbstractBitMap) -> None:
        ...

    def overwrite(self, other: AbstractBitMap) -> None:
        ...

    def clear(self) -> None:
        ...

    def pop(self) -> int:
        ...

    def flip_inplace(self, start: int, end: int) -> None:
        ...

    def add_range(self, range_start: int, range_end: int) -> None:
        ...

    def remove_range(self, range_start: int, range_end: int) -> None:
        ...

class AbstractBitMap64:
    def __init__(self, values: Iterable[int] | None = None, copy_on_write: bool = False, optimize: bool = True) -> None:
        ...

    @property
    def copy_on_write(self) -> bool:
        ...

    def run_optimize(self) -> bool:
        ...

    def shrink_to_fit(self) -> int:
        ...

    def __contains__(self, value: int) -> bool:
        ...

    def __bool__(self) -> bool:
        ...

    def __len__(self) -> int:
        ...

    def __lt__(self, other: AbstractBitMap64) -> bool:
        ...

    def __le__(self, other: AbstractBitMap64) -> bool:
        ...

    def __eq__(self, other: object) -> bool:
        ...

    def __ne__(self, other: object) -> bool:
        ...

    def __gt__(self, other: AbstractBitMap64) -> bool:
        ...

    def __ge__(self, other: AbstractBitMap64) -> bool:
        ...

    def contains_range(self, range_start: int, range_end: int) -> bool:
        ...

    def range_cardinality(self, range_start: int, range_end: int) -> int:
        ...

    def iter_equal_or_larger(self, val: int) -> Iterator[int]:
        ...

    def __iter__(self) -> Iterator[int]:
        ...

    def flip(self, start: int, end: int) -> Self:
        ...

    def shift(self, offset: int) -> Self:
        ...

    def copy(self) -> Self:
        ...

    def isdisjoint(self, other: AbstractBitMap64) -> bool:
        ...

    def issubset(self, other: AbstractBitMap64) -> bool:
        ...

    def issuperset(self, other: AbstractBitMap64) -> bool:
        ...

    def difference(self, *bitmaps: AbstractBitMap64) -> Self:
        ...

    def symmetric_difference(self, other: AbstractBitMap64) -> Self:
        ...

    def union(self, *bitmaps: AbstractBitMap64) -> Self:
        ...

    def intersection(self, *bitmaps: AbstractBitMap64) -> Self:
        ...

    def __or__(self, other: AbstractBitMap64) -> Self:
        ...

    def __and__(self, other: AbstractBitMap64) -> Self:
        ...

    def __xor__(self, other: AbstractBitMap64) -> Self:
        ...

    def __sub__(self, other: AbstractBitMap64) -> Self:
        ...

    def union_cardinality(self, other: AbstractBitMap64) -> int:
        ...

    def intersection_cardinality(self, other: AbstractBitMap64) -> int:
        ...

    def difference_cardinality(self, other: AbstractBitMap64) -> int:
        ...

    def symmetric_difference_cardinality(self, other: AbstractBitMap64) -> int:
        ...

    def intersect(self, other: AbstractBitMap64) -> bool:
        ...

    def jaccard_index(self, other: AbstractBitMap64) -> float:
        ...

    def get_statistics(self) -> _Statistics:
        ...

    def min(self) -> int:
        ...

    def max(self) -> int:
        ...

    def rank(self, value: int) -> int:
        ...

    def next_set_bit(self, value: int) -> int:
        ...

    @overload
    def __getitem__(self, value: int) -> int:
        ...

    @overload
    def __getitem__(self, value: slice) -> Self:
        ...

    def serialize(self) -> bytes:
        ...

    @classmethod
    def deserialize(cls, buff: bytes) -> Self:
        ...

    def __getstate__(self) -> bytes:
        ...

    def __setstate__(self, state: bytes) -> Self:
        ...

    def __sizeof__(self) -> int:
        ...

    def to_array(self) -> array.array[int]:
        ...


class FrozenBitMap64(AbstractBitMap64):
    def __hash__(self) -> int:
        ...


class BitMap64(AbstractBitMap64):
    def add(self, value: int) -> None:
        ...

    def add_checked(self, value: int) -> None:
        ...

    def update(self, *all_values: Iterable[int]) -> None:
        ...

    def discard(self, value: int) -> None:
        ...

    def remove(self, value: int) -> None:
        ...

    def __ior__(self, other: AbstractBitMap64) -> Self:
        ...

    def __iand__(self, other: AbstractBitMap64) -> Self:
        ...

    def __ixor__(self, other: AbstractBitMap64) -> Self:
        ...

    def __isub__(self, other: AbstractBitMap64) -> Self:
        ...

    def intersection_update(self, *all_values: Iterable[int]) -> None:
        ...

    def difference_update(self, *others: AbstractBitMap64) -> None:
        ...

    def symmetric_difference_update(self, other: AbstractBitMap64) -> None:
        ...

    def overwrite(self, other: AbstractBitMap64) -> None:
        ...

    def clear(self) -> None:
        ...

    def pop(self) -> int:
        ...

    def flip_inplace(self, start: int, end: int) -> None:
        ...

    def add_range(self, range_start: int, range_end: int) -> None:
        ...

    def remove_range(self, range_start: int, range_end: int) -> None:
        ...