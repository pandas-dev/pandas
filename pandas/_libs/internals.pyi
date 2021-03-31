from typing import (
    Iterator,
    Sequence,
    overload,
)

import numpy as np

from pandas._typing import ArrayLike

def slice_len(slc: slice, objlen: int = ...) -> int: ...


def get_blkno_indexers(
    blknos: np.ndarray,  # int64_t[:]
    group: bool = ...,
) -> list[tuple[int, slice | np.ndarray]]: ...


def get_blkno_placements(
    blknos: np.ndarray,
    group: bool = ...,
) -> Iterator[tuple[int, BlockPlacement]]: ...


class BlockPlacement:
    def __init__(self, val: int | slice | np.ndarray): ...

    @property
    def indexer(self) -> np.ndarray | slice: ...

    @property
    def as_array(self) -> np.ndarray: ...

    @property
    def is_slice_like(self) -> bool: ...

    @overload
    def __getitem__(self, loc: slice | Sequence[int]) -> BlockPlacement: ...

    @overload
    def __getitem__(self, loc: int) -> int: ...

    def __iter__(self) -> Iterator[int]: ...

    def __len__(self) -> int: ...

    def delete(self, loc) -> BlockPlacement: ...

    def append(self, others: list[BlockPlacement]) -> BlockPlacement: ...


class Block:
    _mgr_locs: BlockPlacement
    ndim: int
    values: ArrayLike

    def __init__(self, values: ArrayLike, placement: BlockPlacement, ndim: int): ...
