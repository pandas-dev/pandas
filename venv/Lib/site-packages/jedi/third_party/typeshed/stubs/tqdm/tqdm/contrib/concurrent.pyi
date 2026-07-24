from _typeshed import SupportsWrite
from collections.abc import Callable, Iterable, Mapping
from typing import Any, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import Unpack

from ..std import tqdm

__all__ = ["thread_map", "process_map"]

_R = TypeVar("_R")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_T3 = TypeVar("_T3")
_T4 = TypeVar("_T4")
_T5 = TypeVar("_T5")

@type_check_only
class _TqdmKwargs(TypedDict, total=False):
    # Concurrent-specific parameters
    tqdm_class: type[tqdm[object]]
    max_workers: int | None
    chunksize: int
    lock_name: str
    # Standard tqdm parameters
    desc: str | None
    total: float | None
    leave: bool | None
    file: SupportsWrite[str] | None
    ncols: int | None
    mininterval: float
    maxinterval: float
    miniters: float | None
    ascii: bool | str | None
    disable: bool | None
    unit: str
    unit_scale: bool | float
    dynamic_ncols: bool
    smoothing: float
    bar_format: str | None
    initial: float
    position: int | None
    postfix: Mapping[str, object] | str | None
    unit_divisor: float
    write_bytes: bool | None
    lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None
    nrows: int | None
    colour: str | None
    delay: float | None

@overload
def thread_map(fn: Callable[[_T1], _R], iter1: Iterable[_T1], **tqdm_kwargs: Unpack[_TqdmKwargs]) -> list[_R]: ...
@overload
def thread_map(
    fn: Callable[[_T1, _T2], _R], iter1: Iterable[_T1], iter2: Iterable[_T2], /, **tqdm_kwargs: Unpack[_TqdmKwargs]
) -> list[_R]: ...
@overload
def thread_map(
    fn: Callable[[_T1, _T2, _T3], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def thread_map(
    fn: Callable[[_T1, _T2, _T3, _T4], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    iter4: Iterable[_T4],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def thread_map(
    fn: Callable[[_T1, _T2, _T3, _T4, _T5], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    iter4: Iterable[_T4],
    iter5: Iterable[_T5],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def thread_map(
    fn: Callable[..., _R],
    iter1: Iterable[Any],
    iter2: Iterable[Any],
    iter3: Iterable[Any],
    iter4: Iterable[Any],
    iter5: Iterable[Any],
    iter6: Iterable[Any],
    *iterables: Iterable[Any],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def process_map(fn: Callable[[_T1], _R], iter1: Iterable[_T1], **tqdm_kwargs: Unpack[_TqdmKwargs]) -> list[_R]: ...
@overload
def process_map(
    fn: Callable[[_T1, _T2], _R], iter1: Iterable[_T1], iter2: Iterable[_T2], **tqdm_kwargs: Unpack[_TqdmKwargs]
) -> list[_R]: ...
@overload
def process_map(
    fn: Callable[[_T1, _T2, _T3], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def process_map(
    fn: Callable[[_T1, _T2, _T3, _T4], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    iter4: Iterable[_T4],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def process_map(
    fn: Callable[[_T1, _T2, _T3, _T4, _T5], _R],
    iter1: Iterable[_T1],
    iter2: Iterable[_T2],
    iter3: Iterable[_T3],
    iter4: Iterable[_T4],
    iter5: Iterable[_T5],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
@overload
def process_map(
    fn: Callable[..., _R],
    iter1: Iterable[Any],
    iter2: Iterable[Any],
    iter3: Iterable[Any],
    iter4: Iterable[Any],
    iter5: Iterable[Any],
    iter6: Iterable[Any],
    *iterables: Iterable[Any],
    **tqdm_kwargs: Unpack[_TqdmKwargs],
) -> list[_R]: ...
