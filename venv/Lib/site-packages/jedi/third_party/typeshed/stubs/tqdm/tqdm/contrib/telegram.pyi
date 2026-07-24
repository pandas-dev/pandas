from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable, Mapping
from typing import NoReturn, TypeVar, overload

from ..auto import tqdm as tqdm_auto
from .utils_worker import MonoWorker

__all__ = ["TelegramIO", "tqdm_telegram", "ttgrange", "tqdm", "trange"]

class TelegramIO(MonoWorker):
    API: str
    token: Incomplete
    chat_id: Incomplete
    session: Incomplete
    text: Incomplete
    def __init__(self, token, chat_id) -> None: ...
    @property
    def message_id(self): ...
    def write(self, s: str) -> Incomplete | None: ...
    def delete(self): ...

_T = TypeVar("_T")

class tqdm_telegram(tqdm_auto[_T]):
    tgio: Incomplete
    @overload
    def __init__(
        self,
        iterable: Iterable[_T],
        desc: str | None = ...,
        total: float | None = ...,
        leave: bool | None = ...,
        file: SupportsWrite[str] | None = ...,
        ncols: int | None = ...,
        mininterval: float = ...,
        maxinterval: float = ...,
        miniters: float | None = ...,
        ascii: bool | str | None = ...,
        disable: bool | None = ...,
        unit: str = ...,
        unit_scale: bool | float = ...,
        dynamic_ncols: bool = ...,
        smoothing: float = ...,
        bar_format: str | None = ...,
        initial: float = ...,
        position: int | None = ...,
        postfix: Mapping[str, object] | str | None = ...,
        unit_divisor: float = ...,
        write_bytes: bool | None = ...,
        lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None = ...,
        nrows: int | None = ...,
        colour: str | None = ...,
        delay: float | None = ...,
        gui: bool = ...,
        token: str = ...,
        chat_id: str = ...,
        **kwargs,
    ) -> None: ...
    @overload
    def __init__(
        self: tqdm_telegram[NoReturn],
        iterable: None = None,
        desc: str | None = ...,
        total: float | None = ...,
        leave: bool | None = ...,
        file: SupportsWrite[str] | None = ...,
        ncols: int | None = ...,
        mininterval: float = ...,
        maxinterval: float = ...,
        miniters: float | None = ...,
        ascii: bool | str | None = ...,
        disable: bool | None = ...,
        unit: str = ...,
        unit_scale: bool | float = ...,
        dynamic_ncols: bool = ...,
        smoothing: float = ...,
        bar_format: str | None = ...,
        initial: float = ...,
        position: int | None = ...,
        postfix: Mapping[str, object] | str | None = ...,
        unit_divisor: float = ...,
        write_bytes: bool | None = ...,
        lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None = ...,
        nrows: int | None = ...,
        colour: str | None = ...,
        delay: float | None = ...,
        gui: bool = ...,
        token: str = ...,
        chat_id: str = ...,
        **kwargs,
    ) -> None: ...
    def display(  # type: ignore[override]
        self, *, msg: str | None = ..., pos: int | None = ..., close: bool = ..., bar_style=..., check_delay: bool = ...
    ) -> None: ...
    def clear(self, *args, **kwargs) -> None: ...
    def close(self) -> None: ...

def ttgrange(*args, **kwargs) -> tqdm_telegram[int]: ...

tqdm = tqdm_telegram
trange = ttgrange
