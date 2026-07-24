from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable, Mapping
from typing import NoReturn, TypeVar, overload

from ..auto import tqdm as tqdm_auto
from .utils_worker import MonoWorker

__all__ = ["SlackIO", "tqdm_slack", "tsrange", "tqdm", "trange"]

class SlackIO(MonoWorker):
    client: Incomplete
    text: Incomplete
    message: Incomplete
    def __init__(self, token, channel) -> None: ...
    def write(self, s): ...

_T = TypeVar("_T")

class tqdm_slack(tqdm_auto[_T]):
    sio: Incomplete
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
        channel: int = ...,
        **kwargs,
    ) -> None: ...
    @overload
    def __init__(
        self: tqdm_slack[NoReturn],
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
        channel: int = ...,
        **kwargs,
    ) -> None: ...
    def display(  # type: ignore[override]
        self, *, msg: str | None = ..., pos: int | None = ..., close: bool = ..., bar_style=..., check_delay: bool = ...
    ) -> None: ...
    def clear(self, *args, **kwargs) -> None: ...

def tsrange(*args, **kwargs) -> tqdm_slack[int]: ...

tqdm = tqdm_slack
trange = tsrange
