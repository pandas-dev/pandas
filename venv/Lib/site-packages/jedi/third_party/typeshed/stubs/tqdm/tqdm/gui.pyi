from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable, Mapping
from typing import NoReturn, TypeVar, overload

from .std import tqdm as std_tqdm

__all__ = ["tqdm_gui", "tgrange", "tqdm", "trange"]

_T = TypeVar("_T")

class tqdm_gui(std_tqdm[_T]):
    mpl: Incomplete
    plt: Incomplete
    toolbar: Incomplete
    mininterval: Incomplete
    xdata: Incomplete
    ydata: Incomplete
    zdata: Incomplete
    hspan: Incomplete
    wasion: Incomplete
    ax: Incomplete
    disable: bool
    def close(self) -> None: ...
    def clear(self, *_, **__) -> None: ...
    def display(self, *_, **__) -> None: ...
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
        **kwargs,
    ) -> None: ...
    @overload
    def __init__(
        self: tqdm_gui[NoReturn],
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
        **kwargs,
    ) -> None: ...

def tgrange(*args, **kwargs) -> tqdm_gui[int]: ...

tqdm = tqdm_gui
trange = tgrange
