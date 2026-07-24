import contextlib
from _typeshed import Incomplete, SupportsRead, SupportsWrite
from collections.abc import Callable, Iterable, Iterator, Mapping, MutableMapping
from types import TracebackType
from typing import Any, ClassVar, Generic, Literal, NoReturn, TypeVar, overload
from typing_extensions import Self

from ._monitor import TMonitor
from .utils import Comparable

__all__ = [
    "tqdm",
    "trange",
    "TqdmTypeError",
    "TqdmKeyError",
    "TqdmWarning",
    "TqdmExperimentalWarning",
    "TqdmDeprecationWarning",
    "TqdmMonitorWarning",
]

class TqdmTypeError(TypeError): ...
class TqdmKeyError(KeyError): ...

class TqdmWarning(Warning):
    def __init__(self, msg, fp_write=None) -> None: ...

class TqdmExperimentalWarning(TqdmWarning, FutureWarning): ...
class TqdmDeprecationWarning(TqdmWarning, DeprecationWarning): ...
class TqdmMonitorWarning(TqdmWarning, RuntimeWarning): ...

_T = TypeVar("_T")
_U = TypeVar("_U")

class tqdm(Comparable, Generic[_T]):
    monitor_interval: ClassVar[int]
    monitor: ClassVar[TMonitor | None]

    @staticmethod
    def format_sizeof(num: float, suffix: str = "", divisor: float = 1000) -> str: ...
    @staticmethod
    def format_interval(t: float) -> str: ...
    @staticmethod
    def format_num(n: float) -> str: ...
    @staticmethod
    def status_printer(file: SupportsWrite[str]) -> Callable[[str], None]: ...
    @staticmethod
    def format_meter(
        n: float,
        total: float,
        elapsed: float,
        ncols: int | None = None,
        prefix: str | None = "",
        ascii: bool | str | None = False,
        unit: str | None = "it",
        unit_scale: bool | float | None = False,
        rate: float | None = None,
        bar_format: str | None = None,
        postfix: str | Mapping[str, object] | None = None,
        unit_divisor: float | None = 1000,
        initial: float | None = 0,
        colour: str | None = None,
    ) -> str: ...
    @overload
    def __init__(
        self,
        iterable: Iterable[_T],
        desc: str | None = None,
        total: float | None = None,
        leave: bool | None = True,
        file: SupportsWrite[str] | None = None,
        ncols: int | None = None,
        mininterval: float = 0.1,
        maxinterval: float = 10.0,
        miniters: float | None = None,
        ascii: bool | str | None = None,
        disable: bool | None = False,
        unit: str = "it",
        unit_scale: bool | float = False,
        dynamic_ncols: bool = False,
        smoothing: float = 0.3,
        bar_format: str | None = None,
        initial: float = 0,
        position: int | None = None,
        postfix: Mapping[str, object] | str | None = None,
        unit_divisor: float = 1000,
        write_bytes: bool = False,
        lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None = None,
        nrows: int | None = None,
        colour: str | None = None,
        delay: float | None = 0,
        gui: bool = False,
        **kwargs,
    ) -> None: ...
    @overload
    def __init__(
        self: tqdm[NoReturn],
        iterable: None = None,
        desc: str | None = None,
        total: float | None = None,
        leave: bool | None = True,
        file: SupportsWrite[str] | None = None,
        ncols: int | None = None,
        mininterval: float = 0.1,
        maxinterval: float = 10.0,
        miniters: float | None = None,
        ascii: bool | str | None = None,
        disable: bool | None = False,
        unit: str = "it",
        unit_scale: bool | float = False,
        dynamic_ncols: bool = False,
        smoothing: float = 0.3,
        bar_format: str | None = None,
        initial: float = 0,
        position: int | None = None,
        postfix: Mapping[str, object] | str | None = None,
        unit_divisor: float = 1000,
        write_bytes: bool | None = False,
        lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None = None,
        nrows: int | None = None,
        colour: str | None = None,
        delay: float | None = 0,
        gui: bool = False,
        **kwargs,
    ) -> None: ...
    def __new__(cls, *_, **__) -> Self: ...
    @classmethod
    def write(cls, s: str, file: SupportsWrite[str] | None = None, end: str = "\n", nolock: bool = False) -> None: ...
    @classmethod
    def external_write_mode(
        cls, file: SupportsWrite[str] | None = None, nolock: bool = False
    ) -> contextlib._GeneratorContextManager[None]: ...
    @classmethod
    def set_lock(cls, lock) -> None: ...
    @classmethod
    def get_lock(cls): ...
    @classmethod
    def pandas(
        cls,
        *,
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
    ) -> None: ...

    iterable: Incomplete
    disable: Incomplete
    pos: Incomplete
    n: Incomplete
    total: Incomplete
    leave: Incomplete
    desc: Incomplete
    fp: Incomplete
    ncols: Incomplete
    nrows: Incomplete
    mininterval: Incomplete
    maxinterval: Incomplete
    miniters: Incomplete
    dynamic_miniters: Incomplete
    ascii: Incomplete
    unit: Incomplete
    unit_scale: Incomplete
    unit_divisor: Incomplete
    initial: Incomplete
    lock_args: Incomplete
    delay: Incomplete
    gui: Incomplete
    dynamic_ncols: Incomplete
    smoothing: Incomplete
    bar_format: Incomplete
    postfix: Incomplete
    colour: Incomplete
    last_print_n: Incomplete
    sp: Incomplete
    last_print_t: Incomplete
    start_t: Incomplete

    def __bool__(self) -> bool: ...
    def __len__(self) -> int: ...
    def __reversed__(self) -> Iterator[_T]: ...
    def __contains__(self, item: object) -> bool: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
    def __del__(self) -> None: ...
    def __hash__(self) -> int: ...
    def __iter__(self) -> Iterator[_T]: ...
    def update(self, n: float | None = 1) -> bool | None: ...
    def close(self) -> None: ...
    def clear(self, nolock: bool = False) -> None: ...
    def refresh(
        self, nolock: bool = False, lock_args: tuple[bool | None, float | None] | tuple[bool | None] | None = None
    ) -> None: ...
    def unpause(self) -> None: ...
    def reset(self, total: float | None = None) -> None: ...
    def set_description(self, desc: str | None = None, refresh: bool | None = True) -> None: ...
    def set_description_str(self, desc: str | None = None, refresh: bool | None = True) -> None: ...
    def set_postfix(self, ordered_dict: Mapping[str, object] | None = None, refresh: bool | None = True, **kwargs) -> None: ...
    def set_postfix_str(self, s: str = "", refresh: bool = True) -> None: ...
    def moveto(self, n) -> None: ...
    @property
    def format_dict(self) -> MutableMapping[str, Any]: ...
    def display(self, msg: str | None = None, pos: int | None = None) -> None: ...
    @overload
    @classmethod
    def wrapattr(
        cls, stream: SupportsRead[_U], method: Literal["read"], total: float | None = None, bytes: bool = True, **tqdm_kwargs
    ) -> contextlib._GeneratorContextManager[SupportsRead[_U]]: ...
    @overload
    @classmethod
    def wrapattr(
        cls, stream: SupportsWrite[_U], method: Literal["write"], total: float | None = None, bytes: bool = True, **tqdm_kwargs
    ) -> contextlib._GeneratorContextManager[SupportsWrite[_U]]: ...

@overload
def trange(
    start: int,
    stop: int,
    step: int | None = ...,
    *,
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
) -> tqdm[int]: ...
@overload
def trange(
    stop: int,
    *,
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
) -> tqdm[int]: ...
