import threading
from _typeshed import SupportsWrite, TraceFunction
from collections.abc import Callable
from pdb import Pdb
from types import FrameType
from typing import Any, Literal, Protocol, TypeVar, type_check_only
from typing_extensions import TypeAlias

from .watch_element import WatchElement

_T = TypeVar("_T")

# Alias used for fields that must always be valid identifiers
# A string `x` counts as a valid identifier if both the following are True
# (1) `x.isidentifier()` evaluates to `True`
# (2) `keyword.iskeyword(x)` evaluates to `False`
_Identifier: TypeAlias = str

class Watch:
    # User-defined callbacks passed to `__call__()` or `config()` set as instance variables have arguments of type `Any` to be
    # compatible with more precisely-annotated signatures.

    custom_printer: Callable[[Any], None] | None
    enable: bool
    file: str | SupportsWrite[str] | None
    pdb: Pdb | None
    pdb_enable: bool
    set_lock: threading.Lock
    stack_limit: int | None
    tracefunc_lock: threading.Lock
    tracefunc_stack: list[TraceFunction | None]
    watch_list: list[WatchElement]

    def __init__(self) -> None: ...
    def __call__(
        self,
        *args: object,
        alias: str = ...,
        callback: Callable[[FrameType, WatchElement, tuple[str, str, int | None]], None] = ...,
        cmp: Callable[[Any, Any], bool] = ...,  # User-defined comparison callback; compares 2 arguments of any type
        copy: Callable[[_T], _T] = ...,
        # User-defined printing callback; writes a string representation of any object to a stream
        custom_printer: Callable[[Any], None] = ...,
        deepcopy: bool = False,
        file: str | SupportsWrite[str] = ...,
        stack_limit: int | None = 5,
        track: Literal["object", "variable"] = ...,
        when: Callable[[Any], bool] = ...,  # User-defined callback for conditional watchpoints
    ) -> None: ...
    def config(
        self,
        *,
        callback: Callable[[FrameType, WatchElement, tuple[str, str, int | None]], None] = ...,
        pdb: Literal[True] = True,
        file: str | SupportsWrite[str] = ...,
        stack_limit: int | None = 5,
        custom_printer: Callable[[Any], None] = ...,  # User-defined printing callback
    ) -> None: ...
    def install(self, func: _Identifier = "watch") -> None: ...
    def restore(self) -> None: ...
    def start_trace(self, frame: FrameType) -> None: ...
    def stop_trace(self, frame: FrameType) -> None: ...
    def tracefunc(self, frame: FrameType, event: str, arg: object) -> _TraceFunc: ...
    def uninstall(self, func: _Identifier = "watch") -> None: ...
    def unwatch(self, *args: object) -> None: ...

@type_check_only
class _TraceFunc(Protocol):
    def __call__(self, frame: FrameType, event: str, arg: object) -> _TraceFunc: ...
