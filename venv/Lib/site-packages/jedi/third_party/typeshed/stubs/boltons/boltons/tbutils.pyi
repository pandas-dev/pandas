from collections.abc import Iterable, Iterator, Mapping
from types import FrameType, TracebackType
from typing import Any, Generic, Literal, TypeVar
from typing_extensions import Self

class Callpoint:
    __slots__ = ("func_name", "lineno", "module_name", "module_path", "lasti", "line")
    func_name: str
    lineno: int
    module_name: str
    module_path: str
    lasti: int
    line: str
    def __init__(
        self, module_name: str, module_path: str, func_name: str, lineno: int, lasti: int, line: str | None = None
    ) -> None: ...
    def to_dict(self) -> dict[str, Any]: ...
    @classmethod
    def from_current(cls, level: int = 1) -> Self: ...
    @classmethod
    def from_frame(cls, frame: FrameType) -> Self: ...
    @classmethod
    def from_tb(cls, tb: TracebackType) -> Self: ...
    def tb_frame_str(self) -> str: ...

_CallpointT_co = TypeVar("_CallpointT_co", bound=Callpoint, covariant=True, default=Callpoint)

class TracebackInfo(Generic[_CallpointT_co]):
    callpoint_type: type[_CallpointT_co]
    frames: list[_CallpointT_co]
    def __init__(self, frames: list[_CallpointT_co]) -> None: ...
    @classmethod
    def from_frame(cls, frame: FrameType | None = None, level: int = 1, limit: int | None = None) -> Self: ...
    @classmethod
    def from_traceback(cls, tb: TracebackType | None = None, limit: int | None = None) -> Self: ...
    @classmethod
    def from_dict(cls, d: Mapping[Literal["frames"], list[_CallpointT_co]]) -> Self: ...
    def to_dict(self) -> dict[str, list[dict[str, _CallpointT_co]]]: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[_CallpointT_co]: ...
    def get_formatted(self) -> str: ...

_TracebackInfoT_co = TypeVar("_TracebackInfoT_co", bound=TracebackInfo, covariant=True, default=TracebackInfo)

class ExceptionInfo(Generic[_TracebackInfoT_co]):
    tb_info_type: type[_TracebackInfoT_co]
    exc_type: str
    exc_msg: str
    tb_info: _TracebackInfoT_co
    def __init__(self, exc_type: str, exc_msg: str, tb_info: _TracebackInfoT_co) -> None: ...
    @classmethod
    def from_exc_info(cls, exc_type: type[BaseException], exc_value: BaseException, traceback: TracebackType) -> Self: ...
    @classmethod
    def from_current(cls) -> Self: ...
    def to_dict(self) -> dict[str, str | dict[str, list[FrameType]]]: ...
    def get_formatted(self) -> str: ...
    def get_formatted_exception_only(self) -> str: ...

class ContextualCallpoint(Callpoint):
    local_reprs: dict[Any, Any]
    pre_lines: list[str]
    post_lines: list[str]
    def __init__(self, *a, **kw) -> None: ...
    @classmethod
    def from_frame(cls, frame: FrameType) -> Self: ...
    @classmethod
    def from_tb(cls, tb: TracebackType) -> Self: ...
    def to_dict(self) -> dict[str, Any]: ...

class ContextualTracebackInfo(TracebackInfo[ContextualCallpoint]):
    callpoint_type: type[ContextualCallpoint]

class ContextualExceptionInfo(ExceptionInfo[ContextualTracebackInfo]):
    tb_info_type: type[ContextualTracebackInfo]

def print_exception(
    etype: type[BaseException] | None,
    value: BaseException | None,
    tb: TracebackType | None,
    limit: int | None = None,
    file: str | None = None,
) -> None: ...

class ParsedException:
    exc_type: str
    exc_msg: str
    frames: list[FrameType]
    def __init__(self, exc_type_name: str, exc_msg: str, frames: Iterable[Mapping[str, Any]] | None = None) -> None: ...
    @property
    def source_file(self) -> str | None: ...
    def to_dict(self) -> dict[str, str | list[FrameType]]: ...
    def to_string(self) -> str: ...
    @classmethod
    def from_string(cls, tb_str: str) -> Self: ...

ParsedTB = ParsedException

__all__ = [
    "ExceptionInfo",
    "TracebackInfo",
    "Callpoint",
    "ContextualExceptionInfo",
    "ContextualTracebackInfo",
    "ContextualCallpoint",
    "print_exception",
    "ParsedException",
]
