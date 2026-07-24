import argparse
from _typeshed import ExcInfo, FileDescriptorOrPath, Incomplete, OptExcInfo, SupportsWrite
from collections.abc import Mapping, Sequence
from typing import Any
from typing_extensions import Self, deprecated

@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
class _CannotLoadFile(Exception): ...

@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
class _Outputter:
    __slots__ = ("_formatter", "_stdout", "_stderr", "__weakref__")
    _formatter: _PlainFormatter | _PrettyFormatter
    _stdout: SupportsWrite[str]
    _stderr: SupportsWrite[str]

    def __init__(
        self, formatter: _PlainFormatter | _PrettyFormatter, stdout: SupportsWrite[str], stderr: SupportsWrite[str]
    ) -> None: ...
    @classmethod
    def from_arguments(
        cls, arguments: Mapping[str, Incomplete], stdout: SupportsWrite[str], stderr: SupportsWrite[str]
    ) -> Self: ...
    def load(self, path: FileDescriptorOrPath) -> Any: ...  # result of json.load()
    def filenotfound_error(self, *, path: FileDescriptorOrPath, exc_info: OptExcInfo) -> None: ...
    def parsing_error(self, *, path: FileDescriptorOrPath, exc_info: ExcInfo) -> None: ...
    def validation_error(self, *, instance_path: FileDescriptorOrPath, error: BaseException) -> None: ...
    def validation_success(self, *, instance_path: FileDescriptorOrPath) -> None: ...

@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
class _PrettyFormatter:
    __slots__ = ("__weakref__",)
    _ERROR_MSG: str
    _SUCCESS_MSG: str

    def __init__(self) -> None: ...
    def filenotfound_error(self, path: FileDescriptorOrPath, exc_info: OptExcInfo) -> str: ...
    def parsing_error(self, path: FileDescriptorOrPath, exc_info: ExcInfo) -> str: ...
    def validation_error(self, instance_path: FileDescriptorOrPath, error: BaseException) -> str: ...
    def validation_success(self, instance_path: FileDescriptorOrPath) -> str: ...

@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
class _PlainFormatter:
    __slots__ = ("_error_format", "__weakref__")
    _error_format: str

    def __init__(self, error_format: str) -> None: ...
    def filenotfound_error(self, path: FileDescriptorOrPath, exc_info: OptExcInfo) -> str: ...
    def parsing_error(self, path: FileDescriptorOrPath, exc_info: ExcInfo) -> str: ...
    def validation_error(self, instance_path: FileDescriptorOrPath, error: BaseException) -> str: ...
    def validation_success(self, instance_path: FileDescriptorOrPath) -> str: ...

parser: argparse.ArgumentParser

@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
def parse_args(args: Sequence[str] | None) -> dict[str, Any]: ...  # result of vars(argparse.Namespace())
@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
def main(args: Sequence[str] = ...) -> None: ...
@deprecated("The jsonschema CLI is deprecated and will be removed in a future version. Please use check-jsonschema instead.")
def run(
    arguments: Mapping[str, Incomplete],
    stdout: SupportsWrite[str] = ...,
    stderr: SupportsWrite[str] = ...,
    stdin: SupportsWrite[str] = ...,
) -> int: ...
