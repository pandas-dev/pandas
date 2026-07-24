from _typeshed import (
    Incomplete,
    OpenBinaryModeReading,
    OpenBinaryModeWriting,
    OpenTextModeReading,
    OpenTextModeWriting,
    SupportsWrite,
    Unused,
)
from re import Pattern
from typing import IO, Any, ClassVar, Final, Generic, Literal, TextIO, TypeVar
from typing_extensions import deprecated

from docutils import TransformSpec, nodes

__docformat__: Final = "reStructuredText"

class InputError(OSError): ...
class OutputError(OSError): ...

def check_encoding(stream: TextIO, encoding: str) -> bool | None: ...
def error_string(err: BaseException) -> str: ...

_S = TypeVar("_S")

class Input(TransformSpec, Generic[_S]):
    component_type: ClassVar[str]
    default_source_path: ClassVar[str | None]
    encoding: str | None
    error_handler: str
    source: _S | None
    source_path: str | None
    successful_encoding: str | None = None
    def __init__(
        self,
        source: _S | None = None,
        source_path: str | None = None,
        encoding: str | None = "utf-8",
        error_handler: str = "strict",
    ) -> None: ...
    def read(self) -> str: ...
    def decode(self, data: str | bytes | bytearray) -> str: ...
    coding_slug: ClassVar[Pattern[bytes]]
    byte_order_marks: ClassVar[tuple[tuple[bytes, str], ...]]
    @deprecated("Deprecated and will be removed in Docutils 1.0.")
    def determine_encoding_from_data(self, data: str | bytes | bytearray) -> str | None: ...
    def isatty(self) -> bool: ...

class Output(TransformSpec):
    component_type: ClassVar[str]
    default_destination_path: ClassVar[str | None]
    encoding: Incomplete
    error_handler: Incomplete
    destination: Incomplete
    destination_path: Incomplete
    def __init__(
        self, destination=None, destination_path=None, encoding: str | None = None, error_handler: str = "strict"
    ) -> None: ...
    def write(self, data: str) -> Any: ...  # returns bytes or str
    def encode(self, data: str) -> Any: ...  # returns bytes or str

class ErrorOutput:
    destination: Incomplete
    encoding: Incomplete
    encoding_errors: Incomplete
    decoding_errors: Incomplete
    def __init__(
        self,
        destination: str | SupportsWrite[str] | SupportsWrite[bytes] | Literal[False] | None = None,
        encoding: str | None = None,
        encoding_errors: str = "backslashreplace",
        decoding_errors: str = "replace",
    ) -> None: ...
    def write(self, data: str | bytes | Exception) -> None: ...
    def close(self) -> None: ...
    def isatty(self) -> bool: ...

class FileInput(Input[IO[str]]):
    autoclose: bool
    def __init__(
        self,
        source=None,
        source_path=None,
        encoding: str | None = "utf-8",
        error_handler: str = "strict",
        autoclose: bool = True,
        mode: OpenTextModeReading | OpenBinaryModeReading = "r",
    ) -> None: ...
    def read(self) -> str: ...
    def readlines(self) -> list[str]: ...
    def close(self) -> None: ...

class FileOutput(Output):
    default_destination_path: ClassVar[str]
    mode: ClassVar[OpenTextModeWriting | OpenBinaryModeWriting]
    opened: bool
    autoclose: Incomplete
    destination: Incomplete
    destination_path: Incomplete
    def __init__(
        self,
        destination=None,
        destination_path=None,
        encoding=None,
        error_handler: str = "strict",
        autoclose: bool = True,
        handle_io_errors=None,
        mode=None,
    ) -> None: ...
    def open(self) -> None: ...
    def write(self, data): ...
    def close(self) -> None: ...

@deprecated("The `BinaryFileOutput` is deprecated by `FileOutput` and will be removed in Docutils 0.24.")
class BinaryFileOutput(FileOutput): ...

class StringInput(Input[str]):
    default_source_path: ClassVar[str]
    def read(self): ...

class StringOutput(Output):
    default_destination_path: ClassVar[str]
    destination: str | bytes  # only defined after call to write()
    def write(self, data): ...

class NullInput(Input[Any]):
    default_source_path: ClassVar[str]
    def read(self) -> str: ...

class NullOutput(Output):
    default_destination_path: ClassVar[str]
    def write(self, data: Unused) -> None: ...

class DocTreeInput(Input[nodes.document]):
    default_source_path: ClassVar[str]
    def read(self): ...
