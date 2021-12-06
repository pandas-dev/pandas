import _compression
import sys
import zlib
from _typeshed import AnyPath, ReadableBuffer
from io import FileIO
from typing import Any, Optional, Protocol, TextIO, Union, overload
from typing_extensions import Literal

_ReadBinaryMode = Literal["r", "rb"]
_WriteBinaryMode = Literal["a", "ab", "w", "wb", "x", "xb"]
_OpenTextMode = Literal["rt", "at", "wt", "xt"]

READ: Literal[1]
WRITE: Literal[2]

class _ReadableFileobj(Protocol):
    def read(self, __n: int) -> bytes: ...
    def seek(self, __n: int) -> Any: ...
    # The following attributes and methods are optional:
    # name: str
    # mode: str
    # def fileno() -> int: ...

class _WritableFileobj(Protocol):
    def write(self, __b: bytes) -> Any: ...
    def flush(self) -> Any: ...
    # The following attributes and methods are optional:
    # name: str
    # mode: str
    # def fileno() -> int: ...

@overload
def open(
    filename: Union[AnyPath, _ReadableFileobj],
    mode: _ReadBinaryMode = ...,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> GzipFile: ...
@overload
def open(
    filename: Union[AnyPath, _WritableFileobj],
    mode: _WriteBinaryMode,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> GzipFile: ...
@overload
def open(
    filename: AnyPath,
    mode: _OpenTextMode,
    compresslevel: int = ...,
    encoding: Optional[str] = ...,
    errors: Optional[str] = ...,
    newline: Optional[str] = ...,
) -> TextIO: ...
@overload
def open(
    filename: Union[AnyPath, _ReadableFileobj, _WritableFileobj],
    mode: str,
    compresslevel: int = ...,
    encoding: Optional[str] = ...,
    errors: Optional[str] = ...,
    newline: Optional[str] = ...,
) -> Union[GzipFile, TextIO]: ...

class _PaddedFile:
    file: _ReadableFileobj
    def __init__(self, f: _ReadableFileobj, prepend: bytes = ...) -> None: ...
    def read(self, size: int) -> bytes: ...
    def prepend(self, prepend: bytes = ...) -> None: ...
    def seek(self, off: int) -> int: ...
    def seekable(self) -> bool: ...

if sys.version_info >= (3, 8):
    class BadGzipFile(OSError): ...

class GzipFile(_compression.BaseStream):
    myfileobj: Optional[FileIO]
    mode: Literal[1, 2]
    name: str
    compress: zlib._Compress
    fileobj: Union[_ReadableFileobj, _WritableFileobj]
    @overload
    def __init__(
        self,
        filename: Optional[AnyPath],
        mode: _ReadBinaryMode,
        compresslevel: int = ...,
        fileobj: Optional[_ReadableFileobj] = ...,
        mtime: Optional[float] = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        mode: _ReadBinaryMode,
        compresslevel: int = ...,
        fileobj: Optional[_ReadableFileobj] = ...,
        mtime: Optional[float] = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        filename: Optional[AnyPath],
        mode: _WriteBinaryMode,
        compresslevel: int = ...,
        fileobj: Optional[_WritableFileobj] = ...,
        mtime: Optional[float] = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        mode: _WriteBinaryMode,
        compresslevel: int = ...,
        fileobj: Optional[_WritableFileobj] = ...,
        mtime: Optional[float] = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        filename: Optional[AnyPath] = ...,
        mode: Optional[str] = ...,
        compresslevel: int = ...,
        fileobj: Union[_ReadableFileobj, _WritableFileobj, None] = ...,
        mtime: Optional[float] = ...,
    ) -> None: ...
    @property
    def filename(self) -> str: ...
    @property
    def mtime(self) -> Optional[int]: ...
    crc: int
    def write(self, data: ReadableBuffer) -> int: ...
    def read(self, size: Optional[int] = ...) -> bytes: ...
    def read1(self, size: int = ...) -> bytes: ...
    def peek(self, n: int) -> bytes: ...
    @property
    def closed(self) -> bool: ...
    def close(self) -> None: ...
    def flush(self, zlib_mode: int = ...) -> None: ...
    def fileno(self) -> int: ...
    def rewind(self) -> None: ...
    def readable(self) -> bool: ...
    def writable(self) -> bool: ...
    def seekable(self) -> bool: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def readline(self, size: Optional[int] = ...) -> bytes: ...

class _GzipReader(_compression.DecompressReader):
    def __init__(self, fp: _ReadableFileobj) -> None: ...
    def read(self, size: int = ...) -> bytes: ...

if sys.version_info >= (3, 8):
    def compress(data: bytes, compresslevel: int = ..., *, mtime: Optional[float] = ...) -> bytes: ...

else:
    def compress(data: bytes, compresslevel: int = ...) -> bytes: ...

def decompress(data: bytes) -> bytes: ...
