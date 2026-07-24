import sys
from typing import Any
from typing_extensions import TypeAlias

from gevent._fileobjectcommon import FileObjectBlock as FileObjectBlock, FileObjectThread as FileObjectThread

if sys.platform != "win32":
    import io
    from _typeshed import (
        FileDescriptorOrPath,
        OpenBinaryMode,
        OpenBinaryModeReading,
        OpenBinaryModeUpdating,
        OpenBinaryModeWriting,
        OpenTextMode,
    )
    from typing import IO, AnyStr, Literal, overload

    from gevent._fileobjectcommon import _IOT, FileObjectBase

    # this is implemented in _fileobjectposix and technically uses an undocumented subclass
    # of RawIOBase, but the interface is the same, so it doesn't seem worth it to add
    # annotations for it. _fileobjectcommon was barely worth it due to the common base class
    # of all three FileObject types
    class FileObjectPosix(FileObjectBase[_IOT, AnyStr]):
        default_bufsize = io.DEFAULT_BUFFER_SIZE
        fileio: io.RawIOBase
        # Text mode: always binds a TextIOWrapper
        @overload
        def __init__(
            self: FileObjectPosix[io.TextIOWrapper, str],
            fobj: FileDescriptorOrPath,
            mode: OpenTextMode = "r",
            bufsize: int | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: int | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...

        # Unbuffered binary mode: binds a FileIO
        @overload
        def __init__(
            self: FileObjectPosix[io.FileIO, bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryMode,
            bufsize: Literal[0],
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: Literal[0] | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...
        @overload
        def __init__(
            self: FileObjectPosix[io.FileIO, bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryMode,
            bufsize: Literal[0] | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            *,
            buffering: Literal[0],
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...

        # Buffering is on: return BufferedRandom, BufferedReader, or BufferedWriter
        @overload
        def __init__(
            self: FileObjectPosix[io.BufferedRandom, bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryModeUpdating,
            bufsize: Literal[-1, 1] | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: Literal[-1, 1] | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...
        @overload
        def __init__(
            self: FileObjectPosix[io.BufferedWriter, bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryModeWriting,
            bufsize: Literal[-1, 1] | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: Literal[-1, 1] | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...
        @overload
        def __init__(
            self: FileObjectPosix[io.BufferedReader, bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryModeReading,
            bufsize: Literal[-1, 1] | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: Literal[-1, 1] | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...

        # Buffering cannot be determined: fall back to BinaryIO
        @overload
        def __init__(
            self: FileObjectPosix[IO[bytes], bytes],
            fobj: FileDescriptorOrPath,
            mode: OpenBinaryMode,
            bufsize: int | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: int | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...

        # Fallback if mode is not specified
        @overload
        def __init__(
            self: FileObjectPosix[IO[Any], Any],
            fobj: FileDescriptorOrPath,
            mode: str,
            bufsize: int | None = None,
            close: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
            buffering: int | None = None,
            closefd: bool | None = None,
            atomic_write: bool = False,
        ) -> None: ...

    _FileObjectType: TypeAlias = type[FileObjectPosix[Any, Any] | FileObjectBlock[Any, Any] | FileObjectThread[Any, Any]]
    __all__ = ["FileObjectPosix", "FileObjectThread", "FileObjectBlock", "FileObject"]
else:
    _FileObjectType: TypeAlias = type[FileObjectBlock[Any, Any] | FileObjectThread[Any, Any]]
    __all__ = ["FileObjectThread", "FileObjectBlock", "FileObject"]

FileObject: _FileObjectType
