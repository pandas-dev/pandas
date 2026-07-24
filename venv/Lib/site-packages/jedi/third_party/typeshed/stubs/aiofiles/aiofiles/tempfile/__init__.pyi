import sys
from _typeshed import (
    BytesPath,
    OpenBinaryMode,
    OpenBinaryModeReading,
    OpenBinaryModeUpdating,
    OpenBinaryModeWriting,
    OpenTextMode,
    StrOrBytesPath,
    StrPath,
)
from asyncio import AbstractEventLoop
from concurrent.futures import Executor
from typing import AnyStr, Literal, overload

from ..base import AiofilesContextManager
from ..threadpool.binary import AsyncBufferedIOBase, AsyncBufferedReader, AsyncFileIO
from ..threadpool.text import AsyncTextIOWrapper

# Text mode: always returns AsyncTextIOWrapper
@overload
def TemporaryFile(
    mode: OpenTextMode,
    buffering: int = -1,
    encoding: str | None = None,
    newline: str | None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncTextIOWrapper]: ...

# Unbuffered binary: returns a FileIO
@overload
def TemporaryFile(
    mode: OpenBinaryMode,
    buffering: Literal[0],
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncFileIO]: ...

# Buffered binary reading/updating: AsyncBufferedReader
@overload
def TemporaryFile(
    mode: OpenBinaryModeReading | OpenBinaryModeUpdating = "w+b",
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedReader]: ...

# Buffered binary writing: AsyncBufferedIOBase
@overload
def TemporaryFile(
    mode: OpenBinaryModeWriting,
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedIOBase]: ...

# 3.12 added `delete_on_close`
if sys.version_info >= (3, 12):
    # Text mode: always returns AsyncTextIOWrapper
    @overload
    def NamedTemporaryFile(
        mode: OpenTextMode,
        buffering: int = -1,
        encoding: str | None = None,
        newline: str | None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        delete_on_close: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncTextIOWrapper]: ...

    # Unbuffered binary: returns a FileIO
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryMode,
        buffering: Literal[0],
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        delete_on_close: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncFileIO]: ...

    # Buffered binary reading/updating: AsyncBufferedReader
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryModeReading | OpenBinaryModeUpdating = "w+b",
        buffering: Literal[-1, 1] = -1,
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        delete_on_close: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncBufferedReader]: ...

    # Buffered binary writing: AsyncBufferedIOBase
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryModeWriting,
        buffering: Literal[-1, 1] = -1,
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        delete_on_close: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncBufferedIOBase]: ...

else:
    # Text mode: always returns AsyncTextIOWrapper
    @overload
    def NamedTemporaryFile(
        mode: OpenTextMode,
        buffering: int = -1,
        encoding: str | None = None,
        newline: str | None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncTextIOWrapper]: ...

    # Unbuffered binary: returns a FileIO
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryMode,
        buffering: Literal[0],
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncFileIO]: ...

    # Buffered binary reading/updating: AsyncBufferedReader
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryModeReading | OpenBinaryModeUpdating = "w+b",
        buffering: Literal[-1, 1] = -1,
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncBufferedReader]: ...

    # Buffered binary writing: AsyncBufferedIOBase
    @overload
    def NamedTemporaryFile(
        mode: OpenBinaryModeWriting,
        buffering: Literal[-1, 1] = -1,
        encoding: None = None,
        newline: None = None,
        suffix: AnyStr | None = None,
        prefix: AnyStr | None = None,
        dir: StrOrBytesPath | None = None,
        delete: bool = True,
        loop: AbstractEventLoop | None = None,
        executor: Executor | None = None,
    ) -> AiofilesContextManager[AsyncBufferedIOBase]: ...

# Text mode: always returns AsyncTextIOWrapper
@overload
def SpooledTemporaryFile(
    max_size: int = 0,
    *,
    mode: OpenTextMode,
    buffering: int = -1,
    encoding: str | None = None,
    newline: str | None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncTextIOWrapper]: ...
@overload
def SpooledTemporaryFile(
    max_size: int,
    mode: OpenTextMode,
    buffering: int = -1,
    encoding: str | None = None,
    newline: str | None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncTextIOWrapper]: ...

# Unbuffered binary: returns a FileIO
@overload
def SpooledTemporaryFile(
    max_size: int = 0,
    mode: OpenBinaryMode = "w+b",
    *,
    buffering: Literal[0],
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncFileIO]: ...
@overload
def SpooledTemporaryFile(
    max_size: int,
    mode: OpenBinaryMode,
    buffering: Literal[0],
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncFileIO]: ...

# Buffered binary reading/updating: AsyncBufferedReader
@overload
def SpooledTemporaryFile(
    max_size: int = 0,
    mode: OpenBinaryModeReading | OpenBinaryModeUpdating = "w+b",
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedReader]: ...

# Buffered binary writing: AsyncBufferedIOBase
@overload
def SpooledTemporaryFile(
    max_size: int = 0,
    *,
    mode: OpenBinaryModeWriting,
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedIOBase]: ...
@overload
def SpooledTemporaryFile(
    max_size: int,
    mode: OpenBinaryModeWriting,
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    newline: None = None,
    suffix: AnyStr | None = None,
    prefix: AnyStr | None = None,
    dir: StrOrBytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedIOBase]: ...
@overload
def TemporaryDirectory(
    suffix: str | None = None,
    prefix: str | None = None,
    dir: StrPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManagerTempDir: ...
@overload
def TemporaryDirectory(
    suffix: bytes | None = None,
    prefix: bytes | None = None,
    dir: BytesPath | None = None,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManagerTempDir: ...

class AiofilesContextManagerTempDir(AiofilesContextManager[str]):
    async def __aenter__(self) -> str: ...

__all__ = ["NamedTemporaryFile", "TemporaryFile", "SpooledTemporaryFile", "TemporaryDirectory"]
