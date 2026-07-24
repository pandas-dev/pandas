from _typeshed import (
    FileDescriptorOrPath,
    OpenBinaryMode,
    OpenBinaryModeReading,
    OpenBinaryModeUpdating,
    OpenBinaryModeWriting,
    OpenTextMode,
)
from asyncio import AbstractEventLoop
from collections.abc import Callable
from concurrent.futures import Executor
from functools import _SingleDispatchCallable
from typing import Any, Literal, overload
from typing_extensions import TypeAlias

from ..base import AiofilesContextManager
from .binary import AsyncBufferedIOBase, AsyncBufferedReader, AsyncFileIO, AsyncIndirectBufferedIOBase, _UnknownAsyncBinaryIO
from .text import AsyncTextIndirectIOWrapper, AsyncTextIOWrapper

_Opener: TypeAlias = Callable[[str, int], int]

# Text mode: always returns AsyncTextIOWrapper
@overload
def open(
    file: FileDescriptorOrPath,
    mode: OpenTextMode = "r",
    buffering: int = -1,
    encoding: str | None = None,
    errors: str | None = None,
    newline: str | None = None,
    closefd: bool = True,
    opener: _Opener | None = None,
    *,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncTextIOWrapper]: ...

# Unbuffered binary: returns a FileIO
@overload
def open(
    file: FileDescriptorOrPath,
    mode: OpenBinaryMode,
    buffering: Literal[0],
    encoding: None = None,
    errors: None = None,
    newline: None = None,
    closefd: bool = True,
    opener: _Opener | None = None,
    *,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncFileIO]: ...

# Buffered binary reading/updating: AsyncBufferedReader
@overload
def open(
    file: FileDescriptorOrPath,
    mode: OpenBinaryModeReading | OpenBinaryModeUpdating,
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    errors: None = None,
    newline: None = None,
    closefd: bool = True,
    opener: _Opener | None = None,
    *,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedReader]: ...

# Buffered binary writing: AsyncBufferedIOBase
@overload
def open(
    file: FileDescriptorOrPath,
    mode: OpenBinaryModeWriting,
    buffering: Literal[-1, 1] = -1,
    encoding: None = None,
    errors: None = None,
    newline: None = None,
    closefd: bool = True,
    opener: _Opener | None = None,
    *,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[AsyncBufferedIOBase]: ...

# Buffering cannot be determined: fall back to _UnknownAsyncBinaryIO
@overload
def open(
    file: FileDescriptorOrPath,
    mode: OpenBinaryMode,
    buffering: int = -1,
    encoding: None = None,
    errors: None = None,
    newline: None = None,
    closefd: bool = True,
    opener: _Opener | None = None,
    *,
    loop: AbstractEventLoop | None = None,
    executor: Executor | None = None,
) -> AiofilesContextManager[_UnknownAsyncBinaryIO]: ...

wrap: _SingleDispatchCallable[Any]

stdin: AsyncTextIndirectIOWrapper
stdout: AsyncTextIndirectIOWrapper
stderr: AsyncTextIndirectIOWrapper
stdin_bytes: AsyncIndirectBufferedIOBase
stdout_bytes: AsyncIndirectBufferedIOBase
stderr_bytes: AsyncIndirectBufferedIOBase

__all__ = ("open", "stdin", "stdout", "stderr", "stdin_bytes", "stdout_bytes", "stderr_bytes")
