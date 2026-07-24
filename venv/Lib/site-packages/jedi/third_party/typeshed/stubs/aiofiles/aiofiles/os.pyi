import sys
from _typeshed import BytesPath, FileDescriptorOrPath, GenericPath, ReadableBuffer, StrOrBytesPath, StrPath
from asyncio.events import AbstractEventLoop
from collections.abc import Sequence
from concurrent.futures import Executor
from os import _ScandirIterator, stat_result
from typing import AnyStr, overload

from aiofiles import ospath
from aiofiles.base import wrap as wrap

__all__ = [
    "path",
    "stat",
    "rename",
    "renames",
    "replace",
    "remove",
    "unlink",
    "mkdir",
    "makedirs",
    "rmdir",
    "removedirs",
    "link",
    "symlink",
    "readlink",
    "listdir",
    "scandir",
    "access",
    "wrap",
    "getcwd",
]

if sys.platform != "win32":
    __all__ += ["statvfs", "sendfile"]

path = ospath

async def stat(
    path: FileDescriptorOrPath,
    *,
    dir_fd: int | None = None,
    follow_symlinks: bool = True,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> stat_result: ...
async def rename(
    src: StrOrBytesPath,
    dst: StrOrBytesPath,
    *,
    src_dir_fd: int | None = None,
    dst_dir_fd: int | None = None,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def renames(
    old: StrOrBytesPath, new: StrOrBytesPath, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> None: ...
async def replace(
    src: StrOrBytesPath,
    dst: StrOrBytesPath,
    *,
    src_dir_fd: int | None = None,
    dst_dir_fd: int | None = None,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def remove(
    path: StrOrBytesPath, *, dir_fd: int | None = None, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> None: ...
async def unlink(
    path: StrOrBytesPath, *, dir_fd: int | None = None, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> None: ...
async def mkdir(
    path: StrOrBytesPath,
    mode: int = 511,
    *,
    dir_fd: int | None = None,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def makedirs(
    name: StrOrBytesPath,
    mode: int = 511,
    exist_ok: bool = False,
    *,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def link(
    src: StrOrBytesPath,
    dst: StrOrBytesPath,
    *,
    src_dir_fd: int | None = None,
    dst_dir_fd: int | None = None,
    follow_symlinks: bool = True,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def symlink(
    src: StrOrBytesPath,
    dst: StrOrBytesPath,
    target_is_directory: bool = False,
    *,
    dir_fd: int | None = None,
    loop: AbstractEventLoop | None = ...,
    executor: Executor | None = ...,
) -> None: ...
async def readlink(
    path: AnyStr, *, dir_fd: int | None = None, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> AnyStr: ...
async def rmdir(
    path: StrOrBytesPath, *, dir_fd: int | None = None, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> None: ...
async def removedirs(name: StrOrBytesPath, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...) -> None: ...
@overload
async def scandir(
    path: None = None, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> _ScandirIterator[str]: ...
@overload
async def scandir(
    path: int, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> _ScandirIterator[str]: ...
@overload
async def scandir(
    path: GenericPath[AnyStr], *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> _ScandirIterator[AnyStr]: ...
@overload
async def listdir(
    path: StrPath | None, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...
) -> list[str]: ...
@overload
async def listdir(path: BytesPath, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...) -> list[bytes]: ...
@overload
async def listdir(path: int, *, loop: AbstractEventLoop | None = ..., executor: Executor | None = ...) -> list[str]: ...
async def access(
    path: FileDescriptorOrPath, mode: int, *, dir_fd: int | None = None, effective_ids: bool = False, follow_symlinks: bool = True
) -> bool: ...
async def getcwd() -> str: ...

if sys.platform != "win32":
    from os import statvfs_result

    @overload
    async def sendfile(
        out_fd: int,
        in_fd: int,
        offset: int | None,
        count: int,
        *,
        loop: AbstractEventLoop | None = ...,
        executor: Executor | None = ...,
    ) -> int: ...
    @overload
    async def sendfile(
        out_fd: int,
        in_fd: int,
        offset: int,
        count: int,
        headers: Sequence[ReadableBuffer] = (),
        trailers: Sequence[ReadableBuffer] = (),
        flags: int = 0,
        *,
        loop: AbstractEventLoop | None = ...,
        executor: Executor | None = ...,
    ) -> int: ...  # FreeBSD and Mac OS X only
    async def statvfs(path: FileDescriptorOrPath) -> statvfs_result: ...  # Unix only
