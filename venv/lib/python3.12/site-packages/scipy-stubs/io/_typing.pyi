# NOTE(scipy-stubs): This private module should not be used outside of scipy-stubs

from typing import IO, Literal, TypeAlias, TypeVar

import optype as op

__all__ = "ByteOrder", "FileLike", "FileName"

ByteOrder: TypeAlias = Literal["S", "<", "little", ">", "big", "=", "native", "|", "I"]

_ByteSOrStr = TypeVar("_ByteSOrStr", bytes, str)
FileName: TypeAlias = str | op.io.CanFSPath[str]
FileLike: TypeAlias = FileName | IO[_ByteSOrStr]
