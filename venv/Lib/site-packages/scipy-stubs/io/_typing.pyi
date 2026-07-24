# NOTE(scipy-stubs): This private module should not be used outside of scipy-stubs

from typing import IO, Literal

import optype as op

__all__ = "ByteOrder", "FileLike", "FileName"

type ByteOrder = Literal["S", "<", "little", ">", "big", "=", "native", "|", "I"]

type FileName = str | op.io.CanFSPath[str]
type FileLike[ByteSOrStr: (bytes, str)] = FileName | IO[ByteSOrStr]
