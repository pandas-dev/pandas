from typing import Any, Final, Literal, Protocol, TypedDict, final, type_check_only

import optype as op

from ._data import _data_matrix

__all__ = ["load_npz", "save_npz"]

@final
@type_check_only
class _CanReadAndSeekBytes(Protocol):
    def read(self, length: int = ..., /) -> bytes: ...
    def seek(self, offset: int, whence: int, /) -> object: ...

@final
@type_check_only
class _PickleKwargs(TypedDict):
    allow_pickle: Literal[False]

###

PICKLE_KWARGS: Final[_PickleKwargs] = ...

def load_npz(file: op.io.ToPath | _CanReadAndSeekBytes) -> _data_matrix | Any: ...
def save_npz(file: op.io.ToPath[str] | op.io.CanWrite[bytes], matrix: _data_matrix, compressed: op.CanBool = True) -> None: ...
