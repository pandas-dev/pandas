from typing import IO, Final, Literal, LiteralString, Self, TypeAlias, overload, type_check_only
from typing_extensions import Protocol

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

from scipy.io._typing import FileLike
from scipy.sparse import csc_array, csc_matrix
from scipy.sparse._base import _spbase

__all__ = ["hb_read", "hb_write"]

_ValueType: TypeAlias = Literal["real", "complex", "pattern", "integer"]
_Structure: TypeAlias = Literal["symmetric", "unsymmetric", "hermitian", "skewsymmetric", "rectangular"]
_Storage: TypeAlias = Literal["assembled", "elemental"]

_Real: TypeAlias = npc.integer | np.float32 | np.float64

@type_check_only
class _HasWidthAndRepeat(Protocol):
    @property
    def width(self, /) -> int: ...
    @property
    def repeat(self, /) -> int | None: ...

###

class MalformedHeader(Exception): ...
class LineOverflow(Warning): ...

class HBInfo:
    title: Final[str]
    key: Final[str]
    total_nlines: Final[int]
    pointer_nlines: Final[int]
    indices_nlines: Final[int]
    values_nlines: Final[int]
    pointer_format: Final[int]
    indices_format: Final[int]
    values_format: Final[int]
    pointer_dtype: Final[int]
    indices_dtype: Final[int]
    values_dtype: Final[int]
    pointer_nbytes_full: Final[int]
    indices_nbytes_full: Final[int]
    values_nbytes_full: Final[int]
    nrows: Final[int]
    ncols: Final[int]
    nnon_zeros: Final[int]
    nelementals: Final[int]
    mxtype: HBMatrixType

    @classmethod
    def from_data(
        cls, m: _spbase, title: str = "Default title", key: str = "0", mxtype: HBMatrixType | None = None, fmt: None = None
    ) -> Self: ...
    @classmethod
    def from_file(cls, fid: IO[str]) -> Self: ...

    #
    def __init__(
        self,
        /,
        title: str,
        key: str,
        total_nlines: int,
        pointer_nlines: int,
        indices_nlines: int,
        values_nlines: int,
        mxtype: HBMatrixType,
        nrows: int,
        ncols: int,
        nnon_zeros: int,
        pointer_format_str: str,
        indices_format_str: str,
        values_format_str: str,
        right_hand_sides_nlines: int = 0,
        nelementals: int = 0,
    ) -> None: ...
    def dump(self, /) -> str: ...

class HBMatrixType:
    value_type: Final[_ValueType]
    structure: Final[_Structure]
    storage: Final[_Storage]

    @property
    def fortran_format(self, /) -> LiteralString: ...
    @classmethod
    def from_fortran(cls, fmt: str) -> Self: ...

    #
    def __init__(self, /, value_type: _ValueType, structure: _Structure, storage: _Storage = "assembled") -> None: ...

class HBFile:
    @property
    def title(self, /) -> str: ...
    @property
    def key(self, /) -> str: ...
    @property
    def type(self, /) -> _ValueType: ...
    @property
    def structure(self, /) -> _Structure: ...
    @property
    def storage(self, /) -> _Storage: ...

    #
    def __init__(self, /, file: IO[str], hb_info: HBMatrixType | None = None) -> None: ...
    def read_matrix(self, /) -> csc_array[_Real]: ...
    def write_matrix(self, /, m: _spbase) -> None: ...

def _nbytes_full(fmt: _HasWidthAndRepeat, nlines: int) -> int: ...
def _expect_int(value: opt.AnyInt, msg: str | None = None) -> int: ...
def _read_hb_data(content: IO[str], header: HBInfo) -> csc_array[_Real]: ...
def _write_data(m: _spbase, fid: IO[str], header: HBInfo) -> None: ...

#
@overload
def hb_read(path_or_open_file: FileLike[str], *, spmatrix: onp.ToTrue = True) -> csc_array[_Real]: ...
@overload
def hb_read(path_or_open_file: FileLike[str], *, spmatrix: onp.ToFalse) -> csc_matrix[_Real]: ...

#
def hb_write(path_or_open_file: FileLike[str], m: _spbase, hb_info: HBInfo | None = None) -> None: ...
