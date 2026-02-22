from typing import ClassVar, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import Unpack

import optype.numpy as onp
import optype.numpy.compat as npc

from ._fast_matrix_market import mminfo, mmread, mmwrite
from ._typing import FileLike
from scipy.sparse import coo_array, coo_matrix, sparray, spmatrix

__all__ = ["MMFile", "mminfo", "mmread", "mmwrite"]

_Format: TypeAlias = Literal["coordinate", "array"]
_Field: TypeAlias = Literal["real", "complex", "pattern", "integer"]
_Symmetry: TypeAlias = Literal["general", "symmetric", "skew-symmetric", "hermitian"]
_Info: TypeAlias = tuple[int, int, int, _Format, _Field, _Symmetry]

@type_check_only
class _MMFileKwargs(TypedDict, total=False):
    rows: int
    cols: int
    entries: int
    format: _Format
    field: _Field
    symmetry: _Symmetry

###

class MMFile:
    __slots__ = "_cols", "_entries", "_field", "_format", "_rows", "_symmetry"

    FORMAT_COORDINATE: ClassVar[str] = "coordinate"
    FORMAT_ARRAY: ClassVar[str] = "array"
    FORMAT_VALUES: ClassVar[tuple[str, ...]] = "coordinate", "array"

    FIELD_INTEGER: ClassVar[str] = "integer"
    FIELD_UNSIGNED: ClassVar[str] = "unsigned-integer"
    FIELD_REAL: ClassVar[str] = "real"
    FIELD_COMPLEX: ClassVar[str] = "complex"
    FIELD_PATTERN: ClassVar[str] = "pattern"
    FIELD_VALUES: ClassVar[tuple[str, ...]] = "integer", "unsigned-integer", "real", "complex", "pattern"

    SYMMETRY_GENERAL: ClassVar[str] = "general"
    SYMMETRY_SYMMETRIC: ClassVar[str] = "symmetric"
    SYMMETRY_SKEW_SYMMETRIC: ClassVar[str] = "skew-symmetric"
    SYMMETRY_HERMITIAN: ClassVar[str] = "hermitian"
    SYMMETRY_VALUES: ClassVar[tuple[str, ...]] = "general", "symmetric", "skew-symmetric", "hermitian"

    DTYPES_BY_FIELD: ClassVar[dict[_Field, Literal["intp", "uint64", "d", "D"]]] = ...

    @property
    def rows(self, /) -> int: ...
    @property
    def cols(self, /) -> int: ...
    @property
    def entries(self, /) -> int: ...
    @property
    def format(self, /) -> _Format: ...
    @property
    def field(self, /) -> _Field: ...
    @property
    def symmetry(self, /) -> _Symmetry: ...
    @property
    def has_symmetry(self, /) -> bool: ...

    #
    def __init__(self, /, **kwargs: Unpack[_MMFileKwargs]) -> None: ...

    #
    @overload
    def read(self, /, source: FileLike[bytes], *, spmatrix: onp.ToTrue = True) -> onp.ArrayND[npc.number] | coo_array: ...
    @overload
    def read(self, /, source: FileLike[bytes], *, spmatrix: onp.ToFalse) -> onp.ArrayND[npc.number] | coo_matrix: ...

    #
    def write(
        self,
        /,
        target: FileLike[bytes],
        a: spmatrix | sparray | onp.ToArrayND,
        comment: str = "",
        field: _Field | None = None,
        precision: int | None = None,
        symmetry: _Symmetry | None = None,
    ) -> None: ...

    #
    @classmethod
    def info(cls, /, source: FileLike[bytes]) -> _Info: ...
    @staticmethod
    def reader() -> None: ...
    @staticmethod
    def writer() -> None: ...

#
def asstr(s: object) -> str: ...
