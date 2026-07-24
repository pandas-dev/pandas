from typing import Any, ClassVar, Literal, TypedDict, Unpack, overload, type_check_only
from typing_extensions import deprecated

import optype as op
import optype.numpy as onp

from ._fast_matrix_market import mminfo, mmread, mmwrite
from ._typing import FileLike
from scipy.sparse import coo_array, coo_matrix, sparray, spmatrix

__all__ = ["MMFile", "mminfo", "mmread", "mmwrite"]

###

type _Format = Literal["coordinate", "array"]
type _Field = Literal["real", "complex", "pattern", "integer"]
type _Symmetry = Literal["general", "symmetric", "skew-symmetric", "hermitian"]
type _Info = tuple[int, int, int, _Format, _Field, _Symmetry]

type _NoValueType = op.JustObject

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

    # dtype is either intp, uint64, float64, or complex128, depending on the field
    @overload
    @deprecated("The default value for `spmatrix` is changing to False in v1.20.")
    def read(self, /, source: FileLike[bytes], *, spmatrix: _NoValueType = ...) -> onp.Array2D[Any] | coo_matrix[Any]: ...
    @overload
    def read(self, /, source: FileLike[bytes], *, spmatrix: Literal[True]) -> onp.Array2D[Any] | coo_matrix[Any]: ...
    @overload
    def read(
        self, /, source: FileLike[bytes], *, spmatrix: Literal[False]
    ) -> onp.Array2D[Any] | coo_array[Any, tuple[int, int]]: ...

    #
    def write(
        self,
        /,
        target: FileLike[bytes],
        a: spmatrix | sparray | onp.ToArray2D,
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
