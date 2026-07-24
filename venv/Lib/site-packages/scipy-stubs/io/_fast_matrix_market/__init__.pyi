import io
from typing import Any, Final, Literal, Unpack, overload, override, type_check_only
from typing_extensions import TypedDict, deprecated

import optype as op
import optype.numpy as onp

from scipy.io._typing import FileLike
from scipy.sparse import coo_array, coo_matrix, sparray, spmatrix

__all__ = ["mminfo", "mmread", "mmwrite"]

type _NoValueType = op.JustObject

type _Format = Literal["coordinate", "array"]
type _Field = Literal["real", "complex", "pattern", "integer"]
type _Symmetry = Literal["general", "symmetric", "skew-symmetric", "hermitian"]
type _Info = tuple[int, int, int, _Format, _Field, _Symmetry]

@type_check_only
class _TextToBytesWrapperKwargs(TypedDict, total=False):
    buffer_size: int

###

PARALLELISM: Final = 0
ALWAYS_FIND_SYMMETRY: Final = False

class _TextToBytesWrapper(io.BufferedReader):
    encoding: Final[str]
    errors: Final[str]

    def __init__(
        self,
        /,
        text_io_buffer: io.TextIOBase,
        encoding: str | None = None,
        errors: str | None = None,
        **kwargs: Unpack[_TextToBytesWrapperKwargs],
    ) -> None: ...
    @override
    def read(self, /, size: int | None = -1) -> bytes: ...
    @override
    def read1(self, /, size: int = -1) -> bytes: ...
    @override
    def peek(self, /, size: int = -1) -> bytes: ...
    @override
    # pyrefly: ignore [bad-override]
    def seek(self, /, offset: int, whence: int = 0) -> None: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]

#
@overload
@deprecated(
    "The default value for `spmatrix` is changing to `False` in v1.20. That means the default return type will be a sparse "
    "array. Unless you use * instead of @, ** for matrix power, or you depend on 2D shapes from e.g. `A.sum(axis=0)` it may not "
    "matter to you. See the spmatrix to sparray migration guide for details. "
    "https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html"
)
def mmread(source: FileLike[bytes], *, spmatrix: _NoValueType = ...) -> onp.Array2D | coo_matrix: ...
@overload
def mmread(source: FileLike[bytes], *, spmatrix: Literal[True]) -> onp.Array2D | coo_matrix: ...
@overload
def mmread(source: FileLike[bytes], *, spmatrix: Literal[False]) -> onp.Array2D | coo_array[Any, tuple[int, int]]: ...

# these defaults are different in `io._mmio`, so we don't specify them to avoid duplicate definitions
def mmwrite(
    target: FileLike[bytes],
    a: onp.ToArray2D | spmatrix | sparray,
    comment: str | None = ...,  # stubdefaulter: ignore[missing-default]
    field: _Field | None = None,
    precision: int | None = None,
    symmetry: _Symmetry | Literal["AUTO"] | None = ...,  # stubdefaulter: ignore[missing-default]
) -> None: ...

#
def mminfo(source: FileLike[bytes]) -> _Info: ...
