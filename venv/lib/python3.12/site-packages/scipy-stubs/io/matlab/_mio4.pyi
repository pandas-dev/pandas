from collections.abc import Mapping, Sequence
from typing import IO, Final, Generic, Literal, LiteralString, Protocol, TypeAlias, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._miobase import MatFileReader
from scipy.sparse import coo_matrix, sparray, spmatrix

__all__ = [
    "SYS_LITTLE_ENDIAN",
    "MatFile4Reader",
    "MatFile4Writer",
    "VarHeader4",
    "VarReader4",
    "VarWriter4",
    "arr_to_2d",
    "mclass_info",
    "mdtypes_template",
    "miDOUBLE",
    "miINT16",
    "miINT32",
    "miSINGLE",
    "miUINT8",
    "miUINT16",
    "mxCHAR_CLASS",
    "mxFULL_CLASS",
    "mxSPARSE_CLASS",
    "np_to_mtypes",
    "order_codes",
]

_OnedAs: TypeAlias = Literal["row", "col"]
_MDType: TypeAlias = Literal[0, 1, 2, 3, 4, 5]
_MClass: TypeAlias = Literal[0, 1, 2]

SYS_LITTLE_ENDIAN: Final[bool] = ...

miDOUBLE: Final = 0
miSINGLE: Final = 1
miINT32: Final = 2
miINT16: Final = 3
miUINT16: Final = 4
miUINT8: Final = 5

mxFULL_CLASS: Final = 0
mxCHAR_CLASS: Final = 1
mxSPARSE_CLASS: Final = 2

mdtypes_template: Final[dict[int | str, str | list[tuple[str, str]]]]
np_to_mtypes: Final[dict[str, int]]
order_codes: Final[dict[int, str]]
mclass_info: Final[dict[int, str]]

_DT = TypeVar("_DT", bound=np.dtype[np.generic])
_DT_co = TypeVar("_DT_co", covariant=True, bound=np.dtype[np.generic], default=np.dtype[np.generic])

@type_check_only
class _SupportsVarHeader(Protocol[_DT_co]):
    @property
    def name(self, /) -> str: ...
    @property
    def dtype(self, /) -> _DT_co: ...
    @property
    def mclass(self, /) -> _MClass: ...
    @property
    def dims(self, /) -> int: ...
    @property
    def is_global(self, /) -> bool: ...

# NOTE: `VarHeader4` is assignable to `_SupportsVarHeader`
class VarHeader4(Generic[_DT_co]):
    name: LiteralString
    dtype: _DT_co
    dims: Final[int]
    mclass: Final[_MClass]
    is_complex: Final[bool]
    is_logical: Final[bool]
    is_global: Final[bool]

    def __init__(self, /, name: str, dtype: _DT_co, mclass: _MClass, dims: int, is_complex: bool) -> None: ...

class VarReader4:
    file_reader: Final[MatFileReader]
    mat_stream: Final[IO[bytes]]
    dtypes: Final[Mapping[str, np.dtype[np.generic]]]
    chars_as_strings: Final[bool]
    squeeze_me: Final[bool]

    def __init__(self, /, file_reader: MatFileReader) -> None: ...
    def read_header(self, /) -> tuple[VarHeader4, int]: ...
    def array_from_header(
        self, /, hdr: _SupportsVarHeader[_DT], process: bool = True
    ) -> np.ndarray[tuple[int, ...], _DT] | coo_matrix: ...
    def read_sub_array(self, /, hdr: _SupportsVarHeader[_DT], copy: bool = True) -> np.ndarray[tuple[int, ...], _DT]: ...
    def read_full_array(self, /, hdr: _SupportsVarHeader[_DT]) -> np.ndarray[tuple[int, ...], _DT]: ...
    def read_char_array(self, /, hdr: _SupportsVarHeader[_DT]) -> np.ndarray[tuple[int, ...], _DT]: ...
    def read_sparse_array(self, /, hdr: _SupportsVarHeader) -> coo_matrix: ...
    def shape_from_header(self, /, hdr: _SupportsVarHeader) -> tuple[int, ...]: ...

class MatFile4Reader(MatFileReader):
    def __init__(self, /, mat_stream: IO[bytes], *args: bool | str, **kwargs: bool | str) -> None: ...
    def initialize_read(self, /) -> None: ...
    def read_var_header(self, /) -> tuple[VarHeader4, int]: ...
    def read_var_array(self, /, header: _SupportsVarHeader[_DT], process: bool = True) -> np.ndarray[tuple[int, ...], _DT]: ...
    def get_variables(self, /, variable_names: str | Sequence[str] | None = None) -> dict[str, onp.ArrayND]: ...
    def list_variables(self, /) -> list[tuple[str, tuple[int, ...], str]]: ...

def arr_to_2d(arr: onp.ArrayND, oned_as: _OnedAs = "row") -> onp.Array2D: ...

class VarWriter4:
    file_stream: IO[bytes]
    oned_as: _OnedAs
    def __init__(self, /, file_writer: MatFile4Writer) -> None: ...
    def write_bytes(self, /, arr: onp.ArrayND) -> None: ...
    def write_string(self, /, s: str) -> None: ...
    def write_header(self, /, name: str, shape: Sequence[int], P: _MDType = 0, T: _MClass = 0, imagf: int = 0) -> None: ...
    def write(self, /, arr: onp.ToArrayND, name: str) -> None: ...
    def write_numeric(self, /, arr: onp.ArrayND[npc.number], name: str) -> None: ...
    def write_char(self, /, arr: onp.ArrayND[np.character], name: str) -> None: ...
    def write_sparse(self, /, arr: spmatrix | sparray, name: str) -> None: ...

class MatFile4Writer:
    file_stream: IO[bytes]
    oned_as: _OnedAs
    def __init__(self, /, file_stream: IO[bytes], oned_as: _OnedAs | None = None) -> None: ...
    def put_variables(self, /, mdict: Mapping[str, onp.ArrayND], write_header: bool | None = None) -> None: ...
