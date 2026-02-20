from collections.abc import Mapping, Sequence
from typing import IO, Any, Final, Generic, Literal, Self, SupportsIndex, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import numpy.typing as npt
import optype.numpy as onp

from ._typing import FileLike
from scipy._typing import ExitMixin

__all__ = ["netcdf_file", "netcdf_variable"]

###

_ShapeT_co = TypeVar("_ShapeT_co", covariant=True, bound=tuple[int, ...], default=tuple[Any, ...])
_ScalarT = TypeVar("_ScalarT", bound=np.generic)
_ScalarT_co = TypeVar("_ScalarT_co", covariant=True, bound=np.generic, default=Any)

_FileModeRWA: TypeAlias = Literal["r", "w", "a"]
_TypeCode: TypeAlias = Literal["b", "c", "h", "i", "f", "d"]
_TypeSize: TypeAlias = Literal[1, 2, 4, 8]
_TypeSpec: TypeAlias = tuple[_TypeCode, _TypeSize]
_TypeNC: TypeAlias = Literal[
    b"\x00\x00\x00\x01",
    b"\x00\x00\x00\x02",
    b"\x00\x00\x00\x03",
    b"\x00\x00\x00\x04",
    b"\x00\x00\x00\x05",
    b"\x00\x00\x00\x06",
    b"\x00\x00\x00\n",
    b"\x00\x00\x00\x0b",
    b"\x00\x00\x00\x0c",
]
_TypeFill: TypeAlias = Literal[
    b"\x81",
    b"\x00",
    b"\x80\x01",
    b"\x80\x00\x00\x01",
    b"\x7c\xf0\x00\x00",
    b"\x47\x9e\x00\x00\x00\x00\x00\x00",
]  # fmt: skip

###

IS_PYPY: Final[bool] = ...

ABSENT: Final = b"\x00\x00\x00\x00\x00\x00\x00\x00"
ZERO: Final = b"\x00\x00\x00\x00"
NC_BYTE: Final[_TypeNC] = b"\x00\x00\x00\x01"
NC_CHAR: Final[_TypeNC] = b"\x00\x00\x00\x02"
NC_SHORT: Final[_TypeNC] = b"\x00\x00\x00\x03"
NC_INT: Final[_TypeNC] = b"\x00\x00\x00\x04"
NC_FLOAT: Final[_TypeNC] = b"\x00\x00\x00\x05"
NC_DOUBLE: Final[_TypeNC] = b"\x00\x00\x00\x06"
NC_DIMENSION: Final[_TypeNC] = b"\x00\x00\x00\n"
NC_VARIABLE: Final[_TypeNC] = b"\x00\x00\x00\x0b"
NC_ATTRIBUTE: Final[_TypeNC] = b"\x00\x00\x00\x0c"
FILL_BYTE: Final[_TypeFill] = b"\x81"
FILL_CHAR: Final[_TypeFill] = b"\x00"
FILL_SHORT: Final[_TypeFill] = b"\x80\x01"
FILL_INT: Final[_TypeFill] = b"\x80\x00\x00\x01"
FILL_FLOAT: Final[_TypeFill] = b"\x7c\xf0\x00\x00"
FILL_DOUBLE: Final[_TypeFill] = b"\x47\x9e\x00\x00\x00\x00\x00\x00"

TYPEMAP: Final[dict[_TypeNC, _TypeSpec]] = ...
FILLMAP: Final[dict[_TypeNC, _TypeFill]]
REVERSE: Final[dict[_TypeSpec, _TypeNC]]

class netcdf_file(ExitMixin):
    fp: Final[IO[bytes]]
    filename: Final[str]
    use_mmap: Final[bool]
    mode: Final[_FileModeRWA]
    version_byte: Final[int]
    maskandscale: Final[bool]
    dimensions: Final[dict[str, int]]
    variables: Final[dict[str, netcdf_variable]]

    def __init__(
        self,
        /,
        filename: FileLike[bytes],
        mode: _FileModeRWA = "r",
        mmap: bool | None = None,
        version: int = 1,
        maskandscale: bool = False,
    ) -> None: ...
    def __del__(self, /) -> None: ...
    def __enter__(self, /) -> Self: ...

    #
    def createDimension(self, /, name: str, length: int) -> None: ...

    #
    @overload
    def createVariable(
        self, /, name: str, type: onp.ToDType[_ScalarT], dimensions: Sequence[str]
    ) -> NetCDFVariable[tuple[Any, ...], _ScalarT]: ...
    @overload
    def createVariable(self, /, name: str, type: npt.DTypeLike, dimensions: Sequence[str]) -> NetCDFVariable: ...

    #
    def flush(self, /) -> None: ...
    def sync(self, /) -> None: ...
    def close(self, /) -> None: ...

class netcdf_variable(Generic[_ShapeT_co, _ScalarT_co]):
    data: onp.Array[_ShapeT_co, _ScalarT_co]
    dimensions: Final[Sequence[str]]
    maskandscale: Final[bool]

    #
    @property
    def isrec(self, /) -> bool: ...
    @property
    def shape(self, /) -> _ShapeT_co: ...

    #
    def __init__(
        self,
        /,
        data: onp.Array[_ShapeT_co, _ScalarT_co],
        typecode: str,
        size: int,
        shape: tuple[int, ...] | list[int],
        dimensions: Sequence[str],
        attributes: Mapping[str, object] | None = None,
        maskandscale: bool = False,
    ) -> None: ...

    #
    def __getitem__(
        self, /, index: SupportsIndex | slice | tuple[SupportsIndex | slice, ...]
    ) -> _ScalarT | onp.ArrayND[_ScalarT]: ...
    def __setitem__(
        self: netcdf_variable[tuple[int, ...], _ScalarT], /, index: object, data: _ScalarT | onp.ArrayND[_ScalarT]
    ) -> None: ...

    #
    def assignValue(self, /, value: object) -> None: ...
    def getValue(self, /) -> _ScalarT_co: ...
    def typecode(self, /) -> str: ...
    def itemsize(self, /) -> int: ...

NetCDFFile = netcdf_file
NetCDFVariable = netcdf_variable
