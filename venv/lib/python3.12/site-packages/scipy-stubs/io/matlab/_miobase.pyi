import abc
from collections.abc import Callable, Mapping
from typing import IO, Any, Final, Literal, Protocol, TypeVar, type_check_only

import numpy as np
import numpy.typing as npt
import optype.numpy as onp

from scipy.io._typing import ByteOrder, FileName

__all__ = ["MatReadError", "MatReadWarning", "MatWriteError", "MatWriteWarning"]

_FT = TypeVar("_FT", bound=Callable[..., object])

@type_check_only
class _Decorator(Protocol):
    def __call__(self, f: _FT, /) -> _FT: ...

###

doc_dict: Final[dict[str, str]] = ...
docfiller: Final[_Decorator] = ...

class MatReadError(Exception): ...
class MatWriteError(Exception): ...
class MatReadWarning(UserWarning): ...
class MatWriteWarning(UserWarning): ...

class MatVarReader:
    def __init__(self, /, file_reader: MatFileReader) -> None: ...
    @abc.abstractmethod
    def read_header(self, /) -> dict[str, Any]: ...
    @abc.abstractmethod
    def array_from_header(self, /, header: Mapping[str, object]) -> onp.ArrayND: ...

class MatFileReader:
    mat_stream: Final[IO[bytes]]
    dtypes: Final[Mapping[str, np.dtype[Any]]]
    byte_order: Final[ByteOrder]
    struct_as_record: Final[bool]
    verify_compressed_data_integrity: Final[bool]
    simplify_cells: Final[bool]
    mat_dtype: bool
    squeeze_me: bool
    chars_as_strings: bool

    def __init__(
        self,
        /,
        mat_stream: IO[bytes],
        byte_order: ByteOrder | None = None,
        mat_dtype: bool = False,
        squeeze_me: bool = False,
        chars_as_strings: bool = True,
        matlab_compatible: bool = False,
        struct_as_record: bool = True,
        verify_compressed_data_integrity: bool = True,
        simplify_cells: bool = False,
    ) -> None: ...
    def set_matlab_compatible(self, /) -> None: ...
    def guess_byte_order(self, /) -> ByteOrder: ...
    def end_of_stream(self, /) -> bool: ...

def _get_matfile_version(fileobj: IO[bytes]) -> tuple[Literal[1, 2], int]: ...
def convert_dtypes(dtype_template: Mapping[str, str], order_code: ByteOrder) -> dict[str, np.dtype[Any]]: ...
def read_dtype(mat_stream: IO[bytes], a_dtype: npt.DTypeLike) -> onp.ArrayND: ...
def matfile_version(file_name: FileName, *, appendmat: bool = True) -> tuple[Literal[0, 1, 2], int]: ...
def matdims(arr: onp.ArrayND, oned_as: Literal["column", "row"] = "column") -> tuple[Any, ...]: ...
def arr_dtype_number(arr: onp.HasDType[np.dtype[np.character]], num: int | str) -> np.dtype[np.character]: ...
def arr_to_chars(arr: onp.ArrayND[np.str_]) -> onp.ArrayND[np.str_]: ...

get_matfile_version = matfile_version
