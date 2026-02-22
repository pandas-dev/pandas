from typing import Literal, Self, TypeAlias, TypedDict, final, overload, type_check_only
from typing_extensions import TypeVar, Unpack

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import ExitMixin
from scipy.io._typing import FileLike

__all__ = ["FortranEOFError", "FortranFile", "FortranFormattingError"]

_IntegerT = TypeVar("_IntegerT", bound=npc.integer)
_FloatingT = TypeVar("_FloatingT", bound=npc.floating)

_FileModeRW: TypeAlias = Literal["r", "w"]

@final
@type_check_only
class _DTypeKwargs(TypedDict, total=False):
    dtype: onp.ToDType

###

class FortranEOFError(TypeError, OSError): ...
class FortranFormattingError(TypeError, OSError): ...

class FortranFile(ExitMixin):
    def __init__(
        self,
        /,
        filename: FileLike[bytes],
        mode: _FileModeRW = "r",
        header_dtype: npt.DTypeLike = np.uint32,  # noqa: PYI011
    ) -> None: ...

    #
    def __enter__(self, /) -> Self: ...
    def close(self, /) -> None: ...

    #
    def write_record(self, /, *items: onp.ToArrayND) -> None: ...
    def read_record(self, /, *dtypes: onp.ToDType, **kwargs: Unpack[_DTypeKwargs]) -> onp.Array1D[np.void]: ...

    #
    @overload
    def read_ints(self, /) -> onp.Array1D[np.int32]: ...
    @overload
    def read_ints(self, /, dtype: onp.ToDType[_IntegerT]) -> onp.Array1D[_IntegerT]: ...
    @overload
    def read_ints(self, /, dtype: str | type) -> onp.Array1D: ...

    #
    @overload
    def read_reals(self, /) -> onp.Array1D[np.float64]: ...
    @overload
    def read_reals(self, /, dtype: onp.ToDType[_FloatingT]) -> onp.Array1D[_FloatingT]: ...
    @overload
    def read_reals(self, /, dtype: str | type) -> onp.Array1D: ...
