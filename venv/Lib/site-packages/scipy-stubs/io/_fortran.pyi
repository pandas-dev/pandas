from typing import Any, Literal, Self, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import ExitMixin
from scipy.io._typing import FileLike

__all__ = ["FortranEOFError", "FortranFile", "FortranFormattingError"]

###

class FortranEOFError(TypeError, OSError): ...
class FortranFormattingError(TypeError, OSError): ...

class FortranFile(ExitMixin):
    def __init__(
        self,
        /,
        filename: FileLike[bytes],
        mode: Literal["r", "w"] = "r",
        header_dtype: onp.AnyDType = np.uint32,  # noqa: PYI011
    ) -> None: ...

    #
    def __enter__(self, /) -> Self: ...
    def close(self, /) -> None: ...

    #
    def write_record(self, /, *items: onp.ToArrayND) -> None: ...

    #
    @overload
    def read_record[ScalarT: np.generic](
        self, dtype0: onp.ToDType[ScalarT], /, *, dtype: None = None
    ) -> onp.Array1D[ScalarT]: ...
    @overload
    def read_record[ScalarT: np.generic](self, /, *, dtype: onp.ToDType[ScalarT]) -> onp.Array1D[ScalarT]: ...
    @overload
    def read_record(self, dtype0: str | type, /, *, dtype: None = None) -> onp.Array1D: ...
    @overload
    def read_record(self, /, *, dtype: str | type) -> onp.Array1D: ...
    @overload
    def read_record[T1: np.generic, T2: np.generic](
        self, dtype0: onp.ToDType[T1], dtype1: onp.ToDType[T2], /, *, dtype: None = None
    ) -> tuple[onp.Array1D[T1], onp.Array1D[T2]]: ...
    @overload
    def read_record[T1: np.generic, T2: np.generic](
        self, dtype0: onp.ToDType[T1], /, *, dtype: onp.ToDType[T2]
    ) -> tuple[onp.Array1D[T1], onp.Array1D[T2]]: ...
    @overload
    def read_record(
        self,
        dtype0: np.dtype[Any] | str | type,
        dtype1: np.dtype[Any] | str | type,
        /,
        *dtypes: np.dtype[Any] | str | type,
        dtype: np.dtype[Any] | str | type | None = None,
    ) -> tuple[onp.Array1D, ...]: ...
    @overload
    def read_record(
        self, dtype0: np.dtype[Any] | str | type, /, *, dtype: np.dtype[Any] | str | type
    ) -> tuple[onp.Array1D, ...]: ...

    #
    @overload
    def read_ints(self, /) -> onp.Array1D[np.int32]: ...
    @overload
    def read_ints[IntegerT: npc.integer](self, /, dtype: onp.ToDType[IntegerT]) -> onp.Array1D[IntegerT]: ...
    @overload
    def read_ints(self, /, dtype: str | type) -> onp.Array1D: ...

    #
    @overload
    def read_reals(self, /) -> onp.Array1D[np.float64]: ...
    @overload
    def read_reals[FloatingT: npc.floating](self, /, dtype: onp.ToDType[FloatingT]) -> onp.Array1D[FloatingT]: ...
    @overload
    def read_reals(self, /, dtype: str | type) -> onp.Array1D: ...
