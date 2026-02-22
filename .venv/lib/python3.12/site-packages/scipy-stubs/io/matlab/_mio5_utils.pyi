# defined in scipy/io/matlab/_mio5_utils.pyx

from collections.abc import Iterable
from typing import Any, ClassVar, Final, Literal, Never, Self, TypeAlias
from typing_extensions import CapsuleType

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._miobase import MatVarReader
from ._streams import GenericStream, _FileLike

_bint: TypeAlias = Literal[0, 1] | bool  # noqa: PYI042

###

# NOTE: These apprear to be broken, and will always raise `TypeError: no default __reduce__ due to non-trivial __cinit__`
def __reduce_cython__(self: Never, /) -> Never: ...  # undocumented
def __setstate_cython__(self: Never, pyx_state: Never, /) -> None: ...  # undocumented

swapped_code: Final[Literal[">", "<"]] = ...  # undocumented  # ">" sys.byteorder == "little" else "<"

def byteswap_u4(u4: np.uint32) -> np.uint32: ...  # undocumented

class VarHeader5:  # undocumented
    # cdef readonly object name
    name: Final[object]
    # cdef readonly int mclass
    mclass: Final[int]
    # cdef readonly object dims
    dims: Final[Iterable[int | npc.integer]]
    # cdef readonly int is_logical
    is_logical: Final[_bint]
    # cdef public int is_global
    is_global: _bint
    # cdef readonly size_t nzmax
    nzmax: Final[int]

    def __reduce_cython__(self) -> tuple[Any, ...]: ...
    def __setstate_cython__(self, /, state: tuple[object, ...]) -> None: ...

    #
    def set_dims(self, /, dims: object) -> None: ...

class VarReader5:  # undocumented
    __pyx_vtable__: ClassVar[CapsuleType] = ...

    # cdef public int is_swapped, little_endian
    is_swapped: _bint
    little_endian: _bint

    def __new__(cls, preader: MatVarReader) -> Self: ...
    def set_stream(self, /, fobj: GenericStream | _FileLike) -> None: ...
    def read_tag(self, /) -> tuple[int, int, str | None]: ...
    def read_numeric(self, /, copy: _bint = True, nnz: int = -1) -> onp.Array1D[Any]: ...
    def read_full_tag(self, /) -> tuple[np.uint32, np.uint32]: ...
    def read_header(self, /, check_stream_limit: _bint) -> VarHeader5: ...
    def array_from_header(self, /, header: VarHeader5, process: _bint = 1) -> Any: ...  # squeezed ndarray or sparse csc_array
    def shape_from_header(self, /, header: VarHeader5) -> tuple[int, ...]: ...
    def read_real_complex(self, /, header: VarHeader5) -> onp.ArrayND[np.float64 | np.complex128]: ...
    def read_char(self, /, header: VarHeader5) -> onp.ArrayND[np.str_]: ...
    def read_cells(self, /, header: VarHeader5) -> onp.ArrayND[np.object_]: ...
    def read_fieldnames(self, /) -> list[str]: ...
    def read_struct(self, /, header: VarHeader5) -> onp.ArrayND[np.object_]: ...
    def read_opaque(self, /, hdr: VarHeader5) -> onp.Array1D[np.void]: ...
