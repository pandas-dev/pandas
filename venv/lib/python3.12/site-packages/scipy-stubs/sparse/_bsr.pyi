from collections.abc import Sequence
from typing import Any, ClassVar, Generic, Literal, Never, TypeAlias, overload, type_check_only
from typing_extensions import TypeIs, TypeVar, override

import numpy as np
import numpy.typing as npt
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._compressed import _cs_matrix
from ._coo import coo_array, coo_matrix
from ._data import _minmax_mixin
from ._matrix import spmatrix
from ._typing import _ToShape2D

__all__ = ["bsr_array", "bsr_matrix", "isspmatrix_bsr"]

###

_T = TypeVar("_T")
_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool_, default=Any, covariant=True)

_ToMatrixPy: TypeAlias = Sequence[_T] | Sequence[Sequence[_T]]
_ToMatrix: TypeAlias = _spbase[_ScalarT] | onp.CanArrayND[_ScalarT] | Sequence[onp.CanArrayND[_ScalarT]] | _ToMatrixPy[_ScalarT]

_ToData2: TypeAlias = tuple[onp.ArrayND[_ScalarT], onp.ArrayND[npc.integer]]
_ToData3: TypeAlias = tuple[onp.ArrayND[_ScalarT], onp.ArrayND[npc.integer], onp.ArrayND[npc.integer]]
_ToData: TypeAlias = _ToData2[_ScalarT] | _ToData3[_ScalarT]

###

class _bsr_base(_cs_matrix[_ScalarT_co, tuple[int, int]], _minmax_mixin[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    _format: ClassVar = "bsr"

    data: onp.Array3D[_ScalarT_co]

    @property
    @override
    def format(self, /) -> Literal["bsr"]: ...
    @property
    @override
    def ndim(self, /) -> Literal[2]: ...
    @property
    @override
    def shape(self, /) -> tuple[int, int]: ...
    @property
    def blocksize(self, /) -> tuple[int, int]: ...

    #
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    # These raise `NotImplementedError` at runtime (so why inherit from `IndexMixin` then...?)
    @override
    def __getitem__(self, key: Never, /) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @override
    def __setitem__(self, key: Never, val: object, /) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]

#
class bsr_array(_bsr_base[_ScalarT_co], sparray[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_array[_ScalarT_co, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_array[_ScalarT, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> bsr_array[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> bsr_array[np.float64]: ...

    # NOTE: keep in sync with `bsr_matrix.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co] | _ToData[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like | None
    def __init__(
        self: bsr_array[np.float64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (positional)
    def __init__(
        self: bsr_array[np.bool_],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (keyword)
    def __init__(
        self: bsr_array[np.bool_],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (positional)
    def __init__(
        self: bsr_array[np.int64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (keyword)
    def __init__(
        self: bsr_array[np.int64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: bsr_array[np.complex128],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (keyword)
    def __init__(
        self: bsr_array[np.complex128],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: bsr_array[np.bool_],
        /,
        arg1: _ToMatrixPy[bool],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.int, dtype: int-like | None
    def __init__(
        self: bsr_array[np.int_],
        /,
        arg1: _ToMatrixPy[op.JustInt],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: bsr_array[np.float64],
        /,
        arg1: _ToMatrixPy[op.JustFloat],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: bsr_array[np.complex128],
        /,
        arg1: _ToMatrixPy[op.JustComplex],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <unknown>
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...

class bsr_matrix(_bsr_base[_ScalarT_co], spmatrix[_ScalarT_co], Generic[_ScalarT_co]):
    # NOTE: These two methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_matrix[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_matrix[_ScalarT]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> bsr_matrix[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> bsr_matrix[np.float64]: ...

    # NOTE: keep in sync with `bsr_array.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co] | _ToData[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: None
    def __init__(
        self: bsr_matrix[np.float64],
        /,
        arg1: _ToShape2D,
        shape: None = None,
        dtype: None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: bsr_matrix[np.bool_],
        /,
        arg1: _ToMatrixPy[bool],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.int, dtype: int-like | None
    def __init__(
        self: bsr_matrix[np.int_],
        /,
        arg1: _ToMatrixPy[op.JustInt],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: bsr_matrix[np.float64],
        /,
        arg1: _ToMatrixPy[op.JustFloat],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: bsr_matrix[np.complex128],
        /,
        arg1: _ToMatrixPy[op.JustComplex],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <unknown>
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        blocksize: tuple[int, int] | None = None,
        *,
        maxprint: int | None = None,
    ) -> None: ...

#
def isspmatrix_bsr(x: object) -> TypeIs[bsr_matrix]: ...
