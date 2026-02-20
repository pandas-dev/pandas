# NOTE: Adding `@override` to `@overload`ed methods will crash stubtest (mypy 1.13.0)
# mypy: disable-error-code="misc, override"

from collections.abc import Iterable, Sequence
from typing import Any, ClassVar, Generic, Literal, Never, Self, TypeAlias, overload, type_check_only
from typing_extensions import TypeIs, TypeVar, override

import numpy as np
import numpy.typing as npt
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._coo import coo_array, coo_matrix
from ._index import IndexMixin
from ._matrix import spmatrix
from ._typing import _ToShape1D, _ToShape2D

__all__ = ["dok_array", "dok_matrix", "isspmatrix_dok"]

###

_T = TypeVar("_T")
_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int] | tuple[int, int])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int] | tuple[int, int], default=tuple[int, int], covariant=True)

_1D: TypeAlias = tuple[int]  # noqa: PYI042
_2D: TypeAlias = tuple[int, int]  # noqa: PYI042
# workaround for the typing-spec non-conformance regarding overload behavior of mypy and pyright
_NoD: TypeAlias = tuple[Never] | tuple[Never, Never]
_AnyD: TypeAlias = tuple[Any, ...]

_ToMatrix: TypeAlias = _spbase[_ScalarT] | onp.CanArrayND[_ScalarT] | Sequence[onp.CanArrayND[_ScalarT]] | _ToMatrixPy[_ScalarT]
_ToMatrixPy: TypeAlias = Sequence[_T] | Sequence[Sequence[_T]]

_ToKey1D: TypeAlias = onp.ToJustInt | tuple[onp.ToJustInt]
_ToKey2D: TypeAlias = tuple[onp.ToJustInt, onp.ToJustInt]

_ToKeys1: TypeAlias = Iterable[_ToKey1D]
_ToKeys2: TypeAlias = Iterable[_ToKey2D]
_ToKeys: TypeAlias = Iterable[_ToKey1D | _ToKey2D]

_C2T = TypeVar("_C2T", bound=_dok_base[np.float64, _2D])

###

class _dok_base(  # pyright: ignore[reportIncompatibleMethodOverride]
    _spbase[_ScalarT_co, _ShapeT_co],
    IndexMixin[_ScalarT_co, _ShapeT_co],
    dict[tuple[Any, ...], _ScalarT_co | Any],
    Generic[_ScalarT_co, _ShapeT_co],
):
    _format: ClassVar = "dok"
    _allow_nd: ClassVar = 1, 2

    dtype: np.dtype[_ScalarT_co]

    @property
    @override
    def format(self, /) -> Literal["dok"]: ...
    @property
    @override
    def ndim(self, /) -> Literal[1, 2]: ...
    @property
    @override
    def shape(self, /) -> _ShapeT_co: ...

    #
    def __init__(
        self,
        /,
        arg1: onp.ToComplexND,
        shape: _ToShape1D | _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def __len__(self, /) -> int: ...
    @override
    def __delitem__(self: _dok_base[Any, _ShapeT], key: _ShapeT, /) -> None: ...  # pyrefly: ignore[bad-override]

    #
    @override
    @overload
    def setdefault(  # pyrefly: ignore[bad-override]
        self: _dok_base[_ScalarT, _ShapeT], key: _ShapeT, default: _T, /
    ) -> _ScalarT | _T: ...
    @overload
    def setdefault(self: _dok_base[_ScalarT, _ShapeT], key: _ShapeT, default: None = None, /) -> _ScalarT | None: ...  # pyright: ignore[reportIncompatibleMethodOverride]

    #
    @override
    @overload
    def get(  # pyrefly: ignore[bad-override]
        self: _dok_base[_ScalarT, _ShapeT], /, key: _ShapeT, default: _T
    ) -> _ScalarT | _T: ...
    @overload
    def get(self: _dok_base[_ScalarT, _ShapeT], /, key: _ShapeT, default: float = 0.0) -> _ScalarT | float: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

    #
    @override
    def __or__(self, other: Never, /) -> Never: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
    @override
    def __ror__(self, other: Never, /) -> Never: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
    @override
    def __ior__(self, other: Never, /) -> Self: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
    @override
    def update(self, /, val: Never) -> Never: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]

    #
    @override
    def count_nonzero(self, /, axis: None = None) -> int: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]

    #
    @override
    @overload
    @classmethod
    def fromkeys(  # pyrefly: ignore[bad-override]
        cls: type[_dok_base[np.bool_, _2D]], iterable: _ToKeys2, v: onp.ToBool, /
    ) -> _dok_base[np.bool_, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[np.bool_, _1D]], iterable: _ToKeys1, v: onp.ToBool, /) -> _dok_base[np.bool_, _1D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[_ScalarT, _2D]], iterable: _ToKeys2, v: _ScalarT, /) -> _dok_base[_ScalarT, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[_ScalarT, _1D]], iterable: _ToKeys1, v: _ScalarT, /) -> _dok_base[_ScalarT, _1D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[np.int_, _2D]], iterable: _ToKeys2, v: op.JustInt = 1, /) -> _dok_base[np.int_, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[np.int_, _1D]], iterable: _ToKeys1, v: op.JustInt = 1, /) -> _dok_base[np.int_, _1D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[np.float64, _2D]], iterable: _ToKeys2, v: op.JustFloat, /) -> _dok_base[np.float64, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_dok_base[np.float64, _1D]], iterable: _ToKeys1, v: op.JustFloat, /) -> _dok_base[np.float64, _1D]: ...
    @overload
    @classmethod
    def fromkeys(
        cls: type[_dok_base[np.complex128, _NoD]], iterable: _ToKeys, v: op.JustComplex, /
    ) -> _dok_base[np.complex128, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[_C2T], iterable: _ToKeys2, v: op.JustComplex, /) -> _C2T: ...
    @overload
    @classmethod
    def fromkeys(  # pyright: ignore[reportIncompatibleMethodOverride]
        cls: type[_dok_base[np.complex128, _1D]], iterable: _ToKeys1, v: op.JustComplex, /
    ) -> _dok_base[np.complex128, _1D]: ...

#
class dok_array(_dok_base[_ScalarT_co, _ShapeT_co], sparray[_ScalarT_co, _ShapeT_co], Generic[_ScalarT_co, _ShapeT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_array[_ScalarT_co, _2D]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_array[_ScalarT, _2D]: ...
    @type_check_only
    def __assoc_as_float32__(self, /) -> dok_array[np.float32, _ShapeT_co]: ...
    @type_check_only
    def __assoc_as_float64__(self, /) -> dok_array[np.float64, _ShapeT_co]: ...

    # NOTE: keep the 2d overloads in sync with `dok_matrix.__init__`
    # TODO(jorenham): Overloads for specific shape types.
    @overload  # sparse or dense (know dtype & shape), dtype: None
    def __init__(
        self,
        /,
        arg1: _spbase[_ScalarT_co, _ShapeT_co] | onp.CanArrayND[_ScalarT_co, _ShapeT_co],
        shape: _ShapeT_co | None = None,
        dtype: onp.ToDType[_ScalarT_co] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d array-like (know dtype), dtype: None
    def __init__(
        self: dok_array[_ScalarT, _1D],
        /,
        arg1: Sequence[_ScalarT],
        shape: _ToShape1D | None = None,
        dtype: onp.ToDType[_ScalarT] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like (know dtype), dtype: None
    def __init__(
        self: dok_array[_ScalarT, _2D],
        /,
        arg1: Sequence[Sequence[_ScalarT] | onp.CanArrayND[_ScalarT]],  # assumes max. 2-d
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self: dok_array[_ScalarT, _AnyD],
        /,
        arg1: _ToMatrix[_ScalarT],
        shape: _ToShape1D | _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: float64-like | None
    def __init__(
        self: dok_array[np.float64, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like | None
    def __init__(
        self: dok_array[np.float64, _2D],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: bool-like (positional)
    def __init__(
        self: dok_array[np.bool_, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: bool-like (keyword)
    def __init__(
        self: dok_array[np.bool_, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (positional)
    def __init__(
        self: dok_array[np.bool_, _2D],
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
        self: dok_array[np.bool_, _2D],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: int-like (positional)
    def __init__(
        self: dok_array[np.int64, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: int-like (keyword)
    def __init__(
        self: dok_array[np.int64, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (positional)
    def __init__(
        self: dok_array[np.int64, _2D],
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
        self: dok_array[np.int64, _2D],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: dok_array[np.complex128, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d shape-like, dtype: complex128-like (keyword)
    def __init__(
        self: dok_array[np.complex128, _1D],
        /,
        arg1: _ToShape1D,
        shape: _ToShape1D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: dok_array[np.complex128, _2D],
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
        self: dok_array[np.complex128, _2D],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d array-like bool, dtype: bool-like | None
    def __init__(
        self: dok_array[np.bool_, _1D],
        /,
        arg1: onp.ToJustBoolStrict1D,
        shape: _ToShape1D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like bool, dtype: bool-like | None
    def __init__(
        self: dok_array[np.bool_, _2D],
        /,
        arg1: onp.ToJustBoolStrict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d array-like ~int, dtype: int-like | None
    def __init__(
        self: dok_array[np.int64, _1D],
        /,
        arg1: onp.ToJustInt64Strict1D,
        shape: _ToShape1D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~int, dtype: int-like | None
    def __init__(
        self: dok_array[np.int64, _2D],
        /,
        arg1: onp.ToJustInt64Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d array-like ~float, dtype: float64-like | None
    def __init__(
        self: dok_array[np.float64, _1D],
        /,
        arg1: onp.ToJustFloat64Strict1D,
        shape: _ToShape1D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~float, dtype: float64-like | None
    def __init__(
        self: dok_array[np.float64, _2D],
        /,
        arg1: onp.ToJustFloat64Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-d array-like ~complex, dtype: complex128-like | None
    def __init__(
        self: dok_array[np.complex128, _1D],
        /,
        arg1: onp.ToJustComplex128Strict1D,
        shape: _ToShape1D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~complex, dtype: complex128-like | None
    def __init__(
        self: dok_array[np.complex128, _2D],
        /,
        arg1: onp.ToJustComplex128Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-D, dtype: <known> (positional)
    def __init__(
        self: dok_array[_ScalarT, _1D],
        /,
        arg1: onp.ToComplexStrict1D,
        shape: _ToShape1D | None,
        dtype: onp.ToDType[_ScalarT],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 1-D, dtype: <known> (keyword)
    def __init__(
        self: dok_array[_ScalarT, _1D],
        /,
        arg1: onp.ToComplexStrict1D,
        shape: _ToShape1D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (positional)
    def __init__(
        self: dok_array[_ScalarT, _2D],
        /,
        arg1: onp.ToComplexStrict2D,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (keyword)
    def __init__(
        self: dok_array[_ScalarT, _2D],
        /,
        arg1: onp.ToComplexStrict2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <unknown>
    def __init__(
        self: dok_array[Any, _AnyD],
        /,
        arg1: onp.ToComplex1D | onp.ToComplex2D,
        shape: _ToShape1D | _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @overload
    # pyrefly: ignore[bad-override]
    def __getitem__(self, key: onp.CanArrayND[np.bool_ | npc.integer] | list[int] | slice, /) -> Self: ...
    @overload
    def __getitem__(self: dok_array[_ScalarT, _ShapeT], key: _spbase[np.bool_, _ShapeT], /) -> dok_array[_ScalarT, _ShapeT]: ...
    @overload
    def __getitem__(self: dok_array[_ScalarT, _NoD], key: _ToKey1D, /) -> Any: ...
    @overload
    def __getitem__(self: dok_array[_ScalarT, _2D], key: _ToKey2D, /) -> _ScalarT: ...
    @overload
    def __getitem__(self: dok_array[_ScalarT, _1D], key: _ToKey1D, /) -> _ScalarT: ...
    @overload
    def __getitem__(self: dok_array[_ScalarT, _2D], key: _ToKey1D, /) -> coo_array[_ScalarT, _1D]: ...  # pyright: ignore[reportIncompatibleMethodOverride]

    # NOTE: This horrible code duplication is required due to the lack of higher-kinded typing (HKT) support.
    # https://github.com/python/typing/issues/548
    @override
    @overload
    @classmethod
    def fromkeys(  # pyrefly: ignore[bad-override]
        cls: type[dok_array[np.bool_, _NoD]], iterable: _ToKeys, v: onp.ToBool, /
    ) -> dok_array[np.bool_, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.bool_, _2D]], iterable: _ToKeys2, v: onp.ToBool, /) -> dok_array[np.bool_, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.bool_, _1D]], iterable: _ToKeys1, v: onp.ToBool, /) -> dok_array[np.bool_, _1D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[_ScalarT, _NoD]], iterable: _ToKeys, v: _ScalarT, /) -> dok_array[_ScalarT, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[_ScalarT, _2D]], iterable: _ToKeys2, v: _ScalarT, /) -> dok_array[_ScalarT, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[_ScalarT, _1D]], iterable: _ToKeys1, v: _ScalarT, /) -> dok_array[_ScalarT, _1D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.int_, _NoD]], iterable: _ToKeys, v: op.JustInt = 1, /) -> dok_array[np.int_, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.int_, _2D]], iterable: _ToKeys2, v: op.JustInt = 1, /) -> dok_array[np.int_, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.int_, _1D]], iterable: _ToKeys1, v: op.JustInt = 1, /) -> dok_array[np.int_, _1D]: ...
    @overload
    @classmethod
    def fromkeys(
        cls: type[dok_array[np.float64, _NoD]], iterable: _ToKeys, v: op.JustFloat, /
    ) -> dok_array[np.float64, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.float64, _2D]], iterable: _ToKeys2, v: op.JustFloat, /) -> dok_array[np.float64, _2D]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_array[np.float64, _1D]], iterable: _ToKeys1, v: op.JustFloat, /) -> dok_array[np.float64, _1D]: ...
    @overload
    @classmethod
    def fromkeys(
        cls: type[dok_array[np.complex128, _NoD]], iterable: _ToKeys, v: op.JustComplex, /
    ) -> dok_array[np.complex128, _AnyD]: ...
    @overload
    @classmethod
    def fromkeys(
        cls: type[dok_array[np.complex128, _2D]], iterable: _ToKeys2, v: op.JustComplex, /
    ) -> dok_array[np.complex128, _2D]: ...
    @overload
    @classmethod
    def fromkeys(  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]
        cls: type[dok_array[np.complex128, _1D]], iterable: _ToKeys1, v: op.JustComplex, /
    ) -> dok_array[np.complex128, _1D]: ...

#
class dok_matrix(_dok_base[_ScalarT_co, _2D], spmatrix[_ScalarT_co], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_matrix[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_matrix[_ScalarT]: ...
    @type_check_only
    def __assoc_as_float32__(self, /) -> dok_matrix[np.float32]: ...
    @type_check_only
    def __assoc_as_float64__(self, /) -> dok_matrix[np.float64]: ...

    #
    @property
    @override
    def ndim(self, /) -> Literal[2]: ...

    # NOTE: keep the in sync with `dok_array.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self: dok_matrix[_ScalarT],  # this self annotation works around a mypy bug
        /,
        arg1: _ToMatrix[_ScalarT],
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: None
    def __init__(
        self: dok_matrix[np.float64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: dok_matrix[np.bool_],
        /,
        arg1: onp.ToJustBoolStrict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.int, dtype: int-like | None
    def __init__(
        self: dok_matrix[np.int64],
        /,
        arg1: onp.ToJustInt64Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: dok_matrix[np.float64],
        /,
        arg1: onp.ToJustFloat64Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: dok_matrix[np.complex128],
        /,
        arg1: onp.ToJustComplex128Strict2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
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
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @overload
    def __getitem__(  # pyrefly: ignore[bad-override]
        self, key: _ToKey1D | onp.CanArrayND[np.bool_ | npc.integer] | _spbase[np.bool_, _2D] | list[int] | slice, /
    ) -> Self: ...
    @overload
    def __getitem__(self, key: _ToKey2D, /) -> _ScalarT_co: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

    #
    @overload
    @classmethod
    def fromkeys(  # pyrefly: ignore[bad-override]
        cls: type[dok_matrix[np.bool_]], iterable: _ToKeys2, v: onp.ToBool, /
    ) -> dok_matrix[np.bool_]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_matrix[_ScalarT]], iterable: _ToKeys2, v: _ScalarT, /) -> dok_matrix[_ScalarT]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_matrix[np.int_]], iterable: _ToKeys2, v: op.JustInt = 1, /) -> dok_matrix[np.int_]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_matrix[np.float64]], iterable: _ToKeys2, v: op.JustFloat, /) -> dok_matrix[np.float64]: ...
    @overload
    @classmethod
    def fromkeys(cls: type[dok_matrix[np.complex128]], iterable: _ToKeys2, v: op.JustComplex, /) -> dok_matrix[np.complex128]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

#
def isspmatrix_dok(x: object) -> TypeIs[dok_matrix]: ...
