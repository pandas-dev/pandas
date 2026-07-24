from collections.abc import Callable
from typing import Any, TypeAlias, overload, type_check_only

from _typeshed import Incomplete
from llvmlite.ir import IRBuilder
from typing_extensions import (
    Concatenate,
    Generic,
    Never,
    ParamSpec,
    Self,
    TypeVar,
    override,
)

from numba.core import serialize
from numba.core.base import BaseContext
from numba.core.compiler import CompileResult
from numba.core.dispatcher import Dispatcher
from numba.core.types import Type
from numba.core.typing.templates import Signature
from numba.np import npyimpl

from .ufunc_base import (
    UfuncBase,
    UfuncLowererBase,
    _AccumulateT_co,
    _IdentityT_co,
    _OuterT_co,
    _ReduceAtT_co,
    _ReduceT_co,
)

###

_T = TypeVar("_T")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_T3 = TypeVar("_T3")
_Tss = ParamSpec("_Tss")
_CallT_co = TypeVar(
    "_CallT_co",
    bound=Callable[Concatenate[Never, ...], object],
    default=Callable[Concatenate[Any, ...], Any],
    covariant=True,
)
_AtT_co = TypeVar(
    "_AtT_co",
    bound=Callable[[Never, Never, Never], None],
    default=Callable[[Any, Any, Any | None], None],
    covariant=True,
)
_DUFuncT = TypeVar("_DUFuncT", bound=DUFunc)
_DUFuncT_co = TypeVar("_DUFuncT_co", bound=DUFunc, covariant=True)

_ToSignature: TypeAlias = Signature | tuple[Type | Signature, ...] | str

@type_check_only
class DUFuncKernel(npyimpl._Kernel, Generic[_DUFuncT_co]):
    __name__: str

    dufunc: _DUFuncT_co

    inner_sig: Signature | None
    cres: CompileResult | None

###

# undocumented
class UfuncAtIterator:
    ufunc: DUFunc
    a: Incomplete
    a_ty: Type
    indices: Incomplete
    indices_ty: Type
    b: Incomplete
    b_ty: Type

    def __init__(
        self,
        ufunc: DUFunc,
        a: Incomplete,
        a_ty: Type,
        indices: Incomplete,
        indices_ty: Type,
        b: Incomplete | None = None,
        b_ty: Type | None = None,
    ) -> None: ...
    def run(self, context: BaseContext, builder: IRBuilder) -> None: ...
    def need_advanced_indexing(self) -> None: ...

def make_dufunc_kernel(_dufunc: _DUFuncT) -> type[DUFuncKernel[_DUFuncT]]: ...

class DUFuncLowerer(
    UfuncLowererBase[_DUFuncT_co, DUFuncKernel],
    Generic[_DUFuncT_co],
):
    @override
    def __init__(self, dufunc: _DUFuncT_co) -> None: ...

#
class DUFunc(
    serialize.ReduceMixin,
    UfuncBase[
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
        None,
    ],
    Generic[
        _CallT_co,
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
    ],
):
    __name__: str
    __doc__: str | None
    reorderable: bool

    def __init__(
        self,
        py_func: Callable[..., object],
        identity: _IdentityT_co | None = None,
        cache: bool = False,
        targetoptions: dict[str, Any] | None = None,
    ) -> None: ...
    def __call__(
        self: DUFunc[Callable[_Tss, _T]], *args: _Tss.args, **kws: _Tss.kwargs
    ) -> _T: ...

    #
    def build_ufunc(self) -> Self: ...
    def add(self, sig: _ToSignature) -> CompileResult: ...
    def match_signature(
        self, ewise_types: tuple[Type, ...], sig: Signature
    ) -> bool: ...

    #
    @property
    def targetoptions(self) -> dict[str, Any]: ...

    #
    @override  # type: ignore[override]
    @overload
    def at(
        self: DUFunc[Any, Callable[[_T1, _T2, None], Any]],
        a: _T1,
        indices: _T2,
        b: None = None,
    ) -> None: ...
    @overload
    def at(
        self: DUFunc[Any, Callable[[_T1, _T2, _T3], Any]],
        a: _T1,
        indices: _T2,
        b: _T3,
    ) -> None: ...

    # `ReduceMixin` methods
    @override
    def _reduce_states(self) -> dict[str, Any]: ...
    @classmethod
    @override
    def _rebuild(  # type: ignore[override]
        cls,
        dispatcher: Dispatcher,
        identity: _T,
        frozen: bool,
        siglist: list[Signature],
    ) -> DUFunc[
        _CallT_co,
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _T,
    ]: ...
