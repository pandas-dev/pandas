from collections.abc import Callable
from typing import Any, TypeAlias, type_check_only

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
from numba.core.compiler import CompileResult
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
_GUFuncT = TypeVar("_GUFuncT", bound=GUFunc)
_GUFuncT_co = TypeVar("_GUFuncT_co", bound=GUFunc, covariant=True)

_ToSignature: TypeAlias = Signature | tuple[Type | Signature, ...] | str

@type_check_only
class GUFuncKernel(npyimpl._Kernel, Generic[_GUFuncT_co]):
    __name__: str

    dufunc: _GUFuncT_co  # not a typo

    inner_sig: Signature | None
    cres: CompileResult | None

###

def make_gufunc_kernel(_dufunc: _GUFuncT) -> type[GUFuncKernel[_GUFuncT]]: ...

class GUFuncLowerer(
    UfuncLowererBase[_GUFuncT_co, GUFuncKernel],
    Generic[_GUFuncT_co],
):
    @override
    def __init__(self, gufunc: _GUFuncT_co) -> None: ...

class GUFunc(
    serialize.ReduceMixin,
    UfuncBase[
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
        str,
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

    def __init__(
        self,
        py_func: Callable[..., object],
        signature: str,
        identity: _IdentityT_co | None = None,
        cache: bool | None = None,
        is_dynamic: bool = False,
        targetoptions: dict[str, Any] | None = None,
        writable_args: tuple[int | str, ...] = ...,
    ) -> None: ...
    def __call__(
        self: GUFunc[Callable[_Tss, _T]], *args: _Tss.args, **kws: _Tss.kwargs
    ) -> _T: ...

    #
    def add(self, fty: _ToSignature) -> None: ...
    def expected_ndims(self) -> tuple[tuple[int, ...], tuple[int, ...]]: ...
    def build_ufunc(self) -> Self: ...
    def match_signature(
        self, ewise_types: tuple[Type, ...], sig: Signature
    ) -> bool: ...

    #
    @property
    def is_dynamic(self) -> bool: ...

    # `ReduceMixin` methods
    @override
    def _reduce_states(self) -> dict[str, Any]: ...
    @classmethod
    @override
    def _rebuild(  # type: ignore[override]
        cls,
        py_func: Callable[..., object],
        signature: str,
        identity: _T,
        cache: bool,
        is_dynamic: bool,
        targetoptions: dict[str, Any] | None,
        writable_args: tuple[int | str, ...],
        typesigs: list[Signature],
        frozen: bool,
    ) -> GUFunc[
        _CallT_co,
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _T,
    ]: ...
