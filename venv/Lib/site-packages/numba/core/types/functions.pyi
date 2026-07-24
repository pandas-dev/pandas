import weakref
from collections.abc import Callable as _Callable
from collections.abc import Mapping, Sequence
from typing import Any, Final, Hashable, Iterable, NamedTuple, Self, TypeAlias

from _typeshed import Incomplete
from typing_extensions import Generic, TypeVar, override

from numba.core import dispatcher
from numba.core.typing.context import Context
from numba.core.typing.templates import FunctionTemplate, Signature

from .abstract import Callable, DTypeSpec, Dummy, Literal, Type
from .common import Opaque

_T_co = TypeVar("_T_co", covariant=True)
_ThisT = TypeVar("_ThisT", bound=Type)
_ThisT_co = TypeVar("_ThisT_co", bound=Type, default=Type, covariant=True)
_TemplateTypeT = TypeVar("_TemplateTypeT", bound=_TemplateType)
_DispatcherT_co = TypeVar(
    "_DispatcherT_co",
    bound=dispatcher.Dispatcher,
    default=dispatcher.Dispatcher,
    covariant=True,
)
_TupleT_co = TypeVar(
    "_TupleT_co",
    bound=tuple[Any, ...],
    default=tuple[Any, ...],
    covariant=True,
)

_TemplateType: TypeAlias = type[FunctionTemplate]
_AnyCallable: TypeAlias = _Callable[..., Any]
_BoundFunctionKey: TypeAlias = tuple[Hashable, _ThisT, _Callable | None]
_GetPointerFn: TypeAlias = _Callable[[Any], int]

###

def argsnkwargs_to_str(args: Iterable[object], kwargs: Mapping[str, object]) -> str: ...

class BaseFunction(Callable):
    templates: Final[tuple[_TemplateType, ...]]
    typing_key: Final[Hashable]

    @override  # this free `_TemplateT` is needed because `list` is invariant
    def __init__(
        self,
        template: tuple[_TemplateTypeT, ...] | list[_TemplateTypeT] | _TemplateTypeT,
    ) -> None: ...

    #
    @property
    @override
    def key(self) -> tuple[Incomplete, tuple[_TemplateType, ...]]: ...

    #
    @override
    def get_impl_key(self, sig: Signature) -> Hashable: ...
    @override
    def get_call_type(
        self,
        context: Context,
        args: tuple[Type, ...],
        kws: Mapping[str, Type],
    ) -> Signature | None: ...
    @override
    def get_call_signatures(self) -> tuple[Sequence[Signature], bool]: ...

class Function(BaseFunction, Opaque): ...

class BoundFunction(Callable, Opaque, Generic[_ThisT_co]):
    template: Final[_TemplateType]
    typing_key: Final[Hashable]
    this: _ThisT_co

    @override
    def __init__(self, template: _TemplateType, this: _ThisT_co) -> None: ...

    #
    @override
    def unify(
        self,
        typingctx: Context,
        other: _ThisT,
    ) -> BoundFunction[_ThisT_co | _ThisT] | None: ...

    #
    def copy(self, this: _ThisT) -> BoundFunction[_ThisT]: ...

    #
    @property
    @override
    def key(self) -> _BoundFunctionKey[_ThisT_co]: ...
    @override
    def get_impl_key(self, sig: Signature) -> Hashable: ...
    @override
    def get_call_type(
        self,
        context: Context,
        args: tuple[Type, ...],
        kws: Mapping[str, Type],
    ) -> Signature | None: ...
    @override
    def get_call_signatures(self) -> tuple[Sequence[Signature], bool]: ...

class MakeFunctionLiteral(Literal, Opaque): ...

class WeakType(Type, Generic[_T_co]):
    @property
    @override
    def key(self) -> weakref.ref[_T_co]: ...

    #
    @override
    def __eq__(self, other: Self) -> bool: ...  # type: ignore[override]
    @override
    def __hash__(self) -> int: ...

class Dispatcher(WeakType[_DispatcherT_co], Callable, Dummy, Generic[_DispatcherT_co]):
    @override
    def __init__(self, dispatcher: _DispatcherT_co) -> None: ...

    #
    @property
    def dispatcher(self) -> _DispatcherT_co: ...

    #
    @override
    def get_call_type(
        self,
        context: Context,
        args: tuple[Type, ...],
        kws: Mapping[str, Type],
    ) -> Signature | None: ...
    @override
    def get_call_signatures(self) -> tuple[Sequence[Signature], bool]: ...
    @override
    def get_impl_key(self, sig: Signature) -> _AnyCallable: ...
    def get_overload(self, sig: Signature) -> _AnyCallable: ...

class ObjModeDispatcher(Dispatcher[_DispatcherT_co], Generic[_DispatcherT_co]): ...

class ExternalFunctionPointer(BaseFunction):
    sig: Final[Signature]
    requires_gil: Final[bool]
    get_pointer: Final[_GetPointerFn]
    cconv: Final[str | None]

    @override
    def __init__(
        self,
        sig: Signature,
        get_pointer: _GetPointerFn,
        cconv: str | None = None,
    ) -> None: ...

    #
    @property
    @override
    def key(self) -> tuple[Signature, str | None, _GetPointerFn]: ...  # type: ignore[override]

class ExternalFunction(Function):
    symbol: Final[str]
    sig: Final[Signature]

    @override
    def __init__(self, symbol: str, sig: Signature) -> None: ...

    #
    @property
    @override
    def key(self) -> tuple[str, Signature]: ...  # type: ignore[override]

class NamedTupleClass(Callable, Opaque, Generic[_TupleT_co]):
    instance_class: type[_TupleT_co]

    @override
    def __init__(self, instance_class: type[_TupleT_co]) -> None: ...

    #
    @property
    @override
    def key(self) -> type[_TupleT_co]: ...

    #
    @override
    def get_call_type(
        self,
        context: Context,
        args: tuple[Type, ...],
        kws: Mapping[str, Type],
    ) -> None: ...
    @override
    def get_call_signatures(self) -> tuple[tuple[()], bool]: ...
    @override
    def get_impl_key(self, sig: Signature) -> type[Self]: ...

class NumberClass(Callable, DTypeSpec, Opaque, Generic[_T_co]):
    instance_type: _T_co

    @override
    def __init__(self, instance_type: _T_co) -> None: ...

    #
    @property
    @override
    def key(self) -> _T_co: ...
    @property
    @override
    def dtype(self) -> _T_co: ...

    #
    @override
    def get_call_type(
        self,
        context: Context,
        args: tuple[Type, ...],
        kws: Mapping[str, Type],
    ) -> None: ...
    @override
    def get_call_signatures(self) -> tuple[tuple[()], bool]: ...
    @override
    def get_impl_key(self, sig: Signature) -> type[Self]: ...

class _RecursiveCallOverloads(NamedTuple):
    qualname: str
    uid: int

class RecursiveCall(Opaque, Generic[_DispatcherT_co]):
    dispatcher_type: Dispatcher[_DispatcherT_co]

    @override
    def __init__(self, dispatcher_type: Dispatcher[_DispatcherT_co]) -> None: ...
    @property
    @override
    def key(self) -> Dispatcher[_DispatcherT_co]: ...

    #
    def add_overloads(
        self,
        args: tuple[Type, ...],
        qualname: str,
        uid: int,
    ) -> None: ...
    def get_overloads(self, args: tuple[Type, ...]) -> _RecursiveCallOverloads: ...
