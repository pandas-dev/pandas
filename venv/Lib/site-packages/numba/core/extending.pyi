import collections
import inspect
import typing
import weakref
from collections.abc import Iterable, Mapping, Sequence
from typing import (
    Callable,
    ClassVar,
    Final,
    Generic,
    Literal,
    ParamSpec,
    Protocol,
    TypeAlias,
    final,
    type_check_only,
)

from _typeshed import Unused
from typing_extensions import NamedTuple, Never, Self, TypeIs, TypeVar, override

from numba._helperlib import _import_cython_function as _import_cython_function

from .datamodel import models as models
from .datamodel import register_default
from .decorators import _FnInline, _JITOptions
from .dispatcher import Dispatcher
from .imputils import lower_builtin as lower_builtin
from .imputils import lower_cast as lower_cast
from .imputils import lower_getattr as lower_getattr
from .imputils import lower_getattr_generic as lower_getattr_generic
from .imputils import lower_setattr as lower_setattr
from .imputils import lower_setattr_generic as lower_setattr_generic
from .pythonapi import NativeValue as NativeValue
from .pythonapi import box as box
from .pythonapi import reflect as reflect
from .pythonapi import unbox as unbox
from .serialize import ReduceMixin as ReduceMixin
from .types import Type
from .typing.asnumbatype import as_numba_type as as_numba_type
from .typing.context import BaseContext
from .typing.templates import infer as infer
from .typing.templates import infer_getattr as infer_getattr
from .typing.typeof import typeof_impl as typeof_impl

register_model = register_default

###
# stubs-only helpers

_T = TypeVar("_T")
_Pss = ParamSpec("_Pss")
_ContextT = TypeVar("_ContextT", bound=BaseContext)
_CallableT = TypeVar("_CallableT", bound=Callable[..., object])
_CallableT_co = TypeVar(
    "_CallableT_co",
    bound=Callable[..., object],
    default=Callable[..., typing.Any],
    covariant=True,
)

_Target: TypeAlias = Literal["generic", "cpu", "cuda", "npyufunc"] | str
_Inline: TypeAlias = Literal["never", "always"] | _FnInline

@type_check_only
class _TypeCallableDecorator(Protocol):
    def __call__(
        self,
        /,
        # We can't use a single typevar for the entire typing_func callable because
        # the `BaseContext` input type is in a contravariant position.
        typing_func: Callable[[_ContextT], _CallableT],
    ) -> Callable[[_ContextT], _CallableT]: ...

@type_check_only
class _OverloadDecorator(Protocol):
    def __call__(self, /, overload_func: _CallableT) -> _CallableT: ...

@type_check_only
class _RegisterJitableDecorator(Protocol):
    def __call__(self, /, fn: _CallableT) -> _CallableT: ...

@type_check_only
class _IntrinsicDecorator(Protocol):
    def __call__(self, /, func: _CallableT) -> _Intrinsic[_CallableT]: ...

###
# stubs

def type_callable(func: Callable[..., object] | str) -> _TypeCallableDecorator: ...
def overload(
    func: Callable[..., object],
    jit_options: _JITOptions = ...,
    strict: bool = True,
    inline: _Inline = "never",
    prefer_literal: bool = False,
    *,
    target: _Target = ...,
) -> _OverloadDecorator: ...

#
@typing.overload
def register_jitable(fn: _CallableT, /) -> _CallableT: ...
@typing.overload
# we can't use `**kwargs: Unpack[_JitOptions]` because mypy doesn't support PEP 728 yet
def register_jitable(**kwargs: object) -> _RegisterJitableDecorator: ...

#
def overload_attribute(
    typ: Type | type[Type],
    attr: str,
    *,
    inline: _Inline = "never",
    prefer_literal: bool = False,
    base: type = ...,
    target: _Target = ...,
) -> _OverloadDecorator: ...
def overload_method(
    typ: Type | type[Type],
    attr: str,
    *,
    inline: _Inline = "never",
    prefer_literal: bool = False,
    target: _Target = ...,
) -> _OverloadDecorator: ...
def overload_classmethod(
    typ: Type | type[Type],
    attr: str,
    *,
    inline: _Inline = "never",
    prefer_literal: bool = False,
    target: _Target = ...,
) -> _OverloadDecorator: ...
def make_attribute_wrapper(
    typeclass: type[Type], struct_attr: str, python_attr: str
) -> None: ...

#
@final
class _Intrinsic(ReduceMixin, Generic[_CallableT_co]):
    _memo: ClassVar[weakref.WeakValueDictionary[str, _Intrinsic]] = ...
    _recent: ClassVar[collections.deque[_Intrinsic]] = ...

    _name: Final[str]
    _defn: _CallableT_co
    _prefer_literal: Final[bool]
    _ctor_kwargs: Final[dict[str, typing.Any]]

    def __init__(
        self,
        name: str,
        defn: _CallableT_co,
        prefer_literal: bool = False,
        **kwargs: object,
    ) -> None: ...

    #
    @property
    def _uuid(self) -> str: ...
    def _set_uuid(self, u: str) -> None: ...
    def _register(self) -> None: ...

    # always raises `NotImplementedError` on call
    def __call__(
        self: _Intrinsic[Callable[_Pss, object]],
        *args: _Pss.args,
        **kwargs: _Pss.kwargs,
    ) -> Never: ...

    #
    def __deepcopy__(self, memo: Unused) -> Self: ...
    def _reduce_states(self) -> dict[str, typing.Any]: ...
    @classmethod
    @override
    def _rebuild(  # type:ignore[override]
        cls, uuid: str, name: str, defn: _CallableT
    ) -> _Intrinsic[_CallableT]: ...

#
@typing.overload
def intrinsic(func: _CallableT, /) -> _Intrinsic[_CallableT]: ...
@typing.overload
def intrinsic(
    *, prefer_literal: bool = False, **kwargs: object
) -> _IntrinsicDecorator: ...

#
def get_cython_function_address(module_name: str, function_name: str) -> int: ...
def include_path() -> str: ...
def sentry_literal_args(
    pysig: inspect.Signature,
    literal_args: Sequence[str],
    args: Iterable[object],
    kwargs: Mapping[str, object],
) -> None: ...

class SentryLiteralArgs(NamedTuple):
    literal_args: Sequence[str]

    def for_function(self, func: Callable[..., object]) -> BoundLiteralArgs: ...
    def for_pysig(self, pysig: inspect.Signature) -> BoundLiteralArgs: ...

class BoundLiteralArgs(NamedTuple):
    pysig: inspect.Signature
    literal_args: Sequence[str]

    def bind(self, *args: object, **kwargs: object) -> None: ...

def is_jitted(
    function: Callable[_Pss, _T],
) -> TypeIs[Dispatcher[Callable[_Pss, _T]]]: ...
