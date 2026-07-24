import abc
import functools
from collections.abc import Hashable as _CanHash, Mapping, Sequence as _Sequence
from typing import (
    Any,
    Final,
    Literal as _Literal,
    Protocol,
    TypeAlias,
    overload,
    type_check_only,
)
from typing_extensions import ClassVar, Generic, Never, Self, TypeVar, override

from numba.core.typeconv import castgraph
from numba.core.typing import context, templates
from .npytypes import Array

###

_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_T_contra = TypeVar("_T_contra", contravariant=True)
_T_co = TypeVar("_T_co", covariant=True)
_TypeT_co = TypeVar("_TypeT_co", bound=Type, default=Type, covariant=True)

@type_check_only
class _CanCastPythonValue(Protocol[_T_contra, _T_co]):
    def cast_python_value(self, value: _T_contra, /) -> _T_co: ...

# the step is `Any` instead of `Literal[1] | None` so we can use it as invariant `list`
# element type
_SpecSlice: TypeAlias = slice[None, None, Any]
_Layout: TypeAlias = _Literal["F", "C", "A"]

###

class _TypeMetaclass(abc.ABCMeta):
    def __init__(
        cls, name: str, bases: tuple[type, ...], orig_vars: dict[str, Any]
    ) -> None: ...
    def __call__(cls, *args: Any, **kwargs: Any) -> Any: ...

class Type(metaclass=_TypeMetaclass):
    mutable: ClassVar[bool] = False
    reflected: ClassVar[bool] = False

    name: Final[str]

    def __init__(self, name: str) -> None: ...

    #
    @overload  # () -> Signature
    def __call__(self, /) -> templates.Signature: ...
    @overload  # (Type, ...) -> Signature
    def __call__(self, arg0: Type, /, *args: Type) -> templates.Signature: ...
    @overload  # (_: object) -> self.cast_python_value(_)
    def __call__(self: _CanCastPythonValue[_T1, _T2], arg0: _T1, /) -> _T2: ...

    #
    def __getitem__(
        self, args: _SpecSlice | tuple[_SpecSlice, ...] | list[_SpecSlice]
    ) -> Array: ...

    #
    @property
    def key(self) -> _CanHash: ...
    @property
    def mangling_args(self) -> tuple[str, tuple[object, ...]]: ...
    @property
    def is_internal(self) -> bool: ...

    #
    def can_convert_to(
        self,
        typingctx: context.Context,
        other: Type,
    ) -> castgraph.Conversion | None: ...
    def can_convert_from(
        self,
        typingctx: context.Context,
        other: Type,
    ) -> castgraph.Conversion | None: ...
    def cast_python_value(self, args: Never) -> Any: ...

    #
    def is_precise(self) -> bool: ...
    def augment(self, other: Type) -> Self | None: ...
    def unify(self, typingctx: context.Context, other: Type) -> Type | None: ...
    def dump(self, tab: str = "") -> None: ...

class Dummy(Type): ...
class Hashable(Type): ...

class Number(Hashable):
    @overload
    def unify(self, typingctx: context.Context, other: Number) -> Number: ...
    @overload
    def unify(self, typingctx: context.Context, other: Type) -> Number | None: ...

class Callable(Type):
    @abc.abstractmethod
    def get_call_type(
        self, context: context.Context, args: tuple[Type, ...], kws: Mapping[str, Type]
    ) -> templates.Signature | None: ...
    @abc.abstractmethod
    def get_call_signatures(self) -> tuple[_Sequence[templates.Signature], bool]: ...
    @abc.abstractmethod
    def get_impl_key(self, sig: templates.Signature) -> _CanHash: ...

class DTypeSpec(Type):
    @property
    @abc.abstractmethod
    def dtype(self) -> Any: ...

class IterableType(Type):
    @property
    @abc.abstractmethod
    def iterator_type(self) -> Type: ...

class Sized(Type): ...

class ConstSized(Sized):
    @abc.abstractmethod
    def __len__(self) -> int: ...

class IteratorType(IterableType):
    def __init__(self, name: str, **kwargs: Never) -> None: ...
    @property
    @abc.abstractmethod
    def yield_type(self) -> Type: ...
    @property
    def iterator_type(self) -> Self: ...

class Container(Sized, IterableType): ...
class Sequence(Container): ...

class MutableSequence(Sequence):
    mutable: ClassVar[bool] = True

class ArrayCompatible(Type):
    array_priority: ClassVar[float] = 0.0

    @property
    @abc.abstractmethod
    def as_array(self) -> ArrayCompatible: ...

    #
    @functools.cached_property
    def ndim(self) -> int: ...
    @functools.cached_property
    def layout(self) -> _Layout: ...
    @functools.cached_property
    def dtype(self) -> Type: ...

class Literal(Type, Generic[_T_co]):
    ctor_map: ClassVar[dict[type, type[Literal[Any]]]] = ...

    def __init__(self, value: _T_co) -> None: ...
    @property
    def literal_value(self) -> _T_co: ...
    @property
    def literal_type(self) -> Type: ...

class TypeRef(Dummy, Generic[_TypeT_co]):
    instance_type: _TypeT_co  # readonly

    def __init__(self, instance_type: _TypeT_co) -> None: ...
    @property
    @override
    def key(self) -> _TypeT_co: ...

class InitialValue(Generic[_T_co]):
    def __init__(self, initial_value: _T_co) -> None: ...
    @property
    def initial_value(self) -> _T_co: ...

class Poison(Type, Generic[_TypeT_co]):  # like `Never`
    def __init__(self, ty: _TypeT_co) -> None: ...
    def __unliteral__(self) -> Poison[Self]: ...
    @override
    def unify(self, typingctx: context.Context, other: Type) -> None: ...
