import sys
import types as _types
from collections.abc import Iterator
from typing import (
    Any,
    Self,
    TypeAliasType,
    _SpecialForm,
    final,
    override,
    type_check_only,
)

if sys.version_info >= (3, 13):
    from typing_extensions import ParamSpec, TypeVar, TypeVarTuple
else:
    from typing import ParamSpec, TypeVar, TypeVarTuple

__all__ = "AnnotatedAlias", "GenericType", "LiteralAlias", "UnionAlias"

type _TypeExpr = type | _types.GenericAlias | GenericType | TypeAliasType
type _TypeParam = TypeVar | ParamSpec | TypeVarTuple | UnpackAlias[tuple[object, ...]]

# represents `typing._GenericAlias`
# NOTE: This is different from `typing.GenericAlias`!
class GenericType:
    @property
    def __origin__(self, /) -> _TypeExpr | _SpecialForm: ...
    @property
    def __args__(self, /) -> tuple[Any, ...]: ...
    @property
    def __parameters__(self, /) -> tuple[_TypeParam, ...]: ...
    @override
    def __init_subclass__(cls, /, *, _root: bool = ...) -> None: ...
    def __init__(
        self,
        origin: _TypeExpr | _SpecialForm,
        args: tuple[object, ...] | object,
        /,
    ) -> None: ...
    def __or__(self, rhs: type | object, /) -> UnionAlias: ...
    def __ror__(self, lhs: type | object, /) -> UnionAlias: ...
    def __getitem__(self, args: type | object, /) -> GenericType: ...
    def copy_with(self, params: object, /) -> GenericType: ...
    def __iter__(self, /) -> Iterator[UnpackAlias[Self]]: ...
    def __call__(self, /, *args: object, **kwargs: object) -> _SpecialForm | object: ...
    def __instancecheck__(self, obj: object, /) -> bool: ...
    def __subclasscheck__(self, obj: type, /) -> bool: ...
    def __mro_entries__(self, bases: tuple[type, ...]) -> tuple[type, ...]: ...

@final
class LiteralAlias(GenericType): ...

@final
class UnionAlias(GenericType): ...

@final
class AnnotatedAlias(GenericType):
    @property
    @override
    def __origin__(self, /) -> _TypeExpr: ...
    @property
    def __metadata__(self, /) -> tuple[object, *tuple[object, ...]]: ...

@type_check_only
class UnpackAlias[Ts_co: tuple[object, ...] | TypeVarTuple | GenericType](GenericType):
    @property
    def __typing_unpacked_tuple_args__(self, /) -> tuple[object, ...] | None: ...
    @property
    def __typing_is_unpacked_typevartuple__(self, /) -> bool: ...
    @property
    @override
    def __origin__(self, /) -> type[_SpecialForm]: ...
    @property
    @override
    def __args__(self, /) -> tuple[Ts_co]: ...
    @property
    @override
    def __parameters__(self, /) -> tuple[()] | tuple[TypeVarTuple]: ...
    @override
    def __init__(self, origin: _SpecialForm, args: tuple[Ts_co], /) -> None: ...  # pyrefly:ignore[bad-override]
