from typing import Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

_T = TypeVar("_T")
_GetterReturnType_co = TypeVar("_GetterReturnType_co", covariant=True)
_SetterValueType_contra = TypeVar("_SetterValueType_contra", contravariant=True)

@type_check_only
class AsymmetricProperty(Protocol[_GetterReturnType_co, _SetterValueType_contra]):
    @overload
    def __get__(self, obj: None, type: type[object] | None = ..., /) -> property: ...
    @overload
    def __get__(self, obj: object, type: type[object] | None = ..., /) -> _GetterReturnType_co: ...
    def __set__(self, obj: object, value: _SetterValueType_contra, /) -> None: ...

@type_check_only
class AsymmetricPropertyWithDelete(
    AsymmetricProperty[_GetterReturnType_co, _SetterValueType_contra], Protocol[_GetterReturnType_co, _SetterValueType_contra]
):
    def __delete__(self, obj: object, /) -> None: ...

SymmetricProperty: TypeAlias = AsymmetricProperty[_T, _T]
SymmetricPropertyWithDelete: TypeAlias = AsymmetricPropertyWithDelete[_T, _T]
