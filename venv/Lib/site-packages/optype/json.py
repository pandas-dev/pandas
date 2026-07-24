"""
Type aliases for `json` standard library.
This assumes that the default encoder and decoder are used.
"""

import sys
from types import MappingProxyType
from typing import TypeAliasType

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

__all__ = "AnyArray", "AnyObject", "AnyValue", "Array", "Object", "_Value"


def __dir__() -> tuple[str, ...]:
    return __all__


###


type _Primitive = bool | int | float | str | None
type _Value = _Primitive | dict[str, _Value] | list[_Value]
type _AnyValue = (
    _Primitive
    # NOTE: `TypedDict` can't be included here, since it's not a sub*type* of
    # `dict[str, Any]` according to the typing docs and typeshed, even though
    # it **literally** is a subclass of `dict`...
    | dict[str, _AnyValue]
    | MappingProxyType[str, _AnyValue]
    | list[_AnyValue]
    | tuple[_AnyValue, ...]
)
_VT = TypeVar("_VT", bound=_Value, default=_Value)
_AVT = TypeVar("_AVT", bound=_AnyValue, default=_AnyValue)


# Return types of `json.load[s]`

Array = TypeAliasType("Array", list[_VT], type_params=(_VT,))  # noqa: UP040
Object = TypeAliasType("Object", dict[str, _VT], type_params=(_VT,))  # noqa: UP040
# ensure that `Value | Array | Object` is equivalent to `Value`
type Value = _Value | Array | Object


# Input types of `json.dumps`

AnyArray = TypeAliasType("AnyArray", list[_AVT] | tuple[_AVT, ...], type_params=(_AVT,))  # noqa: UP040
AnyObject = TypeAliasType("AnyObject", dict[str, _AVT], type_params=(_AVT,))  # noqa: UP040
type AnyValue = _AnyValue | AnyArray | AnyObject | Value
