"""
Type aliases for `json` standard library.
This assumes that the default encoder and decoder are used.
"""

import sys
from types import MappingProxyType
from typing import TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

__all__ = "AnyArray", "AnyObject", "AnyValue", "Array", "Object", "_Value"


def __dir__() -> tuple[str, ...]:
    return __all__


###


_Primitive: TypeAlias = bool | int | float | str | None
_Value: TypeAlias = _Primitive | dict[str, "_Value"] | list["_Value"]
_AnyValue: TypeAlias = (
    _Primitive
    # NOTE: `TypedDict` can't be included here, since it's not a sub*type* of
    # `dict[str, Any]` according to the typing docs and typeshed, even though
    # it **literally** is a subclass of `dict`...
    | dict[str, "_AnyValue"]
    | MappingProxyType[str, "_AnyValue"]
    | list["_AnyValue"]
    | tuple["_AnyValue", ...]
)

_VT = TypeVar("_VT", bound=_Value, default=_Value)
_AVT = TypeVar("_AVT", bound=_AnyValue, default=_AnyValue)


# Return types of `json.load[s]`

Array: TypeAlias = list[_VT]
Object: TypeAlias = dict[str, _VT]
# ensure that `Value | Array | Object` is equivalent to `Value`
Value: TypeAlias = _Value | Array | Object


# Input types of `json.dumps`

AnyArray: TypeAlias = list[_AVT] | tuple[_AVT, ...]
AnyObject: TypeAlias = dict[str, _AVT]
AnyValue: TypeAlias = _AnyValue | AnyArray | AnyObject | Value
