"""
Runtime-protocols for the `dataclasses` standard library.
https://docs.python.org/3/library/dataclasses.html
"""

import dataclasses
import sys
from collections.abc import Mapping
from typing import Any, Protocol, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeVar, runtime_checkable
else:
    from typing_extensions import TypeVar, runtime_checkable

__all__ = ("HasDataclassFields",)


def __dir__() -> tuple[str]:
    return __all__


###

__Field: TypeAlias = dataclasses.Field[Any]
_FieldsT = TypeVar("_FieldsT", bound=Mapping[str, __Field], default=dict[str, __Field])


###


@runtime_checkable
class HasDataclassFields(Protocol[_FieldsT]):
    """Can be used to check whether a type or instance is a dataclass."""

    __dataclass_fields__: _FieldsT
