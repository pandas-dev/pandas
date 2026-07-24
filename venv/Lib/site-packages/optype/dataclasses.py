"""
Runtime-protocols for the `dataclasses` standard library.
https://docs.python.org/3/library/dataclasses.html
"""

import dataclasses
import sys
from typing import Any, ClassVar, Protocol, override

if sys.version_info >= (3, 13):
    from typing import runtime_checkable
else:
    from typing_extensions import runtime_checkable

__all__ = ("HasDataclassFields",)


def __dir__() -> tuple[str]:
    return __all__


###


@runtime_checkable
class HasDataclassFields(Protocol):
    """Can be used to check whether a type or instance is a dataclass."""

    __dataclass_fields__: ClassVar[dict[str, dataclasses.Field[Any]]]

    # Because of https://github.com/python/mypy/issues/3939 just having
    # `__dataclass_fields__` is insufficient for `issubclass` checks.
    @override
    @classmethod
    def __subclasshook__(cls, c: type, /) -> bool:
        """Customize the subclass check."""
        return hasattr(c, "__dataclass_fields__")
