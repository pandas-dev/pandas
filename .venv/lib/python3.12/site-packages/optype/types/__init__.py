import sys
from typing import Literal, LiteralString, Protocol

if sys.version_info >= (3, 13):
    from typing import runtime_checkable
else:
    from typing_extensions import runtime_checkable

from ._typeforms import AnnotatedAlias, GenericType, LiteralAlias, UnionAlias

__all__ = (
    "AnnotatedAlias",
    "Deprecated",
    "GenericType",
    "LiteralAlias",
    "ProtocolType",
    "RuntimeProtocolType",
    "UnionAlias",
    "WrappedFinalType",
)


def __dir__() -> tuple[str, ...]:
    return __all__


@runtime_checkable
class Deprecated(Protocol):
    """
    A runtime-protocol that represents a type or method that's decorated with
    `@typing.deprecated` or `@typing_extensions.deprecated`.
    """

    __deprecated__: str


@runtime_checkable
class WrappedFinalType(Protocol):
    """
    A runtime-protocol that represents a type or method that's decorated with
    `@typing.final` or `@typing_extensions.final`.

    Note that the name `HasFinal` isn't used, because `__final__` is
    undocumented, and therefore not a part of the public API.
    """

    __final__: Literal[True]


@runtime_checkable
class ProtocolType(Protocol):
    _is_protocol: Literal[True]

    if sys.version_info >= (3, 12, 0):
        __protocol_attrs__: set[LiteralString]


@runtime_checkable
class RuntimeProtocolType(ProtocolType, Protocol):
    _is_runtime_protocol: Literal[True]

    if sys.version_info >= (3, 12, 2):
        __non_callable_proto_members__: set[LiteralString]
