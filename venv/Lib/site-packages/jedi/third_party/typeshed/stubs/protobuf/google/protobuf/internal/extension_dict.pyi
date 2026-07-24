from collections.abc import Iterator
from typing import Any, Generic, TypeVar, type_check_only

from google.protobuf.descriptor import FieldDescriptor
from google.protobuf.internal.containers import RepeatedCompositeFieldContainer, RepeatedScalarFieldContainer
from google.protobuf.message import Message

_ContainerMessageT = TypeVar("_ContainerMessageT", bound=Message)
_ExtenderMessageT = TypeVar(
    "_ExtenderMessageT",
    bound=Message | RepeatedScalarFieldContainer[Any] | RepeatedCompositeFieldContainer[Any] | bool | float | str | bytes,
)

@type_check_only
class _ExtensionFieldDescriptor(FieldDescriptor, Generic[_ContainerMessageT, _ExtenderMessageT]): ...

class _ExtensionDict(Generic[_ContainerMessageT]):
    def __init__(self, extended_message: _ContainerMessageT) -> None: ...
    def __getitem__(
        self, extension_handle: _ExtensionFieldDescriptor[_ContainerMessageT, _ExtenderMessageT]
    ) -> _ExtenderMessageT: ...
    def __len__(self) -> int: ...
    def __setitem__(
        self, extension_handle: _ExtensionFieldDescriptor[_ContainerMessageT, _ExtenderMessageT], value: _ExtenderMessageT
    ) -> None: ...
    def __delitem__(self, extension_handle: _ExtensionFieldDescriptor[_ContainerMessageT, _ExtenderMessageT]) -> None: ...
    def __iter__(self) -> Iterator[_ExtensionFieldDescriptor[_ContainerMessageT, Any]]: ...
    def __contains__(self, extension_handle: _ExtensionFieldDescriptor[_ContainerMessageT, _ExtenderMessageT]) -> bool: ...
