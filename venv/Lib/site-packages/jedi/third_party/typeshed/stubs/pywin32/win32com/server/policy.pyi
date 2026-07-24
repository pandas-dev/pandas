from _typeshed import Incomplete
from abc import ABC, abstractmethod
from typing import Any, Final

import _win32typing

__author__: Final[str]
S_OK: Final = 0
IDispatchType: Incomplete
IUnknownType: Incomplete
regSpec: str
regPolicy: str
regDispatcher: str
regAddnPath: str

def CreateInstance(clsid, reqIID: _win32typing.PyIID) -> _win32typing.PyIUnknown: ...

class BasicWrapPolicy(ABC):
    def __init__(self, object) -> None: ...
    def _InvokeEx_(self, dispid, lcid, wFlags, args, kwargs, serviceProvider) -> tuple[Incomplete]: ...
    @abstractmethod
    def _invokeex_(self, dispid, lcid, wFlags, args, kwargs, serviceProvider) -> tuple[Incomplete]: ...

class MappedWrapPolicy(BasicWrapPolicy):
    _dispid_to_func_: dict[int, str]
    def _invokeex_(self, dispid, lcid, wFlags, args, kwargs, serviceProvider) -> tuple[Incomplete]: ...

class DesignatedWrapPolicy(MappedWrapPolicy): ...
class EventHandlerPolicy(DesignatedWrapPolicy): ...

class DynamicPolicy(BasicWrapPolicy):
    def _invokeex_(self, dispid, lcid, wFlags, args, kwargs, serviceProvider) -> tuple[Incomplete]: ...

DefaultPolicy = DesignatedWrapPolicy

# Imports an arbitrary object by it's fully-qualified name.
def resolve_func(spec: str) -> Any: ...

# Imports and calls an arbitrary callable by it's fully-qualified name.
def call_func(spec: str, *args: Any) -> Any: ...

DISPATCH_METHOD: int
DISPATCH_PROPERTYGET: int
DISPATCH_PROPERTYPUT: int
DISPATCH_PROPERTYPUTREF: int
DISPID_EVALUATE: int
DISPID_NEWENUM: int
DISPID_PROPERTYPUT: int
DISPID_STARTENUM: int
DISPID_VALUE: int
