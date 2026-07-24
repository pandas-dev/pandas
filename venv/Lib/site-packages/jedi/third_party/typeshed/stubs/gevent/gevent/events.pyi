from collections.abc import Callable, Mapping, Sequence
from types import ModuleType
from typing import Any, Protocol, TypeVar, type_check_only
from typing_extensions import TypeAlias

from gevent.hub import Hub
from greenlet import greenlet as greenlet_t
from psutil._ntuples import pmem

_T = TypeVar("_T")
# FIXME: While it would be nice to import Interface from zope.interface here so the
#        mypy plugin will work correctly for the people that use it, it causes all
#        sorts of issues to reference a module that is not stubbed in typeshed, so
#        for now we punt and just define an alias for Interface and implementer we
#        can get rid of later
Interface: TypeAlias = Any

def implementer(interface: Interface, /) -> Callable[[_T], _T]: ...

subscribers: list[Callable[[Any], object]]

@type_check_only
class _PeriodicMonitorThread(Protocol):
    def add_monitoring_function(self, function: Callable[[Hub], object], period: float | None) -> object: ...

class IPeriodicMonitorThread(Interface):
    def add_monitoring_function(function: Callable[[Hub], object], period: float | None) -> object: ...

class IPeriodicMonitorThreadStartedEvent(Interface):
    monitor: IPeriodicMonitorThread

@implementer(IPeriodicMonitorThread)
class PeriodicMonitorThreadStartedEvent:
    ENTRY_POINT_NAME: str
    monitor: _PeriodicMonitorThread
    def __init__(self, monitor: _PeriodicMonitorThread) -> None: ...

class IEventLoopBlocked(Interface):
    greenlet: greenlet_t
    blocking_time: float
    info: Sequence[str]
    hub: Hub | None

@implementer(IEventLoopBlocked)
class EventLoopBlocked:
    greenlet: greenlet_t
    blocking_time: float
    info: Sequence[str]
    hub: Hub | None
    def __init__(self, greenlet: greenlet_t, blocking_time: float, info: Sequence[str], *, hub: Hub | None = None) -> None: ...

class IMemoryUsageThresholdExceeded(Interface):
    mem_usage: int
    max_allowed: int
    memory_info: pmem

class _AbstractMemoryEvent:
    mem_usage: int
    max_allowed: int
    memory_info: pmem
    def __init__(self, mem_usage: int, max_allowed: int, memory_info: pmem) -> None: ...

@implementer(IMemoryUsageThresholdExceeded)
class MemoryUsageThresholdExceeded(_AbstractMemoryEvent): ...

class IMemoryUsageUnderThreshold(Interface):
    mem_usage: int
    max_allowed: int
    max_memory_usage: int
    memory_info: pmem

@implementer(IMemoryUsageUnderThreshold)
class MemoryUsageUnderThreshold(_AbstractMemoryEvent):
    max_memory_usage: int
    def __init__(self, mem_usage: int, max_allowed: int, memory_info: pmem, max_usage: int) -> None: ...

class IGeventPatchEvent(Interface):
    source: object
    target: object

@implementer(IGeventPatchEvent)
class GeventPatchEvent:
    source: object
    target: object
    def __init__(self, source: object, target: object) -> None: ...

class IGeventWillPatchEvent(IGeventPatchEvent): ...
class DoNotPatch(BaseException): ...

@implementer(IGeventWillPatchEvent)
class GeventWillPatchEvent(GeventPatchEvent): ...

class IGeventDidPatchEvent(IGeventPatchEvent): ...

@implementer(IGeventWillPatchEvent)
class GeventDidPatchEvent(GeventPatchEvent): ...

class IGeventWillPatchModuleEvent(IGeventWillPatchEvent):
    source: ModuleType
    target: ModuleType
    module_name: str
    target_item_names: list[str]

@implementer(IGeventWillPatchModuleEvent)
class GeventWillPatchModuleEvent(GeventWillPatchEvent):
    ENTRY_POINT_NAME: str
    source: ModuleType
    target: ModuleType
    module_name: str
    target_item_names: list[str]
    def __init__(self, module_name: str, source: ModuleType, target: ModuleType, items: list[str]) -> None: ...

class IGeventDidPatchModuleEvent(IGeventDidPatchEvent):
    source: ModuleType
    target: ModuleType
    module_name: str

@implementer(IGeventDidPatchModuleEvent)
class GeventDidPatchModuleEvent(GeventDidPatchEvent):
    ENTRY_POINT_NAME: str
    source: ModuleType
    target: ModuleType
    module_name: str
    def __init__(self, module_name: str, source: ModuleType, target: ModuleType) -> None: ...

class IGeventWillPatchAllEvent(IGeventWillPatchEvent):
    patch_all_arguments: Mapping[str, Any]
    patch_all_kwargs: Mapping[str, Any]
    def will_patch_module(module_name: str) -> bool: ...

class _PatchAllMixin:
    def __init__(self, patch_all_arguments: Mapping[str, Any], patch_all_kwargs: Mapping[str, Any]) -> None: ...
    @property
    def patch_all_arguments(self) -> dict[str, Any]: ...  # safe to mutate, it's a copy
    @property
    def patch_all_kwargs(self) -> dict[str, Any]: ...  # safe to mutate, it's a copy

@implementer(IGeventWillPatchAllEvent)
class GeventWillPatchAllEvent(_PatchAllMixin, GeventWillPatchEvent):
    ENTRY_POINT_NAME: str
    def will_patch_module(self, module_name: str) -> bool: ...

class IGeventDidPatchBuiltinModulesEvent(IGeventDidPatchEvent):
    patch_all_arguments: Mapping[str, Any]
    patch_all_kwargs: Mapping[str, Any]

@implementer(IGeventDidPatchBuiltinModulesEvent)
class GeventDidPatchBuiltinModulesEvent(_PatchAllMixin, GeventDidPatchEvent):
    ENTRY_POINT_NAME: str

class IGeventDidPatchAllEvent(IGeventDidPatchEvent): ...

@implementer(IGeventDidPatchAllEvent)
class GeventDidPatchAllEvent(_PatchAllMixin, GeventDidPatchEvent):
    ENTRY_POINT_NAME: str

__all__ = [
    "subscribers",
    # monitor thread
    "IEventLoopBlocked",
    "EventLoopBlocked",
    "IMemoryUsageThresholdExceeded",
    "MemoryUsageThresholdExceeded",
    "IMemoryUsageUnderThreshold",
    "MemoryUsageUnderThreshold",
    # Hub
    "IPeriodicMonitorThread",
    "IPeriodicMonitorThreadStartedEvent",
    "PeriodicMonitorThreadStartedEvent",
    # monkey
    "IGeventPatchEvent",
    "GeventPatchEvent",
    "IGeventWillPatchEvent",
    "DoNotPatch",
    "GeventWillPatchEvent",
    "IGeventDidPatchEvent",
    "IGeventWillPatchModuleEvent",
    "GeventWillPatchModuleEvent",
    "IGeventDidPatchModuleEvent",
    "GeventDidPatchModuleEvent",
    "IGeventWillPatchAllEvent",
    "GeventWillPatchAllEvent",
    "IGeventDidPatchBuiltinModulesEvent",
    "GeventDidPatchBuiltinModulesEvent",
    "IGeventDidPatchAllEvent",
    "GeventDidPatchAllEvent",
]
