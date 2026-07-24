from collections.abc import Callable, Sequence
from typing import Any, Generic, NoReturn, Protocol, TypeVar, overload, type_check_only

from gevent._types import _Loop, _Resolver
from gevent.fileobject import _FileObjectType
from gevent.threadpool import ThreadPool

__all__ = ["config"]

_T = TypeVar("_T")

@type_check_only
class _SettingDescriptor(Protocol[_T]):
    @overload
    def __get__(self, obj: None, owner: type[Config]) -> property: ...
    @overload
    def __get__(self, obj: Config, owner: type[Config]) -> _T: ...
    def __set__(self, obj: Config, value: str | _T) -> None: ...

class SettingType(type):
    def fmt_desc(cls, desc: str) -> str: ...

def validate_invalid(value: object) -> NoReturn: ...
def validate_bool(value: str | bool) -> bool: ...
def validate_anything(value: _T) -> _T: ...

convert_str_value_as_is = validate_anything

class Setting(Generic[_T], metaclass=SettingType):
    order: int  # all subclasses have this
    name: str
    environment_key: str
    value: _T
    default: _T
    document: bool
    desc: str
    validate: Callable[[Any], _T]
    def get(self) -> _T: ...
    def set(self, val: str | _T) -> None: ...

class Config:
    settings: dict[str, Setting[Any]]
    def __init__(self) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __setattr__(self, name: str, value: object) -> None: ...
    def set(self, name: str, value: object) -> None: ...
    def __dir__(self) -> list[str]: ...
    def print_help(self) -> None: ...

    # we manually add properties for all the settings in this module
    # SettingType inserts a property into Config for every subclass of Setting
    resolver: _SettingDescriptor[type[_Resolver]]
    threadpool: _SettingDescriptor[type[Threadpool]]
    threadpool_idle_task_timeout: _SettingDescriptor[float]
    loop: _SettingDescriptor[type[_Loop]]
    format_context: _SettingDescriptor[Callable[[Any], str]]
    libev_backend: _SettingDescriptor[str | None]
    fileobject: _SettingDescriptor[_FileObjectType]
    disable_watch_children: _SettingDescriptor[bool]
    track_greenlet_tree: _SettingDescriptor[bool]
    monitor_thread: _SettingDescriptor[bool]
    max_blocking_time: _SettingDescriptor[float]
    memory_monitor_period: _SettingDescriptor[float]
    max_memory_usage: _SettingDescriptor[int | None]
    resolver_nameservers: _SettingDescriptor[Sequence[str] | str | None]
    resolver_timeout: _SettingDescriptor[float | None]
    # these get parsed by gevent.resolver.cares.channel so the Setting does not
    # perform any conversion, but we know at least what types can be valid
    ares_flags: _SettingDescriptor[str | int | None]
    ares_timeout: _SettingDescriptor[str | float | None]
    ares_tries: _SettingDescriptor[str | int | None]
    ares_ndots: _SettingDescriptor[str | int | None]
    ares_udp_port: _SettingDescriptor[str | int | None]
    ares_tcp_port: _SettingDescriptor[str | int | None]
    ares_servers: _SettingDescriptor[Sequence[str] | str | None]
    print_blocking_reports: _SettingDescriptor[bool]

class ImportableSetting(Generic[_T]):
    default: str | Sequence[str]
    shortname_map: dict[str, str]
    def validate(self, value: str | _T) -> _T: ...
    def get_options(self) -> dict[str, _T]: ...

class BoolSettingMixin:
    @staticmethod
    def validate(value: str | bool) -> bool: ...

class IntSettingMixin:
    @staticmethod
    def validate(value: int) -> int: ...

class _PositiveValueMixin(Generic[_T]):
    @staticmethod
    def validate(value: _T) -> _T: ...

class FloatSettingMixin(_PositiveValueMixin[float]): ...
class ByteCountSettingMixin(_PositiveValueMixin[int]): ...

class Resolver(ImportableSetting[type[_Resolver]], Setting[type[_Resolver]]):
    desc: str
    default: list[str]  # type: ignore[assignment]
    shortname_map: dict[str, str]

class Threadpool(ImportableSetting[type[ThreadPool]], Setting[type[ThreadPool]]):
    desc: str
    default: str  # type: ignore[assignment]

class ThreadpoolIdleTaskTimeout(FloatSettingMixin, Setting[float]):
    document: bool
    desc: str
    default: float

class Loop(ImportableSetting[type[_Loop]], Setting[type[_Loop]]):
    desc: str
    default: list[str]  # type: ignore[assignment]
    shortname_map: dict[str, str]

class FormatContext(ImportableSetting[Callable[[Any], str]], Setting[Callable[[Any], str]]):
    default: str  # type: ignore[assignment]

class LibevBackend(Setting[str | None]):
    desc: str
    default: None

class FileObject(ImportableSetting[_FileObjectType], Setting[_FileObjectType]):
    desc: str
    default: list[str]  # type: ignore[assignment]
    shortname_map: dict[str, str]

class WatchChildren(BoolSettingMixin, Setting[bool]):
    desc: str
    default: bool

class TrackGreenletTree(BoolSettingMixin, Setting[bool]):
    default: bool
    desc: str

class MonitorThread(BoolSettingMixin, Setting[bool]):
    default: bool
    desc: str

class MaxBlockingTime(FloatSettingMixin, Setting[float]):
    default: float
    desc: str

class PrintBlockingReports(BoolSettingMixin, Setting[bool]):
    default: bool
    desc: str

class MonitorMemoryPeriod(FloatSettingMixin, Setting[float]):
    default: int
    desc: str

class MonitorMemoryMaxUsage(ByteCountSettingMixin, Setting[int | None]):
    default: None
    desc: str

class AresSettingMixin:
    document: bool
    @property
    def kwarg_name(self) -> str: ...
    validate: Any  # we just want this to mixin without errors

class AresFlags(AresSettingMixin, Setting[str | int | None]):
    default: None

class AresTimeout(AresSettingMixin, Setting[str | float | None]):
    document: bool
    default: None
    desc: str

class AresTries(AresSettingMixin, Setting[str | int | None]):
    default: None

class AresNdots(AresSettingMixin, Setting[str | int | None]):
    default: None

class AresUDPPort(AresSettingMixin, Setting[str | int | None]):
    default: None

class AresTCPPort(AresSettingMixin, Setting[str | int | None]):
    default: None

class AresServers(AresSettingMixin, Setting[Sequence[str] | str | None]):
    document: bool
    default: None
    desc: str

class ResolverNameservers(AresSettingMixin, Setting[Sequence[str] | str | None]):
    document: bool
    default: None
    desc: str
    @property
    def kwarg_name(self) -> str: ...

class ResolverTimeout(FloatSettingMixin, AresSettingMixin, Setting[float | None]):
    document: bool
    desc: str
    @property
    def kwarg_name(self) -> str: ...

config: Config = ...
