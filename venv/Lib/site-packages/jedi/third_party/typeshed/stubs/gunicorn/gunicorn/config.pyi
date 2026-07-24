import argparse
from _typeshed import ConvertibleToInt
from collections.abc import Callable, Container
from ssl import SSLContext, _SSLMethod
from typing import Annotated, Any, ClassVar, Final, overload
from typing_extensions import TypeAlias

from gunicorn.arbiter import Arbiter
from gunicorn.glogging import Logger as GLogger
from gunicorn.http import Request
from gunicorn.http.wsgi import Response
from gunicorn.workers.base import Worker

from ._types import _AddressType, _EnvironType

_ConfigValueType: TypeAlias = Any
# Any: Value type depends on the setting's validator (e.g., validate_bool, validate_pos_int).
# Maybe types: bool, int, str, list[str], dict[str, Any], Callable[..., object], type[object], ssl.PROTOCOL_*.
# Validator ensures type correctness at runtime.

_OnStartingHookType: TypeAlias = Callable[[Arbiter], object]
_OnReloadHookType: TypeAlias = Callable[[Arbiter], object]
_WhenReadyHookType: TypeAlias = Callable[[Arbiter], object]
_PreForkHookType: TypeAlias = Callable[[Arbiter, Worker], object]
_PostForkHookType: TypeAlias = Callable[[Arbiter, Worker], object]
_PostWorkerInitHookType: TypeAlias = Callable[[Worker], object]
_WorkerIntHookType: TypeAlias = Callable[[Worker], object]
_WorkerAbortHookType: TypeAlias = Callable[[Worker], object]
_PreExecHookType: TypeAlias = Callable[[Arbiter], object]
_PreRequestHookType: TypeAlias = Callable[[Worker, Request], object]
_PostRequestHookType: TypeAlias = Callable[[Worker, Request, _EnvironType, Response], object]
_ChildExitHookType: TypeAlias = Callable[[Arbiter, Worker], object]
_WorkerExitHookType: TypeAlias = Callable[[Arbiter, Worker], object]
_NumWorkersChangedHookType: TypeAlias = Callable[[Arbiter, int, int | None], object]
_OnExitHookType: TypeAlias = Callable[[Arbiter], object]
_SSLContextHookType: TypeAlias = Callable[[Config, Callable[[], SSLContext]], SSLContext]
_OnDirtyStartingHookType: TypeAlias = Callable[[Arbiter], object]
_DirtyPostForkHookType: TypeAlias = Callable[[Arbiter, Worker], object]
_DirtyWorkerInitHookType: TypeAlias = Callable[[Worker], object]
_DirtyWorkerExitHookType: TypeAlias = Callable[[Arbiter, Worker], object]

_HookType: TypeAlias = (
    _OnStartingHookType
    | _OnReloadHookType
    | _WhenReadyHookType
    | _PreForkHookType
    | _PostForkHookType
    | _PostWorkerInitHookType
    | _WorkerIntHookType
    | _WorkerAbortHookType
    | _PreExecHookType
    | _PreRequestHookType
    | _PostRequestHookType
    | _ChildExitHookType
    | _WorkerExitHookType
    | _NumWorkersChangedHookType
    | _OnExitHookType
    | _SSLContextHookType
    | _OnDirtyStartingHookType
    | _DirtyPostForkHookType
    | _DirtyWorkerInitHookType
    | _DirtyWorkerExitHookType
)
# Validators
_BoolValidatorType: TypeAlias = Callable[[bool | str | None], bool | None]
_StringValidatorType: TypeAlias = Callable[[str | None], str | None]
_ListStringValidatorType: TypeAlias = Callable[[str | list[str] | None], list[str]]
_IntValidatorType: TypeAlias = Callable[[ConvertibleToInt], int]
_DictValidatorType: TypeAlias = Callable[[dict[str, Any]], dict[str, Any]]
_ClassValidatorType: TypeAlias = Callable[[object | str | None], type[Any] | None]
_UserGroupValidatorType: TypeAlias = Callable[[str | int | None], int]
_AddressValidatorType: TypeAlias = Callable[[str | None], _AddressType | None]
_CallableValidatorType: TypeAlias = Callable[[str | _HookType], _HookType]
_ProxyProtocolValidatorType: TypeAlias = Callable[[str | bool | None], str]
_ASGILoopValidatorType: TypeAlias = Callable[[str | None], str]
_ASGILifespanValidatorType: TypeAlias = Callable[[str | None], str]
_HTTP2FrameSizeValidatorType: TypeAlias = Callable[[ConvertibleToInt], int]
_HTTPProtocolsValidatorType: TypeAlias = Callable[[str | None], list[str]]
_HttpParserValidatorType: TypeAlias = Callable[[str | None], str]

_ValidatorType: TypeAlias = (  # noqa: Y047
    _BoolValidatorType
    | _StringValidatorType
    | _ListStringValidatorType
    | _IntValidatorType
    | _DictValidatorType
    | _ClassValidatorType
    | _UserGroupValidatorType
    | _AddressValidatorType
    | _CallableValidatorType
    | _ProxyProtocolValidatorType
    | _ASGILoopValidatorType
    | _ASGILifespanValidatorType
    | _HTTP2FrameSizeValidatorType
    | _HTTPProtocolsValidatorType
    | _HttpParserValidatorType
)

KNOWN_SETTINGS: list[Setting]
PLATFORM: str

def make_settings(ignore: Container[Setting] | None = None) -> dict[str, Setting]: ...
def auto_int(_: Any, x: str) -> int: ...

class Config:
    settings: dict[str, Setting]
    usage: str | None
    prog: str | None
    env_orig: dict[str, str]

    def __init__(self, usage: str | None = None, prog: str | None = None) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __setattr__(self, name: str, value: Any) -> None: ...
    def set(self, name: str, value: _ConfigValueType) -> None: ...
    def get_cmd_args_from_env(self) -> list[str]: ...
    def parser(self) -> argparse.ArgumentParser: ...
    @property
    def worker_class_str(self) -> str: ...
    @property
    def worker_class(self) -> type[Worker]: ...
    @property
    def address(self) -> list[_AddressType]: ...
    @property
    def uid(self) -> int: ...
    @property
    def gid(self) -> int: ...
    @property
    def proc_name(self) -> str | None: ...
    @property
    def logger_class(self) -> type[GLogger]: ...
    @property
    def is_ssl(self) -> bool: ...
    @property
    def ssl_options(self) -> dict[str, Any]: ...
    @property
    def env(self) -> dict[str, str]: ...
    @property
    def sendfile(self) -> bool: ...
    @property
    def reuse_port(self) -> bool: ...
    @property
    def paste_global_conf(self) -> dict[str, str] | None: ...

class SettingMeta(type):
    def __new__(cls, name: str, bases: tuple[type, ...], attrs: dict[str, Any]) -> SettingMeta: ...
    def fmt_desc(cls, desc: str) -> None: ...

class Setting(metaclass=SettingMeta):
    name: ClassVar[str | None]
    value: _ConfigValueType
    section: ClassVar[str | None]
    cli: ClassVar[list[str] | None]
    validator: ClassVar[Callable[..., Any] | None]  # See `_ValidatorType`
    type: ClassVar[argparse._ActionType | None]
    meta: ClassVar[str | None]
    action: ClassVar[str | None]
    default: ClassVar[Any]
    short: ClassVar[str | None]
    desc: ClassVar[str | None]
    nargs: ClassVar[int | str | None]
    const: ClassVar[bool | str | None]
    order: ClassVar[int]

    def __init__(self) -> None: ...
    def add_option(self, parser: argparse.ArgumentParser) -> None: ...
    def copy(self) -> Setting: ...
    def get(self) -> _ConfigValueType: ...
    def set(self, val: _ConfigValueType) -> None: ...
    def __lt__(self, other: Setting) -> bool: ...

    __cmp__ = __lt__

@overload
def validate_bool(val: bool) -> bool: ...
@overload
def validate_bool(val: None) -> None: ...
@overload
def validate_bool(val: Annotated[str, "Case-insensitive boolean string ('true'/'false' in any case)"]) -> bool: ...
def validate_dict(val: dict[str, Any]) -> dict[str, Any]: ...
def validate_pos_int(val: ConvertibleToInt) -> int: ...
def validate_http2_frame_size(val: ConvertibleToInt) -> int: ...
def validate_ssl_version(val: _SSLMethod) -> _SSLMethod: ...
@overload
def validate_string(val: str) -> str: ...
@overload
def validate_string(val: None) -> None: ...
@overload
def validate_file_exists(val: str) -> str: ...
@overload
def validate_file_exists(val: None) -> None: ...
def validate_list_string(val: str | list[str] | None) -> list[str]: ...
def validate_list_of_existing_files(val: str | list[str] | None) -> list[str]: ...
def validate_string_to_addr_list(val: str | None) -> list[str]: ...
def validate_string_to_list(val: str | None) -> list[str]: ...
@overload
def validate_class(val: str) -> str: ...
@overload
def validate_class(val: None) -> None: ...
@overload
def validate_class(val: object) -> object: ...
def validate_callable(arity: int) -> _CallableValidatorType: ...
def validate_user(val: int | str | None) -> int: ...
def validate_group(val: int | str | None) -> int: ...
def validate_post_request(val: str | _HookType) -> _PostRequestHookType: ...
def validate_chdir(val: str) -> str: ...
@overload
def validate_statsd_address(val: str) -> _AddressType: ...
@overload
def validate_statsd_address(val: None) -> None: ...
def validate_reload_engine(val: str) -> str: ...
@overload
def validate_header_map_behaviour(val: str) -> str: ...
@overload
def validate_header_map_behaviour(val: None) -> None: ...
def get_default_config_file() -> str | None: ...

class ConfigFile(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class WSGIApp(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class Bind(Setting):
    name: ClassVar[str]
    action: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[list[str]]
    desc: ClassVar[str]

class Backlog(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class Workers(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class WorkerClass(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ClassValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class WorkerThreads(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class WorkerConnections(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class MaxRequests(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class MaxRequestsJitter(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class Timeout(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class GracefulTimeout(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class Keepalive(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class LimitRequestLine(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class LimitRequestFields(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class LimitRequestFieldSize(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class Reload(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class ReloadEngine(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[Callable[[str], str]]
    default: ClassVar[str]
    desc: ClassVar[str]

class ReloadExtraFiles(Setting):
    name: ClassVar[str]
    action: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[list[str]]
    desc: ClassVar[str]

class Spew(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class ConfigCheck(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class PrintConfig(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class PreloadApp(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class Sendfile(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    const: ClassVar[bool]
    desc: ClassVar[str]

class ReusePort(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class Chdir(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[Callable[[str], str]]
    default: ClassVar[str]
    default_doc: ClassVar[str]
    desc: ClassVar[str]

class Daemon(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class Env(Setting):
    name: ClassVar[str]
    action: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[list[str]]
    desc: ClassVar[str]

class Pidfile(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class WorkerTmpDir(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class User(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_UserGroupValidatorType]
    default: ClassVar[int]
    default_doc: ClassVar[str]
    desc: ClassVar[str]

class Group(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_UserGroupValidatorType]
    default: ClassVar[int]
    default_doc: ClassVar[str]
    desc: ClassVar[str]

class Umask(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[Callable[[Any, str], int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class Initgroups(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class TmpUploadDir(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class SecureSchemeHeader(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_DictValidatorType]
    default: ClassVar[dict[str, str]]
    desc: ClassVar[str]

class ForwardedAllowIPS(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class AccessLog(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class DisableRedirectAccessToSyslog(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class AccessLogFormat(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class ErrorLog(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class Loglevel(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class CaptureOutput(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class LoggerClass(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ClassValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class LogConfig(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class LogConfigDict(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_DictValidatorType]
    default: ClassVar[dict[str, Any]]
    desc: ClassVar[str]

class LogConfigJson(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class SyslogTo(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]
    default_doc: ClassVar[str]

class Syslog(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class SyslogPrefix(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class SyslogFacility(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class EnableStdioInheritance(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    default: ClassVar[bool]
    action: ClassVar[str]
    desc: ClassVar[str]

class StatsdHost(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    default: ClassVar[None]
    validator: ClassVar[_AddressValidatorType]
    desc: ClassVar[str]

class DogstatsdTags(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    default: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    desc: ClassVar[str]

class StatsdPrefix(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    default: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    desc: ClassVar[str]

class BacklogMetric(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    default: ClassVar[bool]
    action: ClassVar[str]
    desc: ClassVar[str]

class Procname(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class DefaultProcName(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class PythonPath(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class Paste(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class OnStarting(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_OnStartingHookType]
    desc: ClassVar[str]

    def on_starting(server: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class OnReload(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_OnReloadHookType]
    desc: ClassVar[str]

    def on_reload(server: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class WhenReady(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_WhenReadyHookType]
    desc: ClassVar[str]

    def when_ready(server: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class Prefork(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PreForkHookType]
    desc: ClassVar[str]

    def pre_fork(server: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class Postfork(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PostForkHookType]
    desc: ClassVar[str]

    def post_fork(server: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class PostWorkerInit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PostWorkerInitHookType]
    desc: ClassVar[str]

    def post_worker_init(worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class WorkerInt(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_WorkerIntHookType]
    desc: ClassVar[str]

    def worker_int(worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class WorkerAbort(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_WorkerAbortHookType]
    desc: ClassVar[str]

    def worker_abort(worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class PreExec(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PreExecHookType]
    desc: ClassVar[str]

    def pre_exec(server: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class PreRequest(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PreRequestHookType]
    desc: ClassVar[str]

    def pre_request(worker: Worker, req: Request) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class PostRequest(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_PostRequestHookType]
    desc: ClassVar[str]

    def post_request(worker: Worker, req: Request, environ: _EnvironType, resp: Response) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class ChildExit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_ChildExitHookType]
    desc: ClassVar[str]

    def child_exit(server: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class WorkerExit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_WorkerExitHookType]
    desc: ClassVar[str]

    def worker_exit(server: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class NumWorkersChanged(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_NumWorkersChangedHookType]
    desc: ClassVar[str]

    def nworkers_changed(server: Arbiter, new_value: int, old_value: int | None) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class OnExit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    default: ClassVar[_OnExitHookType]
    desc: ClassVar[str]

    def on_exit(server: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class NewSSLContext(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_SSLContextHookType]
    desc: ClassVar[str]

    def ssl_context(config: Config, default_ssl_context_factory: Callable[[], SSLContext]) -> SSLContext: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

def validate_proxy_protocol(val: str | bool | None) -> str: ...

class ProxyProtocol(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ProxyProtocolValidatorType]
    default: ClassVar[str]
    nargs: ClassVar[str]
    const: ClassVar[str]
    desc: ClassVar[str]

class ProxyAllowFrom(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class Protocol(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class UWSGIAllowFrom(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class KeyFile(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class CertFile(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class SSLVersion(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[Callable[[_SSLMethod], _SSLMethod]]
    default: ClassVar[_SSLMethod]
    desc: ClassVar[str]

class CertReqs(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_IntValidatorType]
    default: ClassVar[int]
    desc: ClassVar[str]

class CACerts(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

class SuppressRaggedEOFs(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    action: ClassVar[str]
    default: ClassVar[bool]
    validator: ClassVar[_BoolValidatorType]
    desc: ClassVar[str]

class DoHandshakeOnConnect(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class Ciphers(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[None]
    desc: ClassVar[str]

VALID_HTTP_PROTOCOLS: Final[frozenset[str]]
ALPN_PROTOCOL_MAP: Final[dict[str, str]]

def validate_http_protocols(val: str | None) -> list[str]: ...

class HTTPProtocols(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_HTTPProtocolsValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class HTTP2MaxConcurrentStreams(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class HTTP2InitialWindowSize(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class HTTP2MaxFrameSize(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_HTTP2FrameSizeValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class HTTP2MaxHeaderListSize(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class PasteGlobalConf(Setting):
    name: ClassVar[str]
    action: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[list[str]]
    desc: ClassVar[str]

class PermitObsoleteFolding(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class StripHeaderSpaces(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class PermitUnconventionalHTTPMethod(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class PermitUnconventionalHTTPVersion(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class CasefoldHTTPMethod(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]

class ForwarderHeaders(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class HeaderMap(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

def validate_asgi_loop(val: str | None) -> str: ...
def validate_asgi_lifespan(val: str | None) -> str: ...
def validate_http_parser(val: str | None) -> str: ...

class ASGILoop(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ASGILoopValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class ASGILifespan(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_ASGILifespanValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class ASGIDisconnectGracePeriod(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class HttpParser(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_HttpParserValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class RootPath(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    desc: ClassVar[str]

class DirtyApps(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    action: ClassVar[str]
    meta: ClassVar[str]
    validator: ClassVar[_ListStringValidatorType]
    default: ClassVar[list[str]]
    desc: ClassVar[str]

class DirtyWorkers(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class DirtyTimeout(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class DirtyThreads(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class DirtyGracefulTimeout(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[type[int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class OnDirtyStarting(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_OnDirtyStartingHookType]
    desc: ClassVar[str]

    def on_dirty_starting(arbiter: Arbiter) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class DirtyPostFork(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_DirtyPostForkHookType]
    desc: ClassVar[str]

    def dirty_post_fork(arbiter: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class DirtyWorkerInit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_DirtyWorkerInitHookType]
    desc: ClassVar[str]

    def dirty_worker_init(worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class DirtyWorkerExit(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    validator: ClassVar[_CallableValidatorType]
    type: ClassVar[Callable[..., Any]]
    default: ClassVar[_DirtyWorkerExitHookType]
    desc: ClassVar[str]

    def dirty_worker_exit(arbiter: Arbiter, worker: Worker) -> None: ...  # type: ignore[misc] # pyright: ignore[reportGeneralTypeIssues]

class ControlSocket(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_StringValidatorType]
    default: ClassVar[str]
    default_doc: ClassVar[str]
    desc: ClassVar[str]

class ControlSocketMode(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    meta: ClassVar[str]
    validator: ClassVar[_IntValidatorType]
    type: ClassVar[Callable[[Any, str], int]]
    default: ClassVar[int]
    desc: ClassVar[str]

class ControlSocketDisable(Setting):
    name: ClassVar[str]
    section: ClassVar[str]
    cli: ClassVar[list[str]]
    validator: ClassVar[_BoolValidatorType]
    action: ClassVar[str]
    default: ClassVar[bool]
    desc: ClassVar[str]
