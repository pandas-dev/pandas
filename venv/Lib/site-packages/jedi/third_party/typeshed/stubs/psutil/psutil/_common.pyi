import enum
import io
import sys
import threading
from _typeshed import ConvertibleToFloat, FileDescriptorOrPath, Incomplete, StrOrBytesPath, SupportsWrite
from collections import defaultdict
from collections.abc import Callable
from socket import AF_INET6 as AF_INET6, AddressFamily, SocketKind
from typing import BinaryIO, Final, SupportsIndex, TypeVar, overload
from typing_extensions import ParamSpec

from . import _ntuples as ntp

POSIX: Final[bool]
WINDOWS: Final[bool]
LINUX: Final[bool]
MACOS: Final[bool]
OSX: Final[bool]
FREEBSD: Final[bool]
OPENBSD: Final[bool]
NETBSD: Final[bool]
BSD: Final[bool]
SUNOS: Final[bool]
AIX: Final[bool]

STATUS_RUNNING: Final = "running"
STATUS_SLEEPING: Final = "sleeping"
STATUS_DISK_SLEEP: Final = "disk-sleep"
STATUS_STOPPED: Final = "stopped"
STATUS_TRACING_STOP: Final = "tracing-stop"
STATUS_ZOMBIE: Final = "zombie"
STATUS_DEAD: Final = "dead"
STATUS_WAKE_KILL: Final = "wake-kill"
STATUS_WAKING: Final = "waking"
STATUS_IDLE: Final = "idle"
STATUS_LOCKED: Final = "locked"
STATUS_WAITING: Final = "waiting"
STATUS_SUSPENDED: Final = "suspended"
STATUS_PARKED: Final = "parked"

CONN_ESTABLISHED: Final = "ESTABLISHED"
CONN_SYN_SENT: Final = "SYN_SENT"
CONN_SYN_RECV: Final = "SYN_RECV"
CONN_FIN_WAIT1: Final = "FIN_WAIT1"
CONN_FIN_WAIT2: Final = "FIN_WAIT2"
CONN_TIME_WAIT: Final = "TIME_WAIT"
CONN_CLOSE: Final = "CLOSE"
CONN_CLOSE_WAIT: Final = "CLOSE_WAIT"
CONN_LAST_ACK: Final = "LAST_ACK"
CONN_LISTEN: Final = "LISTEN"
CONN_CLOSING: Final = "CLOSING"
CONN_NONE: Final = "NONE"

class NicDuplex(enum.IntEnum):
    NIC_DUPLEX_FULL = 2
    NIC_DUPLEX_HALF = 1
    NIC_DUPLEX_UNKNOWN = 0

NIC_DUPLEX_FULL: Final = NicDuplex.NIC_DUPLEX_FULL
NIC_DUPLEX_HALF: Final = NicDuplex.NIC_DUPLEX_HALF
NIC_DUPLEX_UNKNOWN: Final = NicDuplex.NIC_DUPLEX_UNKNOWN

class BatteryTime(enum.IntEnum):
    POWER_TIME_UNKNOWN = -1
    POWER_TIME_UNLIMITED = -2

POWER_TIME_UNKNOWN: Final = BatteryTime.POWER_TIME_UNKNOWN
POWER_TIME_UNLIMITED: Final = BatteryTime.POWER_TIME_UNLIMITED

ENCODING: Final[str]
ENCODING_ERRS: Final[str]

conn_tmap: dict[str, tuple[list[AddressFamily], list[SocketKind]]]

class Error(Exception): ...

class NoSuchProcess(Error):
    pid: int
    name: str | None
    msg: str
    def __init__(self, pid: int, name: str | None = None, msg: str | None = None) -> None: ...

class ZombieProcess(NoSuchProcess):
    ppid: int | None
    def __init__(self, pid: int, name: str | None = None, ppid: int | None = None, msg: str | None = None) -> None: ...

class AccessDenied(Error):
    pid: int | None
    name: str | None
    msg: str
    def __init__(self, pid: int | None = None, name: str | None = None, msg: str | None = None) -> None: ...

class TimeoutExpired(Error):
    seconds: float
    pid: int | None
    name: str | None
    msg: str
    def __init__(self, seconds: float, pid: int | None = None, name: str | None = None) -> None: ...

_P = ParamSpec("_P")
_R = TypeVar("_R")
_T = TypeVar("_T")

def usage_percent(used: ConvertibleToFloat, total: float, round_: SupportsIndex | None = None) -> float: ...

# returned function has `cache_clear()` attribute:
def memoize(fun: Callable[_P, _R]) -> Callable[_P, _R]: ...

# returned function has `cache_activate(proc)` and `cache_deactivate(proc)` attributes:
def memoize_when_activated(fun: Callable[_P, _R]) -> Callable[_P, _R]: ...
def isfile_strict(path: StrOrBytesPath) -> bool: ...
def path_exists_strict(path: StrOrBytesPath) -> bool: ...
def supports_ipv6() -> bool: ...
def parse_environ_block(data: str) -> dict[str, str]: ...
def sockfam_to_enum(num: int) -> AddressFamily: ...
def socktype_to_enum(num: int) -> SocketKind: ...
@overload
def conn_to_ntuple(
    fd: int,
    fam: int,
    type_: int,
    laddr: ntp.addr | tuple[str, int] | tuple[()],
    raddr: ntp.addr | tuple[str, int] | tuple[()],
    status: int | str,
    status_map: dict[int, str] | dict[str, str],
    pid: int,
) -> ntp.sconn: ...
@overload
def conn_to_ntuple(
    fd: int,
    fam: int,
    type_: int,
    laddr: ntp.addr | tuple[str, int] | tuple[()],
    raddr: ntp.addr | tuple[str, int] | tuple[()],
    status: int | str,
    status_map: dict[int, str] | dict[str, str],
    pid: None = None,
) -> ntp.pconn: ...
def deprecated_method(replacement: str) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]: ...

class _WrapNumbers:
    lock: threading.Lock
    cache: dict[str, dict[str, tuple[int, ...]]]
    reminders: dict[str, defaultdict[Incomplete, int]]
    reminder_keys: dict[str, defaultdict[Incomplete, set[Incomplete]]]
    def __init__(self) -> None: ...
    def run(self, input_dict: dict[str, tuple[int, ...]], name: str) -> dict[str, tuple[int, ...]]: ...
    def cache_clear(self, name: str | None = None) -> None: ...
    def cache_info(
        self,
    ) -> tuple[
        dict[str, dict[str, tuple[int, ...]]],
        dict[str, defaultdict[Incomplete, int]],
        dict[str, defaultdict[Incomplete, set[Incomplete]]],
    ]: ...

def wrap_numbers(input_dict: dict[str, tuple[int, ...]], name: str) -> dict[str, tuple[int, ...]]: ...
def open_binary(fname: FileDescriptorOrPath) -> BinaryIO: ...
def open_text(fname: FileDescriptorOrPath) -> io.TextIOWrapper: ...
@overload
def cat(fname: FileDescriptorOrPath, _open: Callable[[FileDescriptorOrPath], io.TextIOWrapper] = ...) -> str: ...
@overload
def cat(
    fname: FileDescriptorOrPath, fallback: _T = ..., _open: Callable[[FileDescriptorOrPath], io.TextIOWrapper] = ...
) -> str | _T: ...
@overload
def bcat(fname: FileDescriptorOrPath) -> str: ...
@overload
def bcat(fname: FileDescriptorOrPath, fallback: _T = ...) -> str | _T: ...
def bytes2human(n: int, format: str = "%(value).1f%(symbol)s") -> str: ...
def get_procfs_path() -> str: ...
def decode(s: bytes) -> str: ...
def term_supports_colors(file: SupportsWrite[str] = sys.stdout) -> bool: ...
def hilite(s: str, color: str | None = None, bold: bool = False) -> str: ...
def print_color(s: str, color: str | None = None, bold: bool = False, file: SupportsWrite[str] = sys.stdout) -> None: ...
def debug(msg: str | Exception) -> None: ...

__all__ = [
    # OS constants
    "FREEBSD",
    "BSD",
    "LINUX",
    "NETBSD",
    "OPENBSD",
    "MACOS",
    "OSX",
    "POSIX",
    "SUNOS",
    "WINDOWS",
    # connection constants
    "CONN_CLOSE",
    "CONN_CLOSE_WAIT",
    "CONN_CLOSING",
    "CONN_ESTABLISHED",
    "CONN_FIN_WAIT1",
    "CONN_FIN_WAIT2",
    "CONN_LAST_ACK",
    "CONN_LISTEN",
    "CONN_NONE",
    "CONN_SYN_RECV",
    "CONN_SYN_SENT",
    "CONN_TIME_WAIT",
    # net constants
    "NIC_DUPLEX_FULL",
    "NIC_DUPLEX_HALF",
    "NIC_DUPLEX_UNKNOWN",
    # process status constants
    "STATUS_DEAD",
    "STATUS_DISK_SLEEP",
    "STATUS_IDLE",
    "STATUS_LOCKED",
    "STATUS_RUNNING",
    "STATUS_SLEEPING",
    "STATUS_STOPPED",
    "STATUS_SUSPENDED",
    "STATUS_TRACING_STOP",
    "STATUS_WAITING",
    "STATUS_WAKE_KILL",
    "STATUS_WAKING",
    "STATUS_ZOMBIE",
    "STATUS_PARKED",
    # other constants
    "ENCODING",
    "ENCODING_ERRS",
    "AF_INET6",
    # utility functions
    "conn_tmap",
    "deprecated_method",
    "isfile_strict",
    "memoize",
    "parse_environ_block",
    "path_exists_strict",
    "usage_percent",
    "supports_ipv6",
    "sockfam_to_enum",
    "socktype_to_enum",
    "wrap_numbers",
    "open_text",
    "open_binary",
    "cat",
    "bcat",
    "bytes2human",
    "conn_to_ntuple",
    "debug",
    # shell utils
    "hilite",
    "term_supports_colors",
    "print_color",
]
