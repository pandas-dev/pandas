import sys
from _typeshed import Incomplete, StrOrBytesPath
from collections.abc import Callable, Collection, Iterable, Iterator
from contextlib import AbstractContextManager
from subprocess import _CMD, _ENV, _FILE
from types import TracebackType
from typing import Any, Literal, Protocol, overload, type_check_only
from typing_extensions import Self, TypeAlias, deprecated

from psutil._common import (
    AIX as AIX,
    BSD as BSD,
    CONN_CLOSE as CONN_CLOSE,
    CONN_CLOSE_WAIT as CONN_CLOSE_WAIT,
    CONN_CLOSING as CONN_CLOSING,
    CONN_ESTABLISHED as CONN_ESTABLISHED,
    CONN_FIN_WAIT1 as CONN_FIN_WAIT1,
    CONN_FIN_WAIT2 as CONN_FIN_WAIT2,
    CONN_LAST_ACK as CONN_LAST_ACK,
    CONN_LISTEN as CONN_LISTEN,
    CONN_NONE as CONN_NONE,
    CONN_SYN_RECV as CONN_SYN_RECV,
    CONN_SYN_SENT as CONN_SYN_SENT,
    CONN_TIME_WAIT as CONN_TIME_WAIT,
    FREEBSD as FREEBSD,
    LINUX as LINUX,
    MACOS as MACOS,
    NETBSD as NETBSD,
    NIC_DUPLEX_FULL as NIC_DUPLEX_FULL,
    NIC_DUPLEX_HALF as NIC_DUPLEX_HALF,
    NIC_DUPLEX_UNKNOWN as NIC_DUPLEX_UNKNOWN,
    OPENBSD as OPENBSD,
    OSX as OSX,
    POSIX as POSIX,
    POWER_TIME_UNKNOWN as POWER_TIME_UNKNOWN,
    POWER_TIME_UNLIMITED as POWER_TIME_UNLIMITED,
    STATUS_DEAD as STATUS_DEAD,
    STATUS_DISK_SLEEP as STATUS_DISK_SLEEP,
    STATUS_IDLE as STATUS_IDLE,
    STATUS_LOCKED as STATUS_LOCKED,
    STATUS_PARKED as STATUS_PARKED,
    STATUS_RUNNING as STATUS_RUNNING,
    STATUS_SLEEPING as STATUS_SLEEPING,
    STATUS_STOPPED as STATUS_STOPPED,
    STATUS_TRACING_STOP as STATUS_TRACING_STOP,
    STATUS_WAITING as STATUS_WAITING,
    STATUS_WAKING as STATUS_WAKING,
    STATUS_ZOMBIE as STATUS_ZOMBIE,
    SUNOS as SUNOS,
    WINDOWS as WINDOWS,
    AccessDenied as AccessDenied,
    Error as Error,
    NoSuchProcess as NoSuchProcess,
    TimeoutExpired as TimeoutExpired,
    ZombieProcess as ZombieProcess,
)

from . import _ntuples as _ntp

if sys.platform == "linux":
    from ._pslinux import (
        IOPRIO_CLASS_BE as IOPRIO_CLASS_BE,
        IOPRIO_CLASS_IDLE as IOPRIO_CLASS_IDLE,
        IOPRIO_CLASS_NONE as IOPRIO_CLASS_NONE,
        IOPRIO_CLASS_RT as IOPRIO_CLASS_RT,
    )
    def sensors_temperatures(fahrenheit: bool = False) -> dict[str, list[_ntp.shwtemp]]: ...
    def sensors_fans() -> dict[str, list[_ntp.sfan]]: ...
    PROCFS_PATH: str
    RLIMIT_AS: int
    RLIMIT_CORE: int
    RLIMIT_CPU: int
    RLIMIT_DATA: int
    RLIMIT_FSIZE: int
    RLIMIT_LOCKS: int
    RLIMIT_MEMLOCK: int
    RLIMIT_MSGQUEUE: int
    RLIMIT_NICE: int
    RLIMIT_NOFILE: int
    RLIMIT_NPROC: int
    RLIMIT_RSS: int
    RLIMIT_RTPRIO: int
    RLIMIT_RTTIME: int
    RLIMIT_SIGPENDING: int
    RLIMIT_STACK: int
    RLIM_INFINITY: int
if sys.platform == "win32":
    from ._psutil_windows import (
        ABOVE_NORMAL_PRIORITY_CLASS as ABOVE_NORMAL_PRIORITY_CLASS,
        BELOW_NORMAL_PRIORITY_CLASS as BELOW_NORMAL_PRIORITY_CLASS,
        HIGH_PRIORITY_CLASS as HIGH_PRIORITY_CLASS,
        IDLE_PRIORITY_CLASS as IDLE_PRIORITY_CLASS,
        NORMAL_PRIORITY_CLASS as NORMAL_PRIORITY_CLASS,
        REALTIME_PRIORITY_CLASS as REALTIME_PRIORITY_CLASS,
    )
    from ._pswindows import (
        CONN_DELETE_TCB as CONN_DELETE_TCB,
        IOPRIO_HIGH as IOPRIO_HIGH,
        IOPRIO_LOW as IOPRIO_LOW,
        IOPRIO_NORMAL as IOPRIO_NORMAL,
        IOPRIO_VERYLOW as IOPRIO_VERYLOW,
        win_service_get as win_service_get,
        win_service_iter as win_service_iter,
    )

# Linux + glibc, Windows, macOS, FreeBSD, NetBSD:
def heap_info() -> _ntp.pheap: ...
def heap_trim() -> None: ...

if sys.platform == "linux":
    from ._pslinux import sensors_battery as sensors_battery
elif sys.platform == "darwin":
    from ._psosx import sensors_battery as sensors_battery
elif sys.platform == "win32":
    from ._pswindows import sensors_battery as sensors_battery
else:
    def sensors_battery(): ...

AF_LINK: int
version_info: tuple[int, int, int]
__version__: str
__author__: str

_Status: TypeAlias = Literal[
    "running",
    "sleeping",
    "disk-sleep",
    "stopped",
    "tracing-stop",
    "zombie",
    "dead",
    "wake-kill",
    "waking",
    "idle",
    "locked",
    "waiting",
    "suspended",
    "parked",
]

class Process:
    def __init__(self, pid: int | None = None) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    @property
    def pid(self) -> int: ...
    # Only present if attrs argument is passed to process_iter
    info: dict[str, Incomplete]
    def oneshot(self) -> AbstractContextManager[None]: ...
    def as_dict(
        self, attrs: list[str] | tuple[str, ...] | set[str] | frozenset[str] | None = None, ad_value=None
    ) -> dict[str, Incomplete]: ...
    def parent(self) -> Process | None: ...
    def parents(self) -> list[Process]: ...
    def is_running(self) -> bool: ...
    def ppid(self) -> int: ...
    def name(self) -> str: ...
    def exe(self) -> str: ...
    def cmdline(self) -> list[str]: ...
    def status(self) -> _Status: ...
    def username(self) -> str: ...
    def create_time(self) -> float: ...
    def cwd(self) -> str: ...
    def nice(self, value: int | None = None) -> int: ...
    if sys.platform != "win32":
        def uids(self) -> _ntp.puids: ...
        def gids(self) -> _ntp.pgids: ...
        def terminal(self) -> str: ...
        def num_fds(self) -> int: ...
    if sys.platform != "darwin":
        def io_counters(self) -> _ntp.pio: ...
        def ionice(self, ioclass: int | None = None, value: int | None = None) -> _ntp.pionice: ...
        @overload
        def cpu_affinity(self, cpus: None = None) -> list[int]: ...
        @overload
        def cpu_affinity(self, cpus: list[int]) -> None: ...
        def memory_maps(self, grouped: bool = True) -> list[Incomplete]: ...
    if sys.platform == "linux":
        def rlimit(self, resource: int, limits: tuple[int, int] | None = None) -> tuple[int, int]: ...
        def cpu_num(self) -> int: ...

    def environ(self) -> dict[str, str]: ...
    if sys.platform == "win32":
        def num_handles(self) -> int: ...

    def num_ctx_switches(self) -> _ntp.pctxsw: ...
    def num_threads(self) -> int: ...
    def threads(self) -> list[_ntp.pthread]: ...
    def children(self, recursive: bool = False) -> list[Process]: ...
    def cpu_percent(self, interval: float | None = None) -> float: ...
    def cpu_times(self) -> _ntp.pcputimes: ...
    def memory_info(self) -> _ntp.pmem: ...
    def memory_full_info(self) -> _ntp.pfullmem: ...
    def memory_percent(self, memtype: str = "rss") -> float: ...
    def open_files(self) -> list[_ntp.popenfile]: ...
    @deprecated('use "net_connections" method instead')
    def connections(self, kind: str = "inet") -> list[_ntp.pconn]: ...
    def send_signal(self, sig: int) -> None: ...
    def suspend(self) -> None: ...
    def resume(self) -> None: ...
    def terminate(self) -> None: ...
    def kill(self) -> None: ...
    def wait(self, timeout: float | None = None) -> int: ...
    def net_connections(self, kind: str = "inet") -> list[_ntp.pconn]: ...

class Popen(Process):
    # sync with subprocess.Popen.__init__:
    if sys.version_info >= (3, 11):
        def __init__(
            self,
            args: _CMD,
            bufsize: int = -1,
            executable: StrOrBytesPath | None = None,
            stdin: _FILE | None = None,
            stdout: _FILE | None = None,
            stderr: _FILE | None = None,
            preexec_fn: Callable[[], object] | None = None,
            close_fds: bool = True,
            shell: bool = False,
            cwd: StrOrBytesPath | None = None,
            env: _ENV | None = None,
            universal_newlines: bool | None = None,
            startupinfo: Any | None = None,
            creationflags: int = 0,
            restore_signals: bool = True,
            start_new_session: bool = False,
            pass_fds: Collection[int] = (),
            *,
            text: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            user: str | int | None = None,
            group: str | int | None = None,
            extra_groups: Iterable[str | int] | None = None,
            umask: int = -1,
            pipesize: int = -1,
            process_group: int | None = None,
        ) -> None: ...
    elif sys.version_info >= (3, 10):
        def __init__(
            self,
            args: _CMD,
            bufsize: int = -1,
            executable: StrOrBytesPath | None = None,
            stdin: _FILE | None = None,
            stdout: _FILE | None = None,
            stderr: _FILE | None = None,
            preexec_fn: Callable[[], object] | None = None,
            close_fds: bool = True,
            shell: bool = False,
            cwd: StrOrBytesPath | None = None,
            env: _ENV | None = None,
            universal_newlines: bool | None = None,
            startupinfo: Any | None = None,
            creationflags: int = 0,
            restore_signals: bool = True,
            start_new_session: bool = False,
            pass_fds: Collection[int] = (),
            *,
            text: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            user: str | int | None = None,
            group: str | int | None = None,
            extra_groups: Iterable[str | int] | None = None,
            umask: int = -1,
            pipesize: int = -1,
        ) -> None: ...
    else:
        def __init__(
            self,
            args: _CMD,
            bufsize: int = -1,
            executable: StrOrBytesPath | None = None,
            stdin: _FILE | None = None,
            stdout: _FILE | None = None,
            stderr: _FILE | None = None,
            preexec_fn: Callable[[], object] | None = None,
            close_fds: bool = True,
            shell: bool = False,
            cwd: StrOrBytesPath | None = None,
            env: _ENV | None = None,
            universal_newlines: bool | None = None,
            startupinfo: Any | None = None,
            creationflags: int = 0,
            restore_signals: bool = True,
            start_new_session: bool = False,
            pass_fds: Collection[int] = (),
            *,
            text: bool | None = None,
            encoding: str | None = None,
            errors: str | None = None,
            user: str | int | None = None,
            group: str | int | None = None,
            extra_groups: Iterable[str | int] | None = None,
            umask: int = -1,
        ) -> None: ...

    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
    def __getattribute__(self, name: str) -> Any: ...
    def __dir__(self) -> list[str]: ...

@type_check_only
class _ProcessIterCallable(Protocol):
    def __call__(
        self, attrs: list[str] | tuple[str, ...] | set[str] | frozenset[str] | None = None, ad_value=None
    ) -> Iterator[Process]: ...
    def cache_clear(self) -> None: ...

def pids() -> list[int]: ...
def pid_exists(pid: int) -> bool: ...

process_iter: _ProcessIterCallable

def wait_procs(
    procs: Iterable[Process], timeout: float | None = None, callback: Callable[[Process], object] | None = None
) -> tuple[list[Process], list[Process]]: ...
def cpu_count(logical: bool = True) -> int | None: ...
@overload
def cpu_freq(percpu: Literal[False] = False) -> _ntp.scpufreq: ...
@overload
def cpu_freq(percpu: Literal[True]) -> list[_ntp.scpufreq]: ...
@overload
def cpu_times(percpu: Literal[False] = False) -> _ntp.scputimes: ...
@overload
def cpu_times(percpu: Literal[True]) -> list[_ntp.scputimes]: ...
@overload
def cpu_percent(interval: float | None = None, percpu: Literal[False] = False) -> float: ...
@overload
def cpu_percent(interval: float | None, percpu: Literal[True]) -> list[float]: ...
@overload
def cpu_percent(*, percpu: Literal[True]) -> list[float]: ...
@overload
def cpu_times_percent(interval: float | None = None, percpu: Literal[False] = False) -> _ntp.scputimes: ...
@overload
def cpu_times_percent(interval: float | None, percpu: Literal[True]) -> list[_ntp.scputimes]: ...
@overload
def cpu_times_percent(*, percpu: Literal[True]) -> list[_ntp.scputimes]: ...
def cpu_stats() -> _ntp.scpustats: ...
def getloadavg() -> tuple[float, float, float]: ...
def virtual_memory() -> _ntp.svmem: ...
def swap_memory() -> _ntp.sswap: ...
def disk_usage(path: str) -> _ntp.sdiskusage: ...
def disk_partitions(all: bool = False) -> list[_ntp.sdiskpart]: ...

# TODO: Incorrect sdiskio for BSD systems:
@overload
def disk_io_counters(perdisk: Literal[False] = False, nowrap: bool = True) -> _ntp.sdiskio | None: ...
@overload
def disk_io_counters(perdisk: Literal[True], nowrap: bool = True) -> dict[str, _ntp.sdiskio]: ...
@overload
def net_io_counters(pernic: Literal[False] = False, nowrap: bool = True) -> _ntp.snetio: ...
@overload
def net_io_counters(pernic: Literal[True], nowrap: bool = True) -> dict[str, _ntp.snetio]: ...
def net_connections(kind: str = "inet") -> list[_ntp.sconn]: ...
def net_if_addrs() -> dict[str, list[_ntp.snicaddr]]: ...
def net_if_stats() -> dict[str, _ntp.snicstats]: ...
def boot_time() -> float: ...
def users() -> list[_ntp.suser]: ...
