import sys

# sys.platform.startswith(("sunos", "solaris")):
if sys.platform != "linux" and sys.platform != "win32" and sys.platform != "darwin":
    from _typeshed import Incomplete
    from collections.abc import Callable
    from typing import Final, Literal, NamedTuple, TypeVar, overload
    from typing_extensions import ParamSpec

    from psutil._common import (
        AF_INET6 as AF_INET6,
        ENCODING as ENCODING,
        AccessDenied as AccessDenied,
        NoSuchProcess as NoSuchProcess,
        ZombieProcess as ZombieProcess,
        debug as debug,
        get_procfs_path as get_procfs_path,
        isfile_strict as isfile_strict,
        memoize_when_activated as memoize_when_activated,
        sockfam_to_enum as sockfam_to_enum,
        socktype_to_enum as socktype_to_enum,
        usage_percent as usage_percent,
    )

    from . import _ntuples as ntp, _psposix, _psutil_sunos

    _P = ParamSpec("_P")
    _R = TypeVar("_R")

    __extra__all__: Final[list[str]]
    PAGE_SIZE: Final[int]
    AF_LINK: Final[int]
    IS_64_BIT: Final[bool]
    CONN_IDLE: Final = "IDLE"
    CONN_BOUND: Final = "BOUND"
    PROC_STATUSES: Final[dict[int, str]]
    TCP_STATUSES: Final[dict[int, str]]
    proc_info_map: Final[dict[str, int]]

    class scputimes(NamedTuple):
        user: float
        system: float
        idle: float
        iowait: float

    class pcputimes(NamedTuple):
        user: float
        system: float
        children_user: float
        children_system: float

    class svmem(NamedTuple):
        total: int
        available: int
        percent: float
        used: int
        free: int

    class pmem(NamedTuple):
        rss: int
        vms: int

    pfullmem = pmem

    class pmmap_grouped(NamedTuple):
        path: Incomplete
        rss: Incomplete
        anonymous: Incomplete
        locked: Incomplete

    class pmmap_ext(NamedTuple):
        addr: Incomplete
        perms: Incomplete
        path: Incomplete
        rss: Incomplete
        anonymous: Incomplete
        locked: Incomplete

    def virtual_memory() -> svmem: ...
    def swap_memory() -> ntp.sswap: ...
    def cpu_times() -> scputimes: ...
    def per_cpu_times() -> list[scputimes]: ...
    def cpu_count_logical() -> int | None: ...
    def cpu_count_cores() -> int | None: ...
    def cpu_stats() -> ntp.scpustats: ...

    disk_io_counters = _psutil_sunos.disk_io_counters
    disk_usage = _psposix.disk_usage

    def disk_partitions(all: bool = False) -> list[ntp.sdiskpart]: ...

    net_io_counters = _psutil_sunos.net_io_counters
    net_if_addrs = _psutil_sunos.net_if_addrs

    @overload
    def net_connections(kind: str, _pid: Literal[-1] = -1) -> list[ntp.sconn]: ...
    @overload
    def net_connections(kind: str, _pid: int = -1) -> list[ntp.pconn]: ...
    def net_if_stats() -> dict[str, ntp.snicstats]: ...
    def boot_time() -> float: ...
    def users() -> list[ntp.suser]: ...
    def pids() -> list[int]: ...
    def pid_exists(pid: int) -> bool: ...
    def wrap_exceptions(fun: Callable[_P, _R]) -> Callable[_P, _R]: ...

    class Process:
        __slots__ = ["_cache", "_name", "_ppid", "_procfs_path", "pid"]
        pid: int
        def __init__(self, pid: int) -> None: ...
        def oneshot_enter(self) -> None: ...
        def oneshot_exit(self) -> None: ...
        def name(self) -> str: ...
        def exe(self) -> str: ...
        def cmdline(self) -> list[str] | None: ...
        def environ(self) -> dict[str, str]: ...
        def create_time(self) -> float: ...
        def num_threads(self) -> int: ...
        def nice_get(self) -> int: ...
        def nice_set(self, value: int) -> None: ...
        def ppid(self) -> int: ...
        def uids(self) -> ntp.puids: ...
        def gids(self) -> ntp.puids: ...
        def cpu_times(self) -> ntp.pcputimes: ...
        def cpu_num(self) -> int: ...
        def terminal(self) -> str | None: ...
        def cwd(self) -> str: ...
        def memory_info(self) -> pmem: ...
        memory_full_info = memory_info
        def status(self) -> str: ...
        def threads(self) -> list[ntp.pthread]: ...
        def open_files(self) -> list[ntp.popenfile]: ...
        def net_connections(self, kind: str = "inet") -> list[ntp.pconn]: ...

        class nt_mmap_grouped(NamedTuple):
            path: Incomplete
            rss: Incomplete
            anon: Incomplete
            locked: Incomplete

        class nt_mmap_ext(NamedTuple):
            addr: Incomplete
            perms: Incomplete
            path: Incomplete
            rss: Incomplete
            anon: Incomplete
            locked: Incomplete

        def memory_maps(self) -> list[tuple[str, str, str, int, int, int]]: ...
        def num_fds(self) -> int: ...
        def num_ctx_switches(self) -> ntp.pctxsw: ...
        def wait(self, timeout: float | None = None) -> int | None: ...
