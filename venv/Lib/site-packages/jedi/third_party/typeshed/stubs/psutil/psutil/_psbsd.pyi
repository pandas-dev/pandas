import sys

# sys.platform.startswith(("freebsd", "midnightbsd", "openbsd", "netbsd")):
if sys.platform != "linux" and sys.platform != "win32" and sys.platform != "darwin":
    from _typeshed import Incomplete
    from collections import defaultdict
    from collections.abc import Callable
    from contextlib import AbstractContextManager
    from typing import Final, NamedTuple, TypeVar, overload
    from typing_extensions import ParamSpec

    from psutil._common import (
        FREEBSD as FREEBSD,
        NETBSD as NETBSD,
        OPENBSD as OPENBSD,
        AccessDenied as AccessDenied,
        NoSuchProcess as NoSuchProcess,
        ZombieProcess as ZombieProcess,
        conn_tmap as conn_tmap,
        conn_to_ntuple as conn_to_ntuple,
        memoize as memoize,
        usage_percent as usage_percent,
    )

    from . import _ntuples as ntp, _psposix, _psutil_bsd

    _P = ParamSpec("_P")
    _R = TypeVar("_R")

    __extra__all__: Final[list[str]]
    PROC_STATUSES: Final[dict[int, str]]
    TCP_STATUSES: Final[dict[int, str]]
    PAGESIZE: Final[int]
    AF_LINK: Final = _psutil_bsd.AF_LINK
    HAS_PROC_NUM_THREADS: Final[bool]
    kinfo_proc_map: Final[dict[str, int]]

    class svmem(NamedTuple):
        total: int
        available: int
        percent: float
        used: int
        free: int
        active: int
        inactive: int
        buffers: int
        cached: int
        shared: int
        wired: int

    class scputimes(NamedTuple):
        user: float
        nice: float
        system: float
        idle: float
        irq: float

    class pmem(NamedTuple):
        rss: int
        vms: int
        text: int
        data: int
        stack: int

    pfullmem = pmem

    class pcputimes(NamedTuple):
        user: float
        system: float
        children_user: float
        children_system: float

    class pmmap_grouped(NamedTuple):
        path: Incomplete
        rss: Incomplete
        private: Incomplete
        ref_count: Incomplete
        shadow_count: Incomplete

    class pmmap_ext(NamedTuple):
        addr: Incomplete
        perms: Incomplete
        path: Incomplete
        rss: Incomplete
        private: Incomplete
        ref_count: Incomplete
        shadow_count: Incomplete

    class sdiskio(NamedTuple):
        read_count: Incomplete
        write_count: Incomplete
        read_bytes: Incomplete
        write_bytes: Incomplete
        read_time: Incomplete
        write_time: Incomplete
        busy_time: Incomplete

    def virtual_memory() -> svmem: ...
    def swap_memory() -> ntp.sswap: ...
    heap_info = _psutil_bsd.heap_info  # only FreeBSD and NetBSD
    heap_trim = _psutil_bsd.heap_trim  # only FreeBSD and NetBSD
    def cpu_times() -> scputimes: ...
    def per_cpu_times() -> list[scputimes]: ...
    def cpu_count_logical() -> int | None: ...
    def cpu_count_cores() -> int | None: ...
    def cpu_stats() -> ntp.scpustats: ...
    def disk_partitions(all: bool = False) -> list[ntp.sdiskpart]: ...

    disk_usage = _psposix.disk_usage
    disk_io_counters = _psutil_bsd.disk_io_counters
    net_io_counters = _psutil_bsd.net_io_counters
    net_if_addrs = _psutil_bsd.net_if_addrs

    def net_if_stats() -> dict[str, ntp.snicstats]: ...
    def net_connections(kind: str) -> list[ntp.sconn]: ...
    def sensors_battery() -> ntp.sbattery | None: ...  # only FreeBSD
    def sensors_temperatures() -> defaultdict[str, list[ntp.shwtemp]]: ...  # only FreeBSD
    def cpu_freq() -> list[ntp.scpufreq]: ...  # only FreeBSD and OpenBSD
    def boot_time() -> float: ...
    def users() -> list[ntp.suser]: ...

    INIT_BOOT_TIME: Final[float]  # only NetBSD

    def adjust_proc_create_time(ctime: float) -> float: ...  # only NetBSD
    def pids() -> list[int]: ...
    def pid_exists(pid: int) -> bool: ...
    def wrap_exceptions(fun: Callable[_P, _R]) -> Callable[_P, _R]: ...
    def wrap_exceptions_procfs(inst: Process) -> AbstractContextManager[None]: ...

    class Process:
        __slots__ = ["_cache", "_name", "_ppid", "pid"]
        pid: int
        def __init__(self, pid: int) -> None: ...
        def oneshot(
            self,
        ) -> tuple[
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            float,
            int,
            int,
            int,
            int,
            float,
            float,
            float,
            float,
            int,
            int,
            int,
            int,
            int,
            int,
            str,
        ]: ...
        def oneshot_enter(self) -> None: ...
        def oneshot_exit(self) -> None: ...
        def name(self) -> str: ...
        def exe(self) -> str: ...
        def cmdline(self) -> list[str]: ...
        def environ(self) -> dict[str, str]: ...
        def terminal(self) -> str | None: ...
        def ppid(self) -> int: ...
        def uids(self) -> ntp.puids: ...
        def gids(self) -> ntp.pgids: ...
        def cpu_times(self) -> ntp.pcputimes: ...
        def cpu_num(self) -> int: ...  # only FreeBSD
        def memory_info(self) -> pmem: ...
        memory_full_info = memory_info
        def create_time(self, monotonic: bool = False) -> float: ...
        def num_threads(self) -> int: ...
        def num_ctx_switches(self) -> ntp.pctxsw: ...
        def threads(self) -> list[ntp.pthread]: ...
        def net_connections(self, kind: str = "inet") -> list[ntp.pconn]: ...
        def wait(self, timeout: float | None = None) -> int | None: ...
        def nice_get(self) -> int: ...
        def nice_set(self, value: int) -> None: ...
        def status(self) -> str: ...
        def io_counters(self) -> ntp.pio: ...
        def cwd(self) -> str: ...

        class nt_mmap_grouped(NamedTuple):
            path: Incomplete
            rss: Incomplete
            private: Incomplete
            ref_count: Incomplete
            shadow_count: Incomplete

        class nt_mmap_ext(NamedTuple):
            addr: Incomplete
            perms: Incomplete
            path: Incomplete
            rss: Incomplete
            private: Incomplete
            ref_count: Incomplete
            shadow_count: Incomplete

        def open_files(self) -> list[ntp.popenfile]: ...
        def num_fds(self) -> int: ...
        def cpu_affinity_get(self) -> list[int]: ...  # only FreeBSD
        def cpu_affinity_set(self, cpus: list[int]) -> None: ...  # only FreeBSD
        def memory_maps(self) -> list[tuple[str, str, str, int, int, int, int]]: ...  # only FreeBSD
        @overload
        def rlimit(self, resource: int, limits: tuple[int, int]) -> None: ...  # only FreeBSD
        @overload
        def rlimit(self, resource: int, limits: None = None) -> tuple[int, int]: ...  # only FreeBSD
