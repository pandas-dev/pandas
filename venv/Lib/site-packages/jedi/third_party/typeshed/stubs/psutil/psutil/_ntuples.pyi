import sys
from _typeshed import Incomplete
from socket import AddressFamily, SocketKind
from typing import Any, NamedTuple

# All named tuples are defined in this file, but due to the inability to detect some platforms,
# it was decided to store the correct named tuples inside platform-specific files.

class sswap(NamedTuple):
    total: int
    used: int
    free: int
    percent: float
    sin: int
    sout: int

class sdiskusage(NamedTuple):
    total: int
    used: int
    free: int
    percent: float

# redefine for linux:
if sys.platform != "linux":
    class sdiskio(NamedTuple):
        read_count: int
        write_count: int
        read_bytes: int
        write_bytes: int
        read_time: int
        write_time: int

class sdiskpart(NamedTuple):
    device: str
    mountpoint: str
    fstype: str
    opts: str

class snetio(NamedTuple):
    bytes_sent: int
    bytes_recv: int
    packets_sent: int
    packets_recv: int
    errin: int
    errout: int
    dropin: int
    dropout: int

class suser(NamedTuple):
    name: str
    terminal: str | None
    host: str | None
    started: float
    pid: str

class sconn(NamedTuple):
    fd: int
    family: AddressFamily
    type: SocketKind
    laddr: addr | tuple[()]
    raddr: addr | tuple[()]
    status: str
    pid: int | None

class snicaddr(NamedTuple):
    family: AddressFamily
    address: str
    netmask: str | None
    broadcast: str | None
    ptp: str | None

class snicstats(NamedTuple):
    isup: bool
    duplex: int
    speed: int
    mtu: int
    flags: str

class scpustats(NamedTuple):
    ctx_switches: int
    interrupts: int
    soft_interrupts: int
    syscalls: int

class scpufreq(NamedTuple):
    current: float
    min: float
    max: float

class shwtemp(NamedTuple):
    label: str
    current: float
    high: float | None
    critical: float | None

class sbattery(NamedTuple):
    percent: int
    secsleft: int
    power_plugged: bool

class sfan(NamedTuple):
    label: str
    current: int

if sys.platform == "win32":
    class pheap(NamedTuple):
        heap_used: Incomplete
        mmap_used: Incomplete
        heap_count: Incomplete

else:
    # if LINUX or MACOS or BSD:
    class pheap(NamedTuple):
        heap_used: Incomplete
        mmap_used: Incomplete

# redefine for linux:
if sys.platform != "linux":
    class pcputimes(NamedTuple):
        user: float
        system: float
        children_user: float
        children_system: float

    class popenfile(NamedTuple):
        path: str
        fd: int

class pthread(NamedTuple):
    id: int
    user_time: float
    system_time: float

class puids(NamedTuple):
    real: int
    effective: int
    saved: int

class pgids(NamedTuple):
    real: int
    effective: int
    saved: int

# redefine for linux and windows:
if sys.platform != "linux" and sys.platform != "win32":
    class pio(NamedTuple):
        read_count: int
        write_count: int
        read_bytes: int
        write_bytes: int

class pionice(NamedTuple):
    ioclass: int
    value: int

class pctxsw(NamedTuple):
    voluntary: int
    involuntary: int

class pconn(NamedTuple):
    fd: int
    family: AddressFamily
    type: SocketKind
    laddr: addr
    raddr: addr
    status: str

class addr(NamedTuple):
    ip: str
    port: int

if sys.platform == "linux":
    class scputimes(NamedTuple):
        # Note: scputimes has different fields depending on exactly how Linux
        # is setup, but we'll include the "complete" set of fields
        user: float
        nice: float
        system: float
        idle: float
        iowait: float
        irq: float
        softirq: float
        steal: float
        guest: float
        guest_nice: float

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
        slab: int

    class sdiskio(NamedTuple):
        read_count: int
        write_count: int
        read_bytes: int
        write_bytes: int
        read_time: int
        write_time: int
        read_merged_count: int
        write_merged_count: int
        busy_time: int

    class popenfile(NamedTuple):
        path: str
        fd: int
        position: int
        mode: str
        flags: int

    class pmem(NamedTuple):
        rss: int
        vms: int
        shared: int
        text: int
        lib: int
        data: int
        dirty: int

    class pfullmem(NamedTuple):
        rss: int
        vms: int
        shared: int
        text: int
        lib: int
        data: int
        dirty: int
        uss: int
        pss: int
        swap: int

    class pmmap_grouped(NamedTuple):
        path: Incomplete
        rss: Incomplete
        size: Incomplete
        pss: Incomplete
        shared_clean: Incomplete
        shared_dirty: Incomplete
        private_clean: Incomplete
        private_dirty: Incomplete
        referenced: Incomplete
        anonymous: Incomplete
        swap: Incomplete

    class pmmap_ext(NamedTuple):
        addr: Incomplete
        perms: Incomplete
        path: Incomplete
        rss: Incomplete
        size: Incomplete
        pss: Incomplete
        shared_clean: Incomplete
        shared_dirty: Incomplete
        private_clean: Incomplete
        private_dirty: Incomplete
        referenced: Incomplete
        anonymous: Incomplete
        swap: Incomplete

    class pio(NamedTuple):
        read_count: int
        write_count: int
        read_bytes: int
        write_bytes: int
        read_chars: int
        write_chars: int

    class pcputimes(NamedTuple):
        user: float
        system: float
        children_user: float
        children_system: float
        iowait: float

elif sys.platform == "win32":
    class scputimes(NamedTuple):
        user: float
        system: float
        idle: float
        interrupt: float
        dpc: float

    class svmem(NamedTuple):
        total: int
        available: int
        percent: float
        used: int
        free: int

    class pmem(NamedTuple):
        rss: int
        vms: int
        num_page_faults: int
        peak_wset: int
        wset: int
        peak_paged_pool: int
        paged_pool: int
        peak_nonpaged_pool: int
        nonpaged_pool: int
        pagefile: int
        peak_pagefile: int
        private: int

    class pfullmem(NamedTuple):
        rss: int
        vms: int
        num_page_faults: int
        peak_wset: int
        wset: int
        peak_paged_pool: int
        paged_pool: int
        peak_nonpaged_pool: int
        nonpaged_pool: int
        pagefile: int
        peak_pagefile: int
        private: int
        uss: int

    class pmmap_grouped(NamedTuple):
        path: Incomplete
        rss: Incomplete

    class pmmap_ext(NamedTuple):
        addr: Incomplete
        perms: Incomplete
        path: Incomplete
        rss: Incomplete

    class pio(NamedTuple):
        read_count: int
        write_count: int
        read_bytes: int
        write_bytes: int
        other_count: int
        other_bytes: int

elif sys.platform == "darwin":
    class scputimes(NamedTuple):
        user: float
        nice: float
        system: float
        idle: float

    class svmem(NamedTuple):
        total: int
        available: int
        percent: float
        used: int
        free: int
        active: int
        inactive: int
        wired: int

    class pmem(NamedTuple):
        rss: int
        vms: int
        pfaults: int
        pageins: int

    class pfullmem(NamedTuple):
        rss: int
        vms: int
        pfaults: int
        pageins: int
        uss: int

else:
    # See _psbsd.pyi, _pssunos.pyi or _psaix.pyi
    # BSD: svmem, scputimes, pmem, pfullmem, pcputimes, pmmap_grouped, pmmap_ext, sdiskio
    # SUNOS: scputimes, pcputimes, svmem, pmem, pfullmem, pmmap_grouped, pmmap_ext
    # AIX: pmem, pfullmem, scputimes, svmem

    scputimes = Incomplete

    class pmem(Any): ...
    class pfullmem(Any): ...
    class svmem(Any): ...
