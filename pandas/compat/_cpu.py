"""
Detection of a CPU's physical core count, and of the CPU allocation actually
available to this process.

Used to pick the default worker count for parallel I/O.  The target is the
number of *physical cores*, efficiency cores included: with the work-queued
parallel read path, efficiency cores contribute real throughput (a slow chunk
just means that worker pulls fewer chunks from the queue).  SMT siblings are
excluded because a hyperthread does not add memory bandwidth to the
bandwidth-bound parsing work its sibling is already doing.

That count is then bounded by the allocation this process actually has, so
that a default :func:`~pandas.read_csv` does not oversubscribe a machine on
which the process has been given only a slice of the CPUs -- a CPU affinity
mask, or a cgroup CPU quota as set by ``docker --cpus`` and Kubernetes CPU
limits.
"""

from __future__ import annotations

import ctypes
import functools
import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

# Widest SMT any mainstream part implements (POWER's SMT8).  Used only as a
# sanity bound on a reported topology, not as a hardware assumption.
_MAX_SMT_THREADS = 8


def _parse_cpu_list(spec: str) -> list[int]:
    """
    Parse a Linux ``sysfs`` cpulist string into the CPU ids it names.

    Handles the comma-and-range format, e.g. ``"0-5"`` -> ``[0, 1, 2, 3, 4, 5]``
    and ``"0,2-3"`` -> ``[0, 2, 3]``.
    """
    ids: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            ids.extend(range(int(lo), int(hi) + 1))
        else:
            ids.append(int(part))
    return ids


def _count_distinct_cores(topology: Iterable[tuple[int, int] | None]) -> int | None:
    """
    Count distinct physical cores from ``(package_id, core_id)`` pairs.

    ``None`` entries (unreadable topology) are ignored; returns ``None`` if no
    pair could be read.
    """
    cores = {pair for pair in topology if pair is not None}
    return len(cores) or None


def _count_processor_core_records(
    buf: ctypes.Array[ctypes.c_byte], length: int
) -> int | None:
    """
    Count physical cores in a ``GetLogicalProcessorInformationEx`` buffer.

    ``buf`` holds a sequence of ``SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX``
    records, one per physical core (SMT siblings share a record); ``Size`` is
    a ``DWORD`` at offset ``+4``.
    """
    addr = ctypes.addressof(buf)
    offset = 0
    count = 0
    # offset + 8, not offset: Size is a DWORD at +4, so a record header needs
    # 8 readable bytes.  Bounding on offset alone over-reads a truncated
    # buffer.  The count is taken after the size check so a malformed record
    # is not counted as a core.
    while offset + 8 <= length:
        size = ctypes.c_uint32.from_address(addr + offset + 4).value
        if size == 0:
            break
        count += 1
        offset += size
    return count or None


def _sysctl_int(name: bytes) -> int | None:
    """Read an integer ``sysctl`` by name on macOS, or None on failure."""
    try:
        libc = ctypes.CDLL("libc.dylib", use_errno=True)
        val = ctypes.c_int(0)
        size = ctypes.c_size_t(ctypes.sizeof(val))
        rc = libc.sysctlbyname(name, ctypes.byref(val), ctypes.byref(size), None, 0)
        if rc == 0 and val.value > 0:
            return val.value
    except (OSError, ValueError):
        pass
    return None


def _physical_cores_darwin() -> int | None:
    """Physical core count on macOS (performance and efficiency clusters)."""
    return _sysctl_int(b"hw.physicalcpu")


def _read_sysfs_int(path: str) -> int | None:
    try:
        with open(path, encoding="utf-8") as fh:
            return int(fh.read().strip())
    except (OSError, ValueError):
        return None


def _read_sysfs_str(path: str) -> str | None:
    # ValueError as well as OSError, matching _read_sysfs_int: a non-UTF-8 byte
    # raises UnicodeDecodeError, and this feeds available_cpu_count(), which
    # _default_n_workers calls unguarded -- so it would surface out of read_csv
    # rather than falling back.
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.read().strip()
    except (OSError, ValueError):
        return None


def _cpu_topology(cpu: int) -> tuple[int, int] | None:
    """``(physical_package_id, core_id)`` for a logical CPU, or None."""
    base = f"/sys/devices/system/cpu/cpu{cpu}/topology"
    pkg = _read_sysfs_int(f"{base}/physical_package_id")
    core = _read_sysfs_int(f"{base}/core_id")
    if pkg is None or core is None:
        return None
    return (pkg, core)


def _physical_cores_linux() -> int | None:
    """
    Physical core count on Linux (SMT siblings collapsed).

    Counted over the CPUs the process may actually run on rather than every
    CPU present.  Counting all present cores would let the availability clamp
    in :func:`available_cpu_count` -- a *logical* count -- decide the worker
    count under an affinity mask, which puts an SMT sibling's worth of extra
    threads on every usable core.
    """
    logical = None
    if hasattr(os, "sched_getaffinity"):
        try:
            logical = sorted(os.sched_getaffinity(0))
        except OSError:
            logical = None
    if not logical:
        spec = _read_sysfs_str("/sys/devices/system/cpu/present")
        logical = _parse_cpu_list(spec) if spec else None
    if not logical:
        return None
    # Topology unreadable -> None, leaving physical_core_count()'s
    # os.cpu_count() as the single fallback path.
    physical = _count_distinct_cores(_cpu_topology(cpu) for cpu in logical)
    if physical is not None and len(logical) > _MAX_SMT_THREADS * physical:
        # Degenerate topology: some hypervisors report physical_package_id=0
        # and core_id=0 for every CPU, which collapses to 1 and would make a
        # 12-CPU guest read serially -- slower than main's flat 4.  No
        # mainstream part exceeds 8 threads per core, so a wider ratio than
        # that means the topology is not describing this machine.
        return None
    return physical


def _physical_cores_windows() -> int | None:
    """
    Physical core count on Windows, or None.

    Not reached by the ``read_csv`` default today: parallel reading is gated
    off on Windows (GH#64347), so ``_default_n_workers`` returns 1 before
    consulting this.  It is kept anyway because ``physical_core_count`` is a
    general ``pandas.compat`` helper whose contract is to answer on every
    platform, and because the Windows gate is a measured, temporary
    limitation rather than a permanent one.
    """
    if sys.platform != "win32":
        return None
    # GetLogicalProcessorInformationEx(RelationProcessorCore) reports one
    # record per physical core (SMT siblings share a record).
    try:
        from ctypes import wintypes

        relation_processor_core = 0
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        get_info = kernel32.GetLogicalProcessorInformationEx
        length = wintypes.DWORD(0)
        get_info(relation_processor_core, None, ctypes.byref(length))
        buf = (ctypes.c_byte * length.value)()
        if not get_info(relation_processor_core, buf, ctypes.byref(length)):
            return None
        return _count_processor_core_records(buf, length.value)
    except Exception:
        return None


def _parse_cgroup_v2_quota(text: str) -> float | None:
    """
    CPU count from a cgroup v2 ``cpu.max`` value.

    The format is ``"<quota> <period>"`` (both in microseconds) or ``"max"`` when
    the CPU time is unlimited; returns ``None`` when unlimited or malformed.
    """
    parts = text.split()
    if not parts or parts[0] == "max":
        return None
    try:
        quota = int(parts[0])
        period = int(parts[1]) if len(parts) > 1 else 100000
    except ValueError:
        return None
    if quota > 0 and period > 0:
        return quota / period
    return None


def _cgroup_cpu_quota() -> float | None:
    """CPUs allowed by the process's cgroup CFS quota, or None if unlimited."""
    # cgroup v2; under a cgroup namespace (the common container case) this is the
    # limit the container itself sees.
    text = _read_sysfs_str("/sys/fs/cgroup/cpu.max")
    if text is not None:
        return _parse_cgroup_v2_quota(text)
    # cgroup v1
    quota = _read_sysfs_int("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    period = _read_sysfs_int("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if quota and quota > 0 and period and period > 0:
        return quota / period
    return None


@functools.lru_cache(maxsize=1)
def available_cpu_count() -> int | None:
    """
    Number of CPUs the process may actually use, or None if unconstrained.

    Takes the tighter of the CPU affinity mask (``taskset``, SLURM cpusets) and
    the cgroup CFS quota (Docker ``--cpus``, Kubernetes CPU limits), so a default
    worker count does not oversubscribe a restricted environment.  Returns
    ``None`` only where neither limit is probed -- in practice macOS and
    Windows.  macOS exposes no affinity mask at all; Windows does (via
    ``GetProcessAffinityMask`` and job-object rate limits) but it is not
    implemented here, since parallel reading is off by default there anyway.
    On Linux ``sched_getaffinity`` always succeeds, so an unconstrained box
    reports its full logical CPU count rather than ``None``.

    Cached for the process lifetime: the CPU allocation is static in practice
    (affinity is normally set once at startup, and a live quota change via
    ``docker update`` or an in-place pod resize is rare), and this is on the
    ``read_csv`` hot path.
    """
    limits = []
    if hasattr(os, "sched_getaffinity"):
        try:
            limits.append(len(os.sched_getaffinity(0)))
        except OSError:
            pass
    quota = _cgroup_cpu_quota()
    if quota is not None:
        limits.append(max(1, int(quota)))
    return min(limits) if limits else None


@functools.lru_cache(maxsize=1)
def physical_core_count() -> int:
    """
    Return the number of physical cores (SMT siblings excluded).

    Efficiency cores count: they contribute real throughput to the
    work-queued parallel read.  Whenever the platform probe is unavailable or
    fails, this falls back to :func:`os.cpu_count`.

    Not uniformly a machine property.  On Linux the count is scoped to the
    process's CPU affinity mask, so ``taskset -c 0-7`` on a 16-physical box
    reports 4, not 16; macOS and Windows report the whole machine, and so
    does the ``os.cpu_count`` fallback on every platform.  Callers wanting a
    hardware inventory rather than "cores this process may use" should not
    use this.

    The result is cached for the lifetime of the process.  As with
    :func:`available_cpu_count`, that makes it stale if affinity changes
    after the first call -- fine for a startup-time allocation, wrong for a
    process that re-pins itself.
    """
    total = os.cpu_count() or 1
    probe = {
        "darwin": _physical_cores_darwin,
        "linux": _physical_cores_linux,
        "win32": _physical_cores_windows,
    }.get(sys.platform)
    try:
        n_physical = probe() if probe is not None else None
    except Exception:
        n_physical = None
    if n_physical is None or n_physical < 1 or n_physical > total:
        return total
    return n_physical
