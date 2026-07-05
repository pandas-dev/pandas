"""
_cpu
====

Detection of a CPU's physical performance-core count.

Used to pick the default worker count for parallel I/O.  The target is the
number of *physical performance cores*: on a heterogeneous CPU (performance /
efficiency clusters, e.g. Apple silicon or recent Intel client parts) that is
the fast cluster, and on a homogeneous CPU it is simply the physical core count.
Efficiency cores and SMT siblings are excluded because they do not speed up the
memory-bandwidth-bound CSV parsing work (and an efficiency-core straggler slows
the parallel gather down).
"""

from __future__ import annotations

import ctypes
import functools
import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


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


def _count_top_efficiency_class(
    buf: ctypes.Array[ctypes.c_byte], length: int
) -> int | None:
    """
    Count performance cores in a ``GetLogicalProcessorInformationEx`` buffer.

    ``buf`` holds a sequence of ``SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX``
    records (one per physical core): ``Size`` is a ``DWORD`` at offset ``+4`` and
    ``PROCESSOR_RELATIONSHIP.EfficiencyClass`` is a ``BYTE`` at offset ``+9``.
    The highest efficiency class is the performance cores.
    """
    addr = ctypes.addressof(buf)
    offset = 0
    eff_classes = []
    while offset < length:
        size = ctypes.c_uint32.from_address(addr + offset + 4).value
        eff_classes.append(ctypes.c_uint8.from_address(addr + offset + 9).value)
        if size == 0:
            break
        offset += size
    if not eff_classes:
        return None
    top = max(eff_classes)
    return sum(1 for cls in eff_classes if cls == top)


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


def _perf_cores_darwin() -> int | None:
    """Physical performance cores on macOS (P-cluster on hybrid, else all)."""
    # perflevel0 is the highest-performance cluster on Apple silicon; on a
    # homogeneous CPU (e.g. an Intel Mac) it is absent, so fall back to the
    # overall physical core count.
    for name in (b"hw.perflevel0.physicalcpu", b"hw.physicalcpu"):
        val = _sysctl_int(name)
        if val is not None:
            return val
    return None


def _read_sysfs_int(path: str) -> int | None:
    try:
        with open(path, encoding="utf-8") as fh:
            return int(fh.read().strip())
    except (OSError, ValueError):
        return None


def _read_sysfs_str(path: str) -> str | None:
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.read().strip()
    except OSError:
        return None


def _cpu_topology(cpu: int) -> tuple[int, int] | None:
    """``(physical_package_id, core_id)`` for a logical CPU, or None."""
    base = f"/sys/devices/system/cpu/cpu{cpu}/topology"
    pkg = _read_sysfs_int(f"{base}/physical_package_id")
    core = _read_sysfs_int(f"{base}/core_id")
    if pkg is None or core is None:
        return None
    return (pkg, core)


def _perf_logical_cpus_linux() -> list[int]:
    """Logical CPU ids that belong to the performance cluster on Linux."""
    # Intel hybrid parts expose the performance cores via a dedicated driver dir.
    spec = _read_sysfs_str("/sys/devices/cpu_core/cpus")
    if spec:
        ids = _parse_cpu_list(spec)
        if ids:
            return ids
    # Homogeneous (or non-Intel): every present CPU is a performance core.
    spec = _read_sysfs_str("/sys/devices/system/cpu/present")
    if spec:
        return _parse_cpu_list(spec)
    return []


def _perf_cores_linux() -> int | None:
    """Physical performance cores on Linux (SMT siblings collapsed)."""
    logical = _perf_logical_cpus_linux()
    if not logical:
        return None
    physical = _count_distinct_cores(_cpu_topology(cpu) for cpu in logical)
    if physical is not None:
        return physical
    # Topology unreadable: fall back to the logical performance-CPU count.
    return len(logical)


def _perf_cores_windows() -> int | None:
    """Physical performance cores on Windows, or None."""
    if sys.platform != "win32":
        return None
    # GetLogicalProcessorInformationEx(RelationProcessorCore) reports an
    # EfficiencyClass per physical core; the highest class is the performance
    # cores.
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
        return _count_top_efficiency_class(buf, length.value)
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
    ``None`` on platforms/setups without either limit (e.g. macOS, Windows, or an
    unconstrained Linux box).

    Cached for the process lifetime: the CPU allocation is effectively static
    (cgroup limits do not change at runtime and affinity is normally set once at
    startup), and this is on the ``read_csv`` hot path.
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
def performance_core_count() -> int:
    """
    Return the number of physical performance cores available to the process.

    On heterogeneous CPUs (e.g. Apple silicon or recent Intel client parts with
    performance/efficiency clusters) this is the size of the fast cluster.  On
    homogeneous CPUs it is the physical core count.  Whenever the platform probe
    is unavailable or fails, this falls back to :func:`os.cpu_count`.

    The result is cached for the lifetime of the process.
    """
    total = os.cpu_count() or 1
    probe = {
        "darwin": _perf_cores_darwin,
        "linux": _perf_cores_linux,
        "win32": _perf_cores_windows,
    }.get(sys.platform)
    try:
        n_perf = probe() if probe is not None else None
    except Exception:
        n_perf = None
    if n_perf is None or n_perf < 1 or n_perf > total:
        return total
    return n_perf
