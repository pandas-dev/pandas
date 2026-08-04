"""
Detection of the CPU allocation actually available to this process.

Used to bound the default worker count for parallel I/O, so that a default
:func:`~pandas.read_csv` does not oversubscribe a machine on which the process
has been given only a slice of the CPUs -- a CPU affinity mask, or a cgroup CPU
quota as set by ``docker --cpus`` and Kubernetes CPU limits.
"""

from __future__ import annotations

import functools
import os


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
