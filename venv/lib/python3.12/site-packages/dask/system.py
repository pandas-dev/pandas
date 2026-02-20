from __future__ import annotations

import math
import os
import sys

try:
    import psutil
except ImportError:
    psutil = None  # type: ignore[assignment]

from dask._compatibility import LINUX

__all__ = ("cpu_count", "CPU_COUNT")


def _try_extract_cgroup_cpu_quota():
    # cgroup v1
    # The directory name isn't standardized across linux distros, check both
    for dirname in ["cpuacct,cpu", "cpu,cpuacct"]:
        try:
            with open(f"/sys/fs/cgroup/{dirname}/cpu.cfs_quota_us") as f:
                quota = int(f.read())
            with open(f"/sys/fs/cgroup/{dirname}/cpu.cfs_period_us") as f:
                period = int(f.read())
            return quota, period
        except Exception:
            pass

    # cgroup v2
    try:
        with open("/proc/self/cgroup") as f:
            group_path = f.read().strip().split(":")[-1]
        if not group_path.endswith("/"):
            group_path = f"{group_path}/"
        with open(f"/sys/fs/cgroup{group_path}cpu.max") as f:
            quota, period = map(int, f.read().split(" "))
            return quota, period
    except Exception:
        pass

    # No cgroup CPU quota found
    return None, None


def cpu_count() -> int:
    """Get the available CPU count for this system.

    Takes the minimum value from the following locations:

    - Total system cpus available on the host.
    - CPU Affinity (if set)
    - Cgroups limit (if set)
    """
    if sys.version_info >= (3, 13):
        # Embeds CPU affinity checks
        count = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):
        # https://docs.python.org/3/library/os.html#interface-to-the-scheduler
        # "only available on some Unix platforms"; neither MacOS nor Windows
        count = len(os.sched_getaffinity(0))
    else:
        # Does not account for CPU affinity.
        # On exotic alternative Python implementations, it may return None.
        count = os.cpu_count() or 1
    assert count

    # Additional CPU affinity check with psutil.
    # NOTE: do not limit this to Python <3.13: on Windows,
    # `psutil.Process().cpu_affinity(value)` does not change the reading of
    # os.process_cpu_count().
    if psutil is not None:
        proc = psutil.Process()
        if hasattr(proc, "cpu_affinity"):
            affinity = proc.cpu_affinity()
            if affinity is not None:
                assert affinity
                count = min(count, len(affinity))

    # Check cgroups if available
    if LINUX:
        quota, period = _try_extract_cgroup_cpu_quota()
        if quota is not None and period is not None:
            # We round up on fractional CPUs
            cgroups_count = math.ceil(quota / period)
            if cgroups_count > 0:
                count = min(count, cgroups_count)

    return count


CPU_COUNT = cpu_count()
