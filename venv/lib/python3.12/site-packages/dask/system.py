from __future__ import annotations

import math
import os
import sys

try:
    import psutil
except ImportError:
    psutil = None  # type: ignore

__all__ = ("cpu_count", "CPU_COUNT")


def _try_extract_cgroup_cpu_quota():
    # cgroup v1
    # The directory name isn't standardized across linux distros, check both
    for dirname in ["cpuacct,cpu", "cpu,cpuacct"]:
        try:
            with open("/sys/fs/cgroup/%s/cpu.cfs_quota_us" % dirname) as f:
                quota = int(f.read())
            with open("/sys/fs/cgroup/%s/cpu.cfs_period_us" % dirname) as f:
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
        with open("/sys/fs/cgroup%scpu.max" % group_path) as f:
            quota, period = map(int, f.read().split(" "))
            return quota, period
    except Exception:
        pass

    # No cgroup CPU quota found
    return None, None


def cpu_count():
    """Get the available CPU count for this system.

    Takes the minimum value from the following locations:

    - Total system cpus available on the host.
    - CPU Affinity (if set)
    - Cgroups limit (if set)
    """
    count = os.cpu_count()

    # Check CPU affinity if available
    if psutil is not None:
        try:
            affinity_count = len(psutil.Process().cpu_affinity())
            if affinity_count > 0:
                count = min(count, affinity_count)
        except Exception:
            pass

    # Check cgroups if available
    if sys.platform == "linux":
        quota, period = _try_extract_cgroup_cpu_quota()
        if quota is not None and period is not None:
            # We round up on fractional CPUs
            cgroups_count = math.ceil(quota / period)
            if cgroups_count > 0:
                count = min(count, cgroups_count)

    return count


CPU_COUNT = cpu_count()
