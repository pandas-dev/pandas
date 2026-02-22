from __future__ import annotations

import builtins
import io
import os
import subprocess
import sys
from contextlib import suppress

import pytest

from dask._compatibility import LINUX, MACOS
from dask.system import CPU_COUNT, cpu_count


def test_cpu_count():
    count = cpu_count()
    assert isinstance(count, int)
    assert count == CPU_COUNT
    assert count <= (os.cpu_count() or 999)
    assert count >= 1


@pytest.mark.skipif(MACOS, reason="No CPU affinity in psutil")
@pytest.mark.parametrize(("affinity"), [{0}, {1}, {0, 1}, {0, 2}])
def test_cpu_affinity_psutil(affinity):
    """Test that cpu_count() respects CPU affinity set by psutil"""
    psutil = pytest.importorskip("psutil")
    proc = psutil.Process()
    prev = proc.cpu_affinity()
    if prev is None:
        pytest.skip("No support for CPU affinity")
    if not affinity.issubset(set(prev)):
        pytest.skip("Not enough CPUs")  # pragma: no cover

    proc.cpu_affinity(affinity)
    try:
        assert cpu_count() == len(affinity)
    finally:
        proc.cpu_affinity(prev)


@pytest.mark.skipif(not LINUX, reason="No command line API for CPU affinity")
@pytest.mark.parametrize(("affinity"), [{0}, {1}, {0, 1}, {0, 2}])
def test_cpu_affinity_taskset(affinity):
    """Test that cpu_count() respects the taskset command line tool"""
    count = len(affinity)
    if (os.cpu_count() or 1) < count:
        raise pytest.skip("Not enough CPUs")  # pragma: no cover

    subprocess.check_call(
        [
            "taskset",
            "-c",
            ",".join(str(i) for i in sorted(affinity)),
            sys.executable,
            "-c",
            f"from dask.system import CPU_COUNT; assert CPU_COUNT == {count}",
        ]
    )


@pytest.mark.skipif(
    sys.implementation.name != "cpython" or sys.version_info < (3, 13),
    reason="-X cpu_count= added in CPython 3.13",
)
@pytest.mark.parametrize("count", [1, 2, 3])
def test_cpu_count_arg(count):
    """Test that cpu_count() respects the python -X cpu_count= parameter"""
    if (os.cpu_count() or 1) < count:
        raise pytest.skip("Not enough CPUs")  # pragma: no cover

    subprocess.check_call(
        [
            sys.executable,
            "-X",
            f"cpu_count={count}",
            "-c",
            f"from dask.system import CPU_COUNT; assert CPU_COUNT == {count}",
        ]
    )


@pytest.fixture
def monkeypatch_cpu_count(monkeypatch):
    def cpu_count():
        return 250  # Absurdly high, unlikely to match real value

    def sched_getaffinity(pid):
        return set(range(250))

    class Process:
        def cpu_affinity(self):
            return list(range(250))

    monkeypatch.setattr(os, "cpu_count", cpu_count)
    if sys.version_info >= (3, 13):
        monkeypatch.setattr(os, "process_cpu_count", cpu_count)
    monkeypatch.setattr(os, "sched_getaffinity", sched_getaffinity)
    with suppress(ImportError):
        import psutil

        monkeypatch.setattr(psutil, "Process", Process)


@pytest.mark.skipif(not LINUX, reason="Control Groups only available on Linux")
@pytest.mark.parametrize("dirname", ["cpuacct,cpu", "cpu,cpuacct", None])
def test_cpu_count_cgroups(dirname, monkeypatch, monkeypatch_cpu_count):
    if dirname:
        paths = {
            f"/sys/fs/cgroup/{dirname}/cpu.cfs_quota_us": io.StringIO("2005"),
            f"/sys/fs/cgroup/{dirname}/cpu.cfs_period_us": io.StringIO("10"),
        }
        builtin_open = builtins.open

        def myopen(path, *args, **kwargs):
            if path in paths:
                return paths.get(path)
            return builtin_open(path, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", myopen)
        monkeypatch.setattr(sys, "platform", "linux")

    count = cpu_count()
    if dirname:
        # Rounds up
        assert count == 201
    else:
        assert count == 250


@pytest.mark.skipif(not LINUX, reason="Control Groups only available on Linux")
@pytest.mark.parametrize("group_name", ["/", "/user.slice", "/user.slice/more.slice"])
@pytest.mark.parametrize("quota", ["max", "2005"])
def test_cpu_count_cgroups_v2(quota, group_name, monkeypatch, monkeypatch_cpu_count):
    if not group_name.endswith("/"):
        group_name = f"{group_name}/"

    paths = {
        "/proc/self/cgroup": io.StringIO(f"0::{group_name}"),
        f"/sys/fs/cgroup{group_name}cpu.max": io.StringIO(f"{quota} 10"),
    }
    builtin_open = builtins.open

    def myopen(path, *args, **kwargs):
        if path in paths:
            return paths.get(path)
        return builtin_open(path, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", myopen)
    monkeypatch.setattr(sys, "platform", "linux")

    count = cpu_count()
    if quota == "max":
        assert count == 250
    else:
        # Rounds up
        assert count == 201
