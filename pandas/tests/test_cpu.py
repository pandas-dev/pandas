"""
Tests for pandas.compat._cpu: detection of the CPU allocation available to
this process (affinity mask and cgroup CFS quota), and the default thread count
derived from it.
"""

import os

import pytest

from pandas.compat import _cpu
from pandas.compat._cpu import (
    _parse_cgroup_v2_quota,
    available_cpu_count,
    get_max_threads,
)

import pandas as pd


@pytest.mark.parametrize(
    "text, expected",
    [
        ("max 100000", None),  # unlimited
        ("max", None),
        ("400000 100000", 4.0),
        ("150000 100000", 1.5),
        ("50000 100000", 0.5),  # docker --cpus=0.5, Kubernetes cpu: 500m
        ("100000", 1.0),  # period defaults to 100000us
        ("0 100000", None),  # zero quota -> ignored
        ("garbage", None),
        ("", None),
    ],
)
def test_parse_cgroup_v2_quota(text, expected):
    assert _parse_cgroup_v2_quota(text) == expected


@pytest.mark.parametrize(
    "affinity, quota, expected",
    [
        ({0, 1, 2, 3, 4, 5, 6, 7}, None, 8),  # affinity only
        ({0, 1, 2, 3, 4, 5, 6, 7}, 4.0, 4),  # cgroup quota tighter
        ({0, 1, 2, 3}, 8.0, 4),  # affinity tighter
        ({0, 1}, 1.5, 1),  # fractional quota floors, min 1
        ({0, 1, 2, 3}, 0.5, 1),  # sub-1-CPU quota still leaves one worker
    ],
)
def test_available_cpu_count(monkeypatch, affinity, quota, expected):
    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: set(affinity), raising=False
    )
    monkeypatch.setattr(_cpu, "_cgroup_cpu_quota", lambda: quota)
    available_cpu_count.cache_clear()
    try:
        assert available_cpu_count() == expected
    finally:
        available_cpu_count.cache_clear()


def test_available_cpu_count_unconstrained(monkeypatch):
    # No affinity limit and no cgroup quota -> None (do not clamp).

    def raise_oserror(pid):
        raise OSError

    monkeypatch.setattr(os, "sched_getaffinity", raise_oserror, raising=False)
    monkeypatch.setattr(_cpu, "_cgroup_cpu_quota", lambda: None)
    available_cpu_count.cache_clear()
    try:
        assert available_cpu_count() is None
    finally:
        available_cpu_count.cache_clear()


def test_available_cpu_count_does_not_raise():
    # get_max_threads calls available_cpu_count() bare, so anything escaping
    # it escapes read_csv.  Every other test mocks both the affinity call and
    # the cgroup reads; this one runs it against the host.
    available_cpu_count.cache_clear()
    try:
        result = available_cpu_count()
    finally:
        available_cpu_count.cache_clear()
    assert result is None or (isinstance(result, int) and result >= 1)


def test_cgroup_cpu_quota_v2(monkeypatch):
    monkeypatch.setattr(
        _cpu,
        "_read_sysfs_str",
        lambda path: "400000 100000" if path == "/sys/fs/cgroup/cpu.max" else None,
    )
    assert _cpu._cgroup_cpu_quota() == 4.0


@pytest.mark.parametrize("quota, expected", [(300000, 3.0), (-1, None)])
def test_cgroup_cpu_quota_v1_fallback(monkeypatch, quota, expected):
    # -1 is the unlimited sentinel cpu.cfs_quota_us holds on an unconstrained
    # cgroup v1 host; reading it as a quota would force every read_csv serial.
    ints = {
        "/sys/fs/cgroup/cpu/cpu.cfs_quota_us": quota,
        "/sys/fs/cgroup/cpu/cpu.cfs_period_us": 100000,
    }
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: None)  # no cgroup v2
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: ints.get(path))
    assert _cpu._cgroup_cpu_quota() == expected


def test_cgroup_cpu_quota_unlimited(monkeypatch):
    monkeypatch.setattr(
        _cpu,
        "_read_sysfs_str",
        lambda path: "max 100000" if path == "/sys/fs/cgroup/cpu.max" else None,
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: None)
    assert _cpu._cgroup_cpu_quota() is None


def test_read_sysfs_str_non_utf8(tmp_path):
    # A non-UTF-8 byte must read as "unknown", not raise out of read_csv.
    path = tmp_path / "cpu.max"
    path.write_bytes(b"\xff\xfe")
    assert _cpu._read_sysfs_str(str(path)) is None


def test_read_sysfs_int_non_numeric(tmp_path):
    # Same for a non-numeric value; every other test replaces _read_sysfs_int
    # with a stub, so this is the only check that it swallows ValueError.
    path = tmp_path / "cpu.cfs_quota_us"
    path.write_text("max", encoding="utf-8")
    assert _cpu._read_sysfs_int(str(path)) is None


class TestGetMaxThreads:
    """Unit tests for the default worker count."""

    @pytest.mark.parametrize("platform_name", ["linux", "darwin", "win32"])
    @pytest.mark.parametrize(
        "cpu_count, available, expected",
        [
            (2, None, 2),  # unconstrained, below the cap -> logical CPU count
            (16, None, 6),  # cap binds
            (16, 1, 1),  # single-CPU container
            (16, 2, 2),  # cgroup/affinity tighter than the cap
            (16, 8, 6),  # allocation looser than the cap -> cap still binds
            (2, 8, 2),  # allocation looser than the machine
        ],
    )
    def test_combines_cap_and_allocation(
        self, monkeypatch, cpu_count, available, expected, platform_name
    ):
        # Default = min(logical CPUs, available CPUs, MAX_DEFAULT_THREADS),
        # on every threaded platform.
        assert _cpu.MAX_DEFAULT_THREADS == 6
        monkeypatch.setattr(_cpu.sys, "platform", platform_name)
        monkeypatch.setattr(_cpu.os, "cpu_count", lambda: cpu_count)
        monkeypatch.setattr(_cpu, "available_cpu_count", lambda: available)
        with pd.option_context("mode.max_threads", None):
            assert get_max_threads() == expected

    def test_wasm_is_serial(self, monkeypatch):
        monkeypatch.setattr(_cpu.sys, "platform", "emscripten")
        # WASM stays serial even if a worker count is requested explicitly.
        with pd.option_context("mode.max_threads", 8):
            assert get_max_threads() == 1

    def test_max_threads_exceeds_cap(self, monkeypatch):
        # MAX_DEFAULT_THREADS and the availability clamp bound the *default*
        # only.  The mode.max_threads docs promise an explicit setting still
        # wins.
        monkeypatch.setattr(_cpu.sys, "platform", "linux")
        monkeypatch.setattr(_cpu.os, "cpu_count", lambda: 2)
        monkeypatch.setattr(_cpu, "available_cpu_count", lambda: 1)
        with pd.option_context("mode.max_threads", 32):
            assert get_max_threads() == 32

    @pytest.mark.parametrize("platform_name", ["linux", "win32"])
    def test_max_threads_wins(self, monkeypatch, platform_name):
        monkeypatch.setattr(_cpu.sys, "platform", platform_name)
        with pd.option_context("mode.max_threads", 3):
            assert get_max_threads() == 3
