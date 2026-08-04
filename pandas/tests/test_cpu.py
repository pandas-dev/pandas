"""
Tests for pandas.compat._cpu: detection of the physical core count, and of the
CPU allocation available to this process (affinity mask and cgroup CFS quota).
"""

import ctypes
import os
import re
import sys

import pytest

from pandas.compat import _cpu
from pandas.compat._cpu import (
    _count_distinct_cores,
    _count_processor_core_records,
    _parse_cgroup_v2_quota,
    _parse_cpu_list,
    available_cpu_count,
    physical_core_count,
)


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
    # Unlike physical_core_count(), which wraps its probe in `except
    # Exception`, _default_n_workers calls available_cpu_count() bare -- so
    # anything escaping it escapes read_csv.  Every other test mocks both the
    # affinity call and the cgroup reads; this one runs it against the host.
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


# ---------------------------------------------------------------------------
# Physical-core detection (GH#66152)
# ---------------------------------------------------------------------------


def test_physical_core_count_within_bounds():
    # Whatever the host topology, the detector returns a positive int no larger
    # than the logical CPU count, and never raises.
    physical_core_count.cache_clear()
    try:
        n_physical = physical_core_count()
    finally:
        physical_core_count.cache_clear()
    assert isinstance(n_physical, int)
    assert 1 <= n_physical <= (os.cpu_count() or 1)


# The probe tests around this pair are fully mocked, and the bounds assertion
# above is also satisfied when the probe returns None and physical_core_count()
# falls back to the logical count -- which is precisely the behaviour GH#66152
# exists to avoid.  These assert the native call actually answers on its own
# platform, so a broken sysctl name or a vanished sysfs layout fails the suite
# instead of silently reverting the default to logical cores.
@pytest.mark.skipif(sys.platform != "darwin", reason="macOS sysctl probe")
def test_physical_cores_darwin_probe_succeeds():
    assert _cpu._physical_cores_darwin() is not None


@pytest.mark.skipif(sys.platform != "linux", reason="Linux sysfs probe")
def test_physical_cores_linux_probe_succeeds():
    if not os.path.exists("/sys/devices/system/cpu/cpu0/topology/core_id"):
        pytest.skip("sysfs CPU topology not exposed here")
    assert _cpu._physical_cores_linux() is not None


@pytest.mark.skipif(sys.platform != "win32", reason="Windows ctypes probe")
def test_physical_cores_windows_probe_succeeds():
    # The synthetic buffers in the _count_processor_core_records tests are
    # built by a helper that assumes the same Size-at-+4 layout the production
    # code assumes, so the two can be wrong together.  Only the real call
    # catches a bad RelationProcessorCore value or argtypes/byref mistake,
    # which _physical_cores_windows' blanket ``except Exception`` would
    # otherwise turn into a silent fallback to the logical count.
    assert _cpu._physical_cores_windows() is not None


@pytest.mark.parametrize(
    "platform_name, probe_name",
    [
        ("darwin", "_physical_cores_darwin"),
        ("linux", "_physical_cores_linux"),
        ("win32", "_physical_cores_windows"),
    ],
)
@pytest.mark.parametrize("failure", ["returns_none", "raises"])
def test_physical_core_count_falls_back(
    monkeypatch, platform_name, probe_name, failure
):
    # A probe that returns None or raises must fall back to os.cpu_count().
    def probe():
        if failure == "raises":
            raise OSError("probe failed")

    monkeypatch.setattr(_cpu.sys, "platform", platform_name)
    monkeypatch.setattr(_cpu, probe_name, probe)
    physical_core_count.cache_clear()
    try:
        assert physical_core_count() == (os.cpu_count() or 1)
    finally:
        physical_core_count.cache_clear()


@pytest.mark.parametrize(
    "probe_result, expected",
    [
        (1, 1),  # a valid subset is used as-is
        (0, "total"),  # nonsensical count -> fall back
        (10**6, "total"),  # more cores than exist -> fall back
    ],
)
def test_physical_core_count_validates_probe(monkeypatch, probe_result, expected):
    total = os.cpu_count() or 1
    monkeypatch.setattr(_cpu.sys, "platform", "darwin")
    monkeypatch.setattr(_cpu, "_physical_cores_darwin", lambda: probe_result)
    physical_core_count.cache_clear()
    try:
        result = physical_core_count()
    finally:
        physical_core_count.cache_clear()
    assert result == (total if expected == "total" else expected)


@pytest.mark.parametrize(
    "spec, expected",
    [
        ("0-5", [0, 1, 2, 3, 4, 5]),
        ("0-3,8", [0, 1, 2, 3, 8]),
        ("0,2,4", [0, 2, 4]),
        ("3", [3]),
        ("", []),
    ],
)
def test_parse_cpu_list(spec, expected):
    assert _parse_cpu_list(spec) == expected


@pytest.mark.parametrize(
    "topology, expected",
    [
        # 4 SMT threads -> 2 physical cores (siblings collapse by (pkg, core)).
        ([(0, 0), (0, 0), (0, 1), (0, 1)], 2),
        # two packages, one core each
        ([(0, 0), (1, 0)], 2),
        # unreadable entries are ignored
        ([(0, 0), None, (0, 0)], 1),
        ([None, None], None),
        ([], None),
    ],
)
def test_count_distinct_cores(topology, expected):
    assert _count_distinct_cores(topology) == expected


def _fake_cpu_topology_reader(cores_per_cpu):
    """Build a fake ``_read_sysfs_int`` returning core/package ids from a map."""

    def read_int(path):
        match = re.search(r"/cpu(\d+)/topology/(\w+)", path)
        if match is None:
            return None
        cpu, field = int(match.group(1)), match.group(2)
        if field == "physical_package_id":
            return 0
        if field == "core_id":
            return cores_per_cpu.get(cpu)
        return None

    return read_int


def test_physical_cores_linux_hybrid_collapses_smt(monkeypatch):
    # Hybrid part: logical CPUs 0-7 are 4 physical P-cores with 2 SMT threads
    # each, and CPUs 8-11 are 4 efficiency cores.  SMT siblings collapse;
    # efficiency cores count.

    strs = {"/sys/devices/system/cpu/present": "0-11"}
    cores = {cpu: cpu // 2 for cpu in range(8)}
    cores.update({cpu: cpu - 4 for cpu in range(8, 12)})
    # Pin the affinity mask: unpatched, the real one decides which CPUs get
    # counted and the assertion depends on the host.
    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: set(range(12)), raising=False
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: strs.get(path))
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 8


def test_physical_cores_linux_respects_affinity(monkeypatch):
    # Same hybrid part, but the process is masked to CPUs 0-1 -- two SMT
    # siblings of a single physical core.  Counting all present cores here
    # would report 8 and let the (logical) availability clamp pick 2 workers
    # for 1 core.

    strs = {"/sys/devices/system/cpu/present": "0-11"}
    cores = {cpu: cpu // 2 for cpu in range(8)}
    cores.update({cpu: cpu - 4 for cpu in range(8, 12)})
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1}, raising=False)
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: strs.get(path))
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 1


def test_physical_cores_linux_homogeneous_uses_present(monkeypatch):
    # No affinity call available (or it fails): fall back to every present
    # CPU, SMT collapsed (0-3 = 2 physical cores).

    strs = {"/sys/devices/system/cpu/present": "0-3"}
    cores = {cpu: cpu // 2 for cpu in range(4)}
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: strs.get(path))
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 2


def test_physical_cores_linux_topology_unreadable(monkeypatch):
    # No topology/ files: report None and let physical_core_count()'s
    # os.cpu_count() fallback be the single fallback path.  Returning the
    # logical count instead would pick the same worker count in every case
    # (available_cpu_count() clamps both identically) -- one path is just
    # simpler than two.

    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: {0, 1, 2, 3}, raising=False
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: "0-3")
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: None)
    assert _cpu._physical_cores_linux() is None


def test_physical_cores_linux_degenerate_topology(monkeypatch):
    # Some hypervisors report physical_package_id=0/core_id=0 for every CPU.
    # Collapsing that to 1 physical core would make a 12-CPU guest read
    # serially -- slower than the previous flat default of 4 -- so an
    # implausible logical:physical ratio is treated as no answer at all.
    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: set(range(12)), raising=False
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: "0-11")
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: 0)
    assert _cpu._physical_cores_linux() is None


def test_physical_cores_linux_accepts_plausible_smt(monkeypatch):
    # The guard must not reject real SMT: 12 logical over 6 physical is 2
    # threads per core, well inside the bound.
    cores = {cpu: cpu // 2 for cpu in range(12)}
    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: set(range(12)), raising=False
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: "0-11")
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 6


def _make_win_processor_buffer(records):
    """Build a synthetic ``GetLogicalProcessorInformationEx`` buffer.

    ``records`` is a list of ``(size, efficiency_class)`` pairs, one per
    physical-core relationship record.  ``Size`` sits at byte offset ``+4`` and
    ``EfficiencyClass`` at ``+9`` within each record.
    """
    total = sum(size for size, _ in records)
    buf = (ctypes.c_byte * total)()
    addr = ctypes.addressof(buf)
    offset = 0
    for size, eff in records:
        ctypes.c_uint32.from_address(addr + offset + 4).value = size
        ctypes.c_uint8.from_address(addr + offset + 9).value = eff
        offset += size
    return buf, total


def test_count_processor_core_records_zero_size_not_counted():
    # A zero Size cannot advance the cursor, so the record is malformed and
    # must not be counted as a core.  Two 32-byte slots, second one zeroed.
    buf = (ctypes.c_byte * 64)()
    ctypes.c_uint32.from_address(ctypes.addressof(buf) + 4).value = 32
    assert _count_processor_core_records(buf, 64) == 1


def test_count_processor_core_records_truncated_header():
    # Size is a DWORD at +4, so a record header needs 8 readable bytes.  A
    # length stopping mid-header must yield None rather than read past it.
    buf, _ = _make_win_processor_buffer([(32, 0)])
    assert _count_processor_core_records(buf, 6) is None


def test_count_processor_core_records_hybrid():
    # 4 performance cores (class 1) + 4 efficiency cores (class 0): every
    # physical core counts, efficiency class does not matter.
    buf, length = _make_win_processor_buffer([(32, 1)] * 4 + [(32, 0)] * 4)
    assert _count_processor_core_records(buf, length) == 8


def test_count_processor_core_records_homogeneous():
    buf, length = _make_win_processor_buffer([(32, 0)] * 8)
    assert _count_processor_core_records(buf, length) == 8


def test_count_processor_core_records_variable_record_size():
    # Records advance by their own Size field, not a fixed stride.
    buf, length = _make_win_processor_buffer([(40, 1), (24, 1), (32, 0)])
    assert _count_processor_core_records(buf, length) == 3


def test_count_processor_core_records_empty():
    buf = (ctypes.c_byte * 0)()
    assert _count_processor_core_records(buf, 0) is None
