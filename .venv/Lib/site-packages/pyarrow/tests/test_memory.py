# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import contextlib
import os
import signal
import subprocess
import sys
import weakref

import pyarrow as pa
from pyarrow.tests import util

import pytest

pytestmark = pytest.mark.processes

possible_backends = ["system", "jemalloc", "mimalloc"]
# Backends which are expected to be present in all builds of PyArrow,
# except if the user manually recompiled Arrow C++.
mandatory_backends = ["system", "mimalloc"]


def backend_factory(backend_name):
    return getattr(pa, f"{backend_name}_memory_pool")


def supported_factories():
    yield pa.default_memory_pool
    for backend_name in pa.supported_memory_backends():
        yield backend_factory(backend_name)


@contextlib.contextmanager
def allocate_bytes(pool, nbytes):
    """
    Temporarily allocate *nbytes* from the given *pool*.
    """
    arr = pa.array([b"x" * nbytes], type=pa.binary(), memory_pool=pool)
    # Fetch the values buffer from the varbinary array and release the rest,
    # to get the desired allocation amount
    buf = arr.buffers()[2]
    arr = None
    assert len(buf) == nbytes
    try:
        yield
    finally:
        buf = None


def check_allocated_bytes(pool):
    """
    Check allocation stats on *pool*.
    """
    allocated_before = pool.bytes_allocated()
    max_mem_before = pool.max_memory()
    num_allocations_before = pool.num_allocations()
    with allocate_bytes(pool, 512):
        assert pool.bytes_allocated() == allocated_before + 512
        new_max_memory = pool.max_memory()
        assert pool.max_memory() >= max_mem_before
        num_allocations_after = pool.num_allocations()
        assert num_allocations_after > num_allocations_before
        assert num_allocations_after < num_allocations_before + 5
    assert pool.bytes_allocated() == allocated_before
    assert pool.max_memory() == new_max_memory
    assert pool.num_allocations() == num_allocations_after


def test_default_allocated_bytes():
    pool = pa.default_memory_pool()
    with allocate_bytes(pool, 1024):
        check_allocated_bytes(pool)
        assert pool.bytes_allocated() == pa.total_allocated_bytes()


def test_proxy_memory_pool():
    pool = pa.proxy_memory_pool(pa.default_memory_pool())
    check_allocated_bytes(pool)
    wr = weakref.ref(pool)
    assert wr() is not None
    del pool
    assert wr() is None


def test_logging_memory_pool(capfd):
    pool = pa.logging_memory_pool(pa.default_memory_pool())
    check_allocated_bytes(pool)
    out, err = capfd.readouterr()
    assert err == ""
    assert out.count("Allocate:") > 0
    assert out.count("Allocate:") == out.count("Free:")


def test_set_memory_pool():
    old_pool = pa.default_memory_pool()
    pool = pa.proxy_memory_pool(old_pool)
    pa.set_memory_pool(pool)
    try:
        allocated_before = pool.bytes_allocated()
        with allocate_bytes(None, 512):
            assert pool.bytes_allocated() == allocated_before + 512
        assert pool.bytes_allocated() == allocated_before
    finally:
        pa.set_memory_pool(old_pool)


def test_default_backend_name():
    pool = pa.default_memory_pool()
    assert pool.backend_name in possible_backends


def test_release_unused():
    pool = pa.default_memory_pool()
    pool.release_unused()


def check_env_var(name, expected, *, expect_warning=False):
    code = f"""if 1:
        import pyarrow as pa

        pool = pa.default_memory_pool()
        assert pool.backend_name in {expected!r}, pool.backend_name
        """
    env = dict(os.environ)
    env['ARROW_DEFAULT_MEMORY_POOL'] = name
    res = subprocess.run([sys.executable, "-c", code], env=env,
                         universal_newlines=True, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr, file=sys.stderr)
        res.check_returncode()  # fail
    errlines = res.stderr.splitlines()
    if expect_warning:
        assert len(errlines) in (1, 2)
        if len(errlines) == 1:
            # ARROW_USE_GLOG=OFF
            assert f"Unsupported backend '{name}'" in errlines[0]
        else:
            # ARROW_USE_GLOG=ON
            assert "InitGoogleLogging()" in errlines[0]
            assert f"Unsupported backend '{name}'" in errlines[1]
    else:
        assert len(errlines) == 0


def test_env_var():
    for backend_name in mandatory_backends:
        check_env_var(backend_name, [backend_name])
    check_env_var("nonexistent", possible_backends, expect_warning=True)


def test_memory_pool_factories():
    def check(factory, name, *, can_fail=False):
        if can_fail:
            try:
                pool = factory()
            except NotImplementedError:
                return
        else:
            pool = factory()
        assert pool.backend_name == name

    for backend_name in possible_backends:
        check(backend_factory(backend_name), backend_name,
              can_fail=backend_name not in mandatory_backends)


def test_supported_memory_backends():
    backends = pa.supported_memory_backends()
    assert set(backends) >= set(mandatory_backends)
    assert set(backends) <= set(possible_backends)


def run_debug_memory_pool(pool_factory, env_value):
    """
    Run a piece of code making an invalid memory write with the
    ARROW_DEBUG_MEMORY_POOL environment variable set to a specific value.
    """
    code = f"""if 1:
        import ctypes
        import pyarrow as pa
        # ARROW-16873: some Python installs enable faulthandler by default,
        # which could dump a spurious stack trace if the following crashes
        import faulthandler
        faulthandler.disable()

        pool = pa.{pool_factory}()
        buf = pa.allocate_buffer(64, memory_pool=pool)

        # Write memory out of bounds
        ptr = ctypes.cast(buf.address, ctypes.POINTER(ctypes.c_ubyte))
        ptr[64] = 0

        del buf
        """
    env = dict(os.environ)
    env['ARROW_DEBUG_MEMORY_POOL'] = env_value
    res = subprocess.run([sys.executable, "-c", code], env=env,
                         universal_newlines=True, stderr=subprocess.PIPE)
    print(res.stderr, file=sys.stderr)
    return res


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_abort(pool_factory):
    res = run_debug_memory_pool(pool_factory.__name__, "abort")
    if os.name == "posix":
        assert res.returncode == -signal.SIGABRT
    else:
        assert res.returncode != 0
    assert "Wrong size on deallocation" in res.stderr


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_trap(pool_factory):
    res = run_debug_memory_pool(pool_factory.__name__, "trap")
    if os.name == "posix":
        assert res.returncode == -signal.SIGTRAP
    else:
        assert res.returncode != 0
    assert "Wrong size on deallocation" in res.stderr


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_warn(pool_factory):
    res = run_debug_memory_pool(pool_factory.__name__, "warn")
    res.check_returncode()
    assert "Wrong size on deallocation" in res.stderr


def check_debug_memory_pool_disabled(pool_factory, env_value, msg):
    if sys.maxsize < 2**32:
        # GH-45011: mimalloc may print warnings in this test on 32-bit Linux, ignore.
        pytest.skip("Test may fail on 32-bit platforms")
    res = run_debug_memory_pool(pool_factory.__name__, env_value)
    # The subprocess either returned successfully or was killed by a signal
    # (due to writing out of bounds), depending on the underlying allocator.
    if os.name == "posix":
        assert res.returncode <= 0
    else:
        res.check_returncode()
    if msg == "":
        assert res.stderr == ""
    else:
        assert msg in res.stderr


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_none(pool_factory):
    check_debug_memory_pool_disabled(pool_factory, "none", "")


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_empty(pool_factory):
    check_debug_memory_pool_disabled(pool_factory, "", "")


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_debug_memory_pool_unknown(pool_factory):
    env_value = "some_arbitrary_value"
    msg = (
        f"Invalid value for ARROW_DEBUG_MEMORY_POOL: '{env_value}'. "
        "Valid values are 'abort', 'trap', 'warn', 'none'."
    )
    check_debug_memory_pool_disabled(pool_factory, env_value, msg)


@pytest.mark.parametrize('pool_factory', supported_factories())
def test_print_stats(pool_factory):
    code = f"""if 1:
        import pyarrow as pa

        pool = pa.{pool_factory.__name__}()
        buf = pa.allocate_buffer(64, memory_pool=pool)
        pool.print_stats()
        """
    res = subprocess.run([sys.executable, "-c", code], check=True,
                         universal_newlines=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    if sys.platform == "linux" and not util.running_on_musllinux():
        # On Linux with glibc at least, all memory pools should emit statistics
        assert res.stderr.strip() != ""
