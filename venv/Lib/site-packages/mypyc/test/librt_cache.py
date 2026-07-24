"""Build and cache librt for use in tests.

This module provides a way to build librt extension modules once and cache
them across test runs, and across different test cases in a single run. The
cache is invalidated when source files or details of the build environment change.

Note: Tests must run in a subprocess to use the cached librt, since importing
this module also triggers the import of the regular installed librt.

Usage:
    from mypyc.test.librt_cache import get_librt_path, run_with_librt

    # Get path to built librt (builds if needed)
    path = get_librt_path()

    # Run a test file in subprocess with built librt
    result = run_with_librt("test_librt.py")
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import sysconfig
from typing import Any

import filelock

from mypyc.build import LIBRT_MODULES, get_cflags, include_dir
from mypyc.common import RUNTIME_C_FILES
from mypyc.test.config import PREFIX


def _librt_build_hash(experimental: bool, opt_level: str) -> str:
    """Compute hash for librt build, including sources and build environment."""
    # Import lazily to ensure mypyc.build has ensured that distutils is correctly set up
    from distutils import ccompiler

    h = hashlib.sha256()
    # Include experimental flag
    h.update(b"exp" if experimental else b"noexp")
    h.update(f"opt={opt_level}".encode())
    # Include full Python version string (includes git hash for dev builds)
    h.update(sys.version.encode())
    # Include debug build status (gettotalrefcount only exists in debug builds)
    is_debug = hasattr(sys, "gettotalrefcount")
    h.update(b"debug" if is_debug else b"release")
    # Include free-threading status (Python 3.13+)
    is_free_threaded = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
    h.update(b"freethreaded" if is_free_threaded else b"gil")
    # Include compiler type (e.g., "unix" or "msvc")
    compiler: Any = ccompiler.new_compiler()
    h.update(compiler.compiler_type.encode())
    # Include environment variables that affect C compilation
    for var in ("CC", "CXX", "CFLAGS", "CPPFLAGS", "LDFLAGS"):
        val = os.environ.get(var, "")
        h.update(f"{var}={val}".encode())
    # Hash runtime files
    for name in RUNTIME_C_FILES:
        path = os.path.join(include_dir(), name)
        h.update(name.encode() + b"|")
        with open(path, "rb") as f:
            h.update(f.read())
    # Hash librt module files
    for mod, files, extra, includes in LIBRT_MODULES:
        for fname in files + extra:
            path = os.path.join(include_dir(), fname)
            h.update(fname.encode() + b"|")
            with open(path, "rb") as f:
                h.update(f.read())
    return h.hexdigest()[:16]


def _generate_setup_py(build_dir: str, experimental: bool, opt_level: str) -> str:
    """Generate setup.py content for building librt directly.

    We inline LIBRT_MODULES/RUNTIME_C_FILES/include_dir/cflags values to avoid
    importing mypyc.build, which recursively imports lots of things.
    """
    lib_rt_dir = include_dir()

    # Get compiler flags using the shared helper
    cflags = get_cflags(opt_level=opt_level, experimental_features=experimental)

    # Serialize values to inline in generated setup.py
    librt_modules_repr = repr(
        [(m.module, m.c_files, m.other_files, m.include_dirs) for m in LIBRT_MODULES]
    )
    runtime_files_repr = repr(RUNTIME_C_FILES)
    cflags_repr = repr(cflags)

    return f"""\
import os
from setuptools import setup, Extension
import build_setup  # noqa: F401  # Monkey-patches compiler for per-file SIMD flags

build_dir = {build_dir!r}
lib_rt_dir = {lib_rt_dir!r}

RUNTIME_C_FILES = {runtime_files_repr}
LIBRT_MODULES = {librt_modules_repr}
CFLAGS = {cflags_repr}

def write_file(path, contents):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(contents)

# Copy runtime C files
for name in RUNTIME_C_FILES:
    src = os.path.join(lib_rt_dir, name)
    dst = os.path.join(build_dir, name)
    with open(src, "rb") as f:
        write_file(dst, f.read())

# Build extensions for each librt module
extensions = []
for mod, file_names, extra_files, includes in LIBRT_MODULES:
    # Copy source files
    for fname in file_names + extra_files:
        src = os.path.join(lib_rt_dir, fname)
        dst = os.path.join(build_dir, fname)
        with open(src, "rb") as f:
            write_file(dst, f.read())

    extensions.append(Extension(
        mod,
        sources=[os.path.join(build_dir, f) for f in file_names + RUNTIME_C_FILES],
        include_dirs=[lib_rt_dir] + [os.path.join(lib_rt_dir, d) for d in includes],
        extra_compile_args=CFLAGS,
    ))

setup(name='librt_cached', ext_modules=extensions)
"""


def get_librt_path(experimental: bool = True, opt_level: str = "0") -> str:
    """Get path to librt built from the repository, building and caching if necessary.

    Uses build/librt-cache/ under the repo root (gitignored). The cache is
    keyed by a hash of sources and build environment, so it auto-invalidates
    when relevant factors change.

    Safe to call from multiple parallel pytest workers - uses file locking.

    Args:
        experimental: Whether to enable experimental features.
        opt_level: Optimization level ("0".."3") used when building librt.

    Returns:
        Path to directory containing built librt modules.
    """
    # Use build/librt-cache/ under the repo root (gitignored)
    cache_root = os.path.join(PREFIX, "build", "librt-cache")
    build_hash = _librt_build_hash(experimental, opt_level)
    build_dir = os.path.join(cache_root, f"librt-{build_hash}")
    lock_file = os.path.join(cache_root, f"librt-{build_hash}.lock")
    marker = os.path.join(build_dir, ".complete")

    os.makedirs(cache_root, exist_ok=True)

    with filelock.FileLock(lock_file, timeout=300):  # 5 min timeout
        if os.path.exists(marker):
            return build_dir

        # Clean up any partial build
        if os.path.exists(build_dir):
            shutil.rmtree(build_dir)

        os.makedirs(build_dir)

        # Create librt package directory for --inplace to copy .so files into
        librt_pkg = os.path.join(build_dir, "librt")
        os.makedirs(librt_pkg)
        with open(os.path.join(librt_pkg, "__init__.py"), "w") as f:
            pass

        # Copy build_setup.py for per-file SIMD compiler flags
        build_setup_src = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build_setup.py"
        )
        build_setup_dst = os.path.join(build_dir, "build_setup.py")
        shutil.copy(build_setup_src, build_setup_dst)

        # Write setup.py
        setup_py = os.path.join(build_dir, "setup.py")
        with open(setup_py, "w") as f:
            f.write(_generate_setup_py(build_dir, experimental, opt_level))

        # Build (parallel builds don't work well because multiple extensions
        # share the same runtime C files, causing race conditions)
        result = subprocess.run(
            [sys.executable, setup_py, "build_ext", "--inplace"],
            cwd=build_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"librt build failed:\n{result.stdout}\n{result.stderr}")

        # Mark complete
        with open(marker, "w") as f:
            f.write("ok")

    return build_dir


def run_with_librt(
    file_path: str, experimental: bool = True, check: bool = True, opt_level: str = "0"
) -> subprocess.CompletedProcess[str]:
    """Run a Python file in a subprocess with built librt available.

    This runs the file in a fresh Python process where the built librt
    is at the front of sys.path, avoiding conflicts with any system librt.

    Args:
        file_path: Path to Python file to execute.
        experimental: Whether to use experimental features.
        check: If True, raise CalledProcessError on non-zero exit.
        opt_level: Optimization level ("0".."3") used when building librt.

    Returns:
        CompletedProcess with stdout, stderr, and returncode.
    """
    librt_path = get_librt_path(experimental, opt_level=opt_level)
    # Prepend librt path to PYTHONPATH
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = librt_path + (os.pathsep + existing if existing else "")

    return subprocess.run(
        [sys.executable, file_path], capture_output=True, text=True, check=check, env=env
    )
