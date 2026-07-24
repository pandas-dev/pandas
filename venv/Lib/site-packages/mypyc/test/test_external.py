"""Test cases that run tests as subprocesses."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from distutils import ccompiler, sysconfig
from typing import Any

from mypyc.build import get_cflags, include_dir

from .config import PREFIX

EXCLUDED_LIB_RT_COMPILE_FILES = ["static_data.c"]


class TestExternal(unittest.TestCase):
    def test_lib_rt_c_files_compile_individually(self) -> None:
        """Compile each top-level lib-rt C file as its own translation unit."""
        lib_rt_dir = include_dir()
        source_names = sorted(
            name
            for name in os.listdir(lib_rt_dir)
            if name.endswith(".c")
            and os.path.isfile(os.path.join(lib_rt_dir, name))
            and name not in EXCLUDED_LIB_RT_COMPILE_FILES
        )
        compiler: Any = ccompiler.new_compiler()
        sysconfig.customize_compiler(compiler)

        include_dirs = [lib_rt_dir]
        for plat_specific in (False, True):
            path = sysconfig.get_python_inc(plat_specific=plat_specific)
            if path and path not in include_dirs:
                include_dirs.append(path)

        with tempfile.TemporaryDirectory() as tmpdir:
            for experimental_features in (False, True):
                cflags = get_cflags(
                    compiler_type=compiler.compiler_type,
                    opt_level="0",
                    experimental_features=experimental_features,
                )
                output_dir = os.path.join(
                    tmpdir, "experimental" if experimental_features else "default"
                )

                for source_name in source_names:
                    source_path = os.path.join(lib_rt_dir, source_name)
                    with self.subTest(source=source_name, experimental=experimental_features):
                        try:
                            compiler.compile(
                                [source_path],
                                output_dir=output_dir,
                                include_dirs=include_dirs,
                                extra_postargs=cflags,
                            )
                        except Exception as err:
                            raise AssertionError(
                                f"failed to compile {source_name} "
                                f"(experimental={experimental_features})"
                            ) from err

    # TODO: Get this to work on Windows.
    # (Or don't. It is probably not a good use of time.)
    @unittest.skipIf(sys.platform.startswith("win"), "rt tests don't work on windows")
    def test_c_unit_test(self) -> None:
        """Run C unit tests in a subprocess."""
        cppflags: list[str] = []
        env = os.environ.copy()
        if sys.platform == "darwin":
            cppflags += ["-O0", "-mmacosx-version-min=10.10", "-stdlib=libc++"]
        elif sys.platform == "linux":
            cppflags += ["-O0"]
        env["CPPFLAGS"] = " ".join(cppflags)
        # Build Python wrapper for C unit tests.

        with tempfile.TemporaryDirectory() as tmpdir:
            status = subprocess.check_call(
                [
                    sys.executable,
                    "setup.py",
                    "build_ext",
                    f"--build-lib={tmpdir}",
                    f"--build-temp={tmpdir}",
                    "--run-capi-tests",
                ],
                env=env,
                cwd=os.path.join(PREFIX, "mypyc", "lib-rt"),
            )
            # Run C unit tests.
            env = os.environ.copy()
            if "GTEST_COLOR" not in os.environ:
                env["GTEST_COLOR"] = "yes"  # Use fancy colors
            status = subprocess.call(
                [sys.executable, "-c", "import sys, test_capi; sys.exit(test_capi.run_tests())"],
                env=env,
                cwd=tmpdir,
            )
            if status != 0:
                raise AssertionError("make test: C unit test failure")
