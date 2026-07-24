"""A basic check to make sure that we are using a mypyc-compiled version when expected."""

from __future__ import annotations

import filecmp
import os
from unittest import TestCase

import mypy
import mypyc


class MypycTest(TestCase):
    def test_using_mypyc(self) -> None:
        if os.getenv("TEST_MYPYC", None) == "1":
            assert not mypy.__file__.endswith(".py"), "Expected to find a mypyc-compiled version"

    def test_shared_files_consistent(self) -> None:
        if os.getenv("TEST_MYPYC", None) != "1":
            mypyc_path = mypyc.__path__[0]
            for f in ["build_setup.py"]:
                assert filecmp.cmp(
                    os.path.join(mypyc_path, f),
                    os.path.join(mypyc_path, "lib-rt", f),
                    shallow=False,
                ), f"Shared files inconsistent, run cp mypyc/{f} mypyc/lib-rt"
