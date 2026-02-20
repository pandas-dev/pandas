"""
Generic test utilities.

Based on scipy._libs._testutils
"""

from __future__ import division, print_function, absolute_import

import os
import sys


__all__ = ["PytestTester"]


class PytestTester(object):
    """
    Pytest test runner entry point.
    """

    def __init__(self, module_name):
        self.module_name = module_name

    def __call__(
        self,
        label="fast",
        verbose=1,
        extra_argv=None,
        doctests=False,
        coverage=False,
        tests=None,
        parallel=None,
    ):
        import pytest

        module = sys.modules[self.module_name]
        module_path = os.path.abspath(module.__path__[0])

        pytest_args = ["-l"]

        if doctests:
            raise ValueError("Doctests not supported")

        if extra_argv:
            pytest_args += list(extra_argv)

        if verbose and int(verbose) > 1:
            pytest_args += ["-" + "v" * (int(verbose) - 1)]

        if coverage:
            pytest_args += ["--cov=" + module_path]

        if label == "fast":
            pytest_args += ["-m", "not slow"]
        elif label != "full":
            pytest_args += ["-m", label]

        if tests is None:
            tests = [self.module_name]

        if parallel is not None and parallel > 1:
            if _pytest_has_xdist():
                pytest_args += ["-n", str(parallel)]
            else:
                import warnings

                warnings.warn(
                    "Could not run tests in parallel because "
                    "pytest-xdist plugin is not available."
                )

        pytest_args += ["--pyargs"] + list(tests)

        try:
            code = pytest.main(pytest_args)
        except SystemExit as exc:
            code = exc.code

        return code == 0


def _pytest_has_xdist():
    """
    Check if the pytest-xdist plugin is installed, providing parallel tests
    """
    # Check xdist exists without importing, otherwise pytests emits warnings
    from importlib.util import find_spec

    return find_spec("xdist") is not None
