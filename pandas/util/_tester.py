"""
Entrypoint for testing from the top-level namespace
"""
import os
import sys

PKG = os.path.dirname(os.path.dirname(__file__))


def test(extra_args=None):
    try:
        # error: No library stub file for module 'pytest'
        import pytest  # type: ignore
    except ImportError:
        raise ImportError("Need pytest>=4.0.2 to run tests")
    try:
        # error: Cannot find module named 'hypothesis'
        import hypothesis  # type:ignore # noqa
    except ImportError:
        raise ImportError("Need hypothesis>=3.58 to run tests")
    cmd = ["--skip-slow", "--skip-network", "--skip-db"]
    if extra_args:
        if not isinstance(extra_args, list):
            extra_args = [extra_args]
        cmd = extra_args
    cmd += [PKG]
    print("running: pytest {}".format(" ".join(cmd)))
    sys.exit(pytest.main(cmd))


__all__ = ["test"]
