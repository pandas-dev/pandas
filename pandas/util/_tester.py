"""
Entrypoint for testing from the top-level namespace
"""
import os

PKG = os.path.dirname(os.path.dirname(__file__))


try:
    import pytest
except ImportError:
    def test():
        raise ImportError("Need pytest>=3.0 to run tests")
else:
    def test(extra_args=None):
        if extra_args:
            cmd = ['-q'] + extra_args + [PKG]
        else:
            cmd = ['-q', PKG]
        pytest.main(cmd)


__all__ = ['test']
