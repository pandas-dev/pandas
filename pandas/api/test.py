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
    def test():
        pytest.main(['-q', PKG])


__all__ = ['test']
