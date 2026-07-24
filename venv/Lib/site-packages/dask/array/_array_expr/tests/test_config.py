import pytest

import dask.config
from dask.array import _array_expr_enabled, array_expr_enabled

pytestmark = pytest.mark.normal_and_array_expr


def test_array_expr_enabled():
    assert _array_expr_enabled() in (True, False)
    assert array_expr_enabled() is _array_expr_enabled()


def test_flip_config():
    """Config changes after the initial import are ignored with a warning"""
    prev = _array_expr_enabled()
    with dask.config.set({"array.query-planning": not prev}):
        with pytest.warns(RuntimeWarning, match="first imported"):
            actual = _array_expr_enabled()
    assert actual == prev


def test_str_config():
    """Prevent accidental use of 'false' string in YAML config,
    which would evaluate to True.
    """
    with dask.config.set({"array.query-planning": "false"}):
        with pytest.raises(TypeError):
            _array_expr_enabled()
