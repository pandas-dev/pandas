from __future__ import annotations

import warnings

import pytest

from dask import utils_test
from dask.highlevelgraph import HighLevelGraph
from dask.utils_test import _check_warning


def test_hlg_layer():
    a = {"x": 1}
    b = {"y": (utils_test.inc, "x")}
    layers = {"a-layer": a, "bee-layer": b}
    dependencies = {"a-layer": set(), "bee-layer": {"a-layer"}}
    hg = HighLevelGraph(layers, dependencies)

    assert utils_test.hlg_layer(hg, "a") is hg.layers["a-layer"]
    assert utils_test.hlg_layer(hg, "b") is hg.layers["bee-layer"]
    with pytest.raises(KeyError, match="No layer starts with"):
        utils_test.hlg_layer(hg, "foo")


def test_hlg_layer_topological():
    a = {"x": 1}
    b = {"y": (utils_test.inc, "x")}
    c = {"z": (utils_test.inc, "x")}
    d = {"r": (sum, ["y", "z"])}
    layers = {"a": a, "b": b, "c": c, "d": d}
    dependencies = {"a": set(), "b": {"a"}, "c": {"a"}, "d": {"b", "c"}}
    hg = HighLevelGraph(layers, dependencies)

    assert utils_test.hlg_layer_topological(hg, -1) is hg.layers["d"]
    assert utils_test.hlg_layer_topological(hg, 0) is hg.layers["a"]
    assert utils_test.hlg_layer_topological(hg, 1) in (hg.layers["b"], hg.layers["c"])


def test__check_warning():
    class MyWarning(Warning):
        pass

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with _check_warning(True, MyWarning, "foo"):
            warnings.warn("foo", MyWarning)

    with pytest.warns(MyWarning, match="foo"):
        with _check_warning(False, MyWarning, "foo"):
            warnings.warn("foo", MyWarning)
