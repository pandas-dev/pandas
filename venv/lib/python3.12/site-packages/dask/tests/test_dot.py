from __future__ import annotations

import copy
import os
import re
import sys
from functools import partial
from operator import add, neg

import pytest

ipycytoscape = pytest.importorskip("ipycytoscape")
from dask.dot import _to_cytoscape_json, cytoscape_graph

if sys.flags.optimize != 2:
    pytest.importorskip("graphviz")
    from dask.dot import dot_graph, label, task_label, to_graphviz
else:
    pytestmark = pytest.mark.skipif(
        True, reason="graphviz exception with Python -OO flag"
    )

from dask import delayed
from dask.utils import ensure_not_exists

try:
    from IPython.display import SVG, Image
except ImportError:
    ipython_not_installed = True
    Image = None
    SVG = None
else:
    ipython_not_installed = False
ipython_not_installed_mark = pytest.mark.skipif(
    ipython_not_installed, reason="IPython not installed"
)


# Since graphviz doesn't store a graph, we need to parse the output
label_re = re.compile(r".*\[label=(.*?) shape=(.*?)\]")


def get_label(line):
    m = label_re.match(line)
    if m:
        return m.group(1)


def get_shape(line):
    m = label_re.match(line)
    if m:
        return m.group(2)


dsk = {
    "a": 1,
    "b": 2,
    "c": (neg, "a"),
    "d": (neg, "b"),
    "e": (add, "c", "d"),
    "f": (sum, ["a", "e"]),
}


def test_task_label():
    assert task_label((partial(add, 1), 1)) == "add"
    assert task_label((add, 1)) == "add"
    assert task_label((add, (add, 1, 2))) == "add(...)"


def test_label():
    assert label("x") == "x"
    assert label("elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487") == "elemwise-#"

    cache = {}
    result = label("elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487", cache=cache)
    assert result == "elemwise-#0"
    # cached
    result = label("elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487", cache=cache)
    assert result == "elemwise-#0"
    assert len(cache) == 1

    result = label("elemwise-e890b510984f344edea9a5e5fe05c0db", cache=cache)
    assert result == "elemwise-#1"
    assert len(cache) == 2

    result = label("elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487", cache=cache)
    assert result == "elemwise-#0"
    assert len(cache) == 2

    assert label("x", cache=cache) == "x"
    assert len(cache) == 2


def test_to_graphviz():
    g = to_graphviz(dsk)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 10  # 10 nodes total
    assert set(labels) == {"c", "d", "e", "f", '""'}
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == {"box", "circle"}


def test_to_graphviz_custom():
    g = to_graphviz(
        dsk,
        data_attributes={"a": {"shape": "square"}},
        function_attributes={"c": {"label": "neg_c", "shape": "ellipse"}},
    )
    labels = set(filter(None, map(get_label, g.body)))
    assert labels == {"neg_c", "d", "e", "f", '""'}
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == {"box", "circle", "square", "ellipse"}


def test_cytoscape_graph_custom(tmp_path):
    # Check that custom styling and layout propagates through
    g = cytoscape_graph(
        dsk,
        filename=os.fsdecode(tmp_path / "mydask.html"),
        rankdir="LR",
        node_sep=20,
        edge_sep=30,
        spacing_factor=2,
        edge_style={"line-color": "red"},
        node_style={"background-color": "green"},
    )
    sty = g.cytoscape_style
    layout = g.cytoscape_layout
    node_sty = next(s for s in sty if s["selector"] == "node")["style"]
    edge_sty = next(s for s in sty if s["selector"] == "edge")["style"]

    assert edge_sty["line-color"] == "red"
    assert node_sty["background-color"] == "green"
    assert layout["rankDir"] == "LR"
    assert layout["nodeSep"] == 20
    assert layout["edgeSep"] == 30
    assert layout["spacingFactor"] == 2


def test_cytoscape_graph_color():
    pytest.importorskip("matplotlib.pyplot")
    from dask.delayed import Delayed

    g = Delayed("f", dsk).visualize(engine="cytoscape")
    init_color = g.graph.nodes[0].data["color"]
    # Check that we changed the semantic color
    g = Delayed("f", dsk).visualize(engine="cytoscape", color="order")
    assert any(n.data["color"] != init_color for n in g.graph.nodes)


def test_to_graphviz_attributes():
    assert to_graphviz(dsk).graph_attr["rankdir"] == "BT"
    assert to_graphviz(dsk, rankdir="LR").graph_attr["rankdir"] == "LR"
    assert to_graphviz(dsk, node_attr={"color": "white"}).node_attr["color"] == "white"
    assert to_graphviz(dsk, edge_attr={"color": "white"}).edge_attr["color"] == "white"


def test_aliases():
    g = to_graphviz({"x": 1, "y": "x"})
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 2
    assert len(g.body) - len(labels) == 1  # Single edge


def test_to_graphviz_verbose():
    g = to_graphviz(dsk, verbose=True)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 10  # 10 nodes total
    assert set(labels) == {"a", "b", "c", "d", "e", "f"}
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == {"box", "circle"}


def test__to_cytoscape_json_verbose():
    data = _to_cytoscape_json(dsk, verbose=True)
    labels = list(map(lambda x: x["data"]["label"], data["nodes"]))
    assert len(labels) == 10
    assert set(labels) == {"a", "b", "c", "d", "e", "f"}
    shapes = list(map(lambda x: x["data"]["shape"], data["nodes"]))
    assert set(shapes) == {"ellipse", "rectangle"}


def test_to_graphviz_collapse_outputs():
    g = to_graphviz(dsk, collapse_outputs=True)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 6  # 6 nodes total
    assert set(labels) == {"c", "d", "e", "f", '""'}
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == {"box", "circle"}


def test__to_cytoscape_json_collapse_outputs():
    data = _to_cytoscape_json(dsk, collapse_outputs=True)
    labels = list(map(lambda x: x["data"]["label"], data["nodes"]))
    assert len(labels) == 6  # 6 nodes total
    assert set(labels) == {"c", "d", "e", "f", ""}
    shapes = list(map(lambda x: x["data"]["shape"], data["nodes"]))
    assert set(shapes) == {"ellipse", "rectangle"}


def test_to_graphviz_collapse_outputs_and_verbose():
    g = to_graphviz(dsk, collapse_outputs=True, verbose=True)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 6  # 6 nodes total
    assert set(labels) == {"a", "b", "c", "d", "e", "f"}
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == {"box", "circle"}


def test__to_cytoscape_json_collapse_outputs_and_verbose():
    data = _to_cytoscape_json(dsk, collapse_outputs=True, verbose=True)
    labels = list(map(lambda x: x["data"]["label"], data["nodes"]))
    assert len(labels) == 6  # 6 nodes total
    assert set(labels) == {"a", "b", "c", "d", "e", "f"}
    shapes = list(map(lambda x: x["data"]["shape"], data["nodes"]))
    assert set(shapes) == {"ellipse", "rectangle"}


def test_to_graphviz_with_unconnected_node():
    dsk2 = dsk.copy()
    dsk2["g"] = 3
    g = to_graphviz(dsk2, verbose=True)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 11  # 11 nodes total
    assert set(labels) == {"a", "b", "c", "d", "e", "f", "g"}

    g = to_graphviz(dsk, verbose=True, collapse_outputs=True)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 6  # 6 nodes total
    assert set(labels) == {"a", "b", "c", "d", "e", "f"}


@pytest.mark.parametrize(
    "format,typ",
    [
        pytest.param("png", Image, marks=ipython_not_installed_mark),
        pytest.param(
            "jpeg",
            Image,
            marks=pytest.mark.xfail(
                reason="jpeg not always supported in dot", strict=False
            ),
        ),
        ("dot", type(None)),
        ("pdf", type(None)),
        pytest.param("svg", SVG, marks=ipython_not_installed_mark),
    ],
)
@pytest.mark.xfail(
    sys.platform == "win32",
    reason="graphviz/pango on conda-forge currently broken for windows",
    strict=False,
)
def test_dot_graph(tmpdir, format, typ):
    # Use a name that the shell would interpret specially to ensure that we're
    # not vulnerable to shell injection when interacting with `dot`.
    filename = str(tmpdir.join("$(touch should_not_get_created.txt)"))

    target = ".".join([filename, format])
    ensure_not_exists(target)
    try:
        result = dot_graph(dsk, filename=filename, format=format)

        assert not os.path.exists("should_not_get_created.txt")
        assert os.path.isfile(target)
        assert isinstance(result, typ)
    finally:
        ensure_not_exists(target)


@pytest.mark.parametrize(
    "format,typ",
    [
        pytest.param("png", Image, marks=ipython_not_installed_mark),
        pytest.param(
            "jpeg",
            Image,
            marks=pytest.mark.xfail(
                reason="jpeg not always supported in dot", strict=False
            ),
        ),
        ("dot", type(None)),
        ("pdf", type(None)),
        pytest.param("svg", SVG, marks=ipython_not_installed_mark),
    ],
)
@pytest.mark.xfail(
    sys.platform == "win32",
    reason="graphviz/pango on conda-forge currently broken for windows",
    strict=False,
)
def test_dot_graph_no_filename(tmpdir, format, typ):
    before = tmpdir.listdir()
    result = dot_graph(dsk, filename=None, format=format)
    # We shouldn't write any files if filename is None.
    after = tmpdir.listdir()
    assert before == after
    assert isinstance(result, typ)


@ipython_not_installed_mark
@pytest.mark.xfail(
    sys.platform == "win32",
    reason="graphviz/pango on conda-forge currently broken for windows",
    strict=False,
)
def test_dot_graph_defaults():
    # Test with default args.
    default_name = "mydask"
    default_format = "png"
    target = ".".join([default_name, default_format])

    ensure_not_exists(target)
    try:
        result = dot_graph(dsk)
        assert os.path.isfile(target)
        assert isinstance(result, Image)
    finally:
        ensure_not_exists(target)


def test_cytoscape_graph(tmpdir):
    target = str(tmpdir.join("mydask.html"))
    ensure_not_exists(target)
    try:
        result = cytoscape_graph(dsk, target)
        assert os.path.isfile(target)
        assert isinstance(result, ipycytoscape.CytoscapeWidget)
    finally:
        ensure_not_exists(target)


@pytest.mark.parametrize(
    "filename,format,target,expected_result_type",
    [
        pytest.param(
            "mydaskpdf", "svg", "mydaskpdf.svg", SVG, marks=ipython_not_installed_mark
        ),
        ("mydask.pdf", None, "mydask.pdf", type(None)),
        pytest.param(
            "mydask.pdf", "svg", "mydask.pdf.svg", SVG, marks=ipython_not_installed_mark
        ),
        pytest.param(
            "mydaskpdf", None, "mydaskpdf.png", Image, marks=ipython_not_installed_mark
        ),
        pytest.param(
            "mydask.pdf.svg",
            None,
            "mydask.pdf.svg",
            SVG,
            marks=ipython_not_installed_mark,
        ),
    ],
)
@pytest.mark.xfail(
    sys.platform == "win32",
    reason="graphviz/pango on conda-forge currently broken for windows",
    strict=False,
)
def test_filenames_and_formats(tmpdir, filename, format, target, expected_result_type):
    result = dot_graph(dsk, filename=str(tmpdir.join(filename)), format=format)
    assert tmpdir.join(target).exists()
    assert isinstance(result, expected_result_type)


def test_delayed_kwargs_apply():
    def f(x, y=True):
        return x + y

    x = delayed(f)(1, y=2)
    label = task_label(x.dask[x.key])
    assert "f" in label
    assert "apply" not in label


@pytest.mark.parametrize("viz_func", [to_graphviz, _to_cytoscape_json])
def test_immutable_attributes(viz_func):
    def inc(x):
        return x + 1

    dsk = {"a": (inc, 1), "b": (inc, 2), "c": (add, "a", "b")}
    attrs_func = {"a": {}}
    attrs_data = {"b": {}}
    attrs_func_test = copy.deepcopy(attrs_func)
    attrs_data_test = copy.deepcopy(attrs_data)

    viz_func(
        dsk,
        function_attributes=attrs_func,
        data_attributes=attrs_data,
    )

    assert attrs_func_test == attrs_func
    assert attrs_data_test == attrs_data
