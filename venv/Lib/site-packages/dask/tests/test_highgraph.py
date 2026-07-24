from __future__ import annotations

import os
import threading
import xml.etree.ElementTree
from collections.abc import Set
from concurrent.futures import ThreadPoolExecutor

import pytest

import dask
from dask.base import collections_to_expr, tokenize
from dask.blockwise import Blockwise
from dask.delayed import Delayed
from dask.highlevelgraph import HighLevelGraph, Layer, MaterializedLayer, to_graphviz
from dask.utils_test import dec, inc


def test_visualize(tmpdir):
    pytest.importorskip("numpy")
    pytest.importorskip("graphviz")
    da = pytest.importorskip("dask.array")
    fn = str(tmpdir)
    a = da.ones(10, chunks=(5,))
    b = a + 1
    c = a + 2
    d = b + c
    d.dask.visualize(fn)
    assert os.path.exists(fn)


def test_basic():
    a = {"x": 1}
    b = {"y": (inc, "x")}
    layers = {"a": a, "b": b}
    dependencies = {"a": set(), "b": {"a"}}
    hg = HighLevelGraph(layers, dependencies)

    assert dict(hg) == {"x": 1, "y": (inc, "x")}
    assert all(isinstance(layer, Layer) for layer in hg.layers.values())


def test_keys_values_items_to_dict_methods():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    a = da.ones(10, chunks=(5,))
    b = a + 1
    c = a + 2
    d = b + c
    hg = d.dask

    keys, values, items = hg.keys(), hg.values(), hg.items()
    assert isinstance(keys, Set)
    assert list(keys) == list(hg)
    assert list(values) == [hg[i] for i in hg]
    assert list(items) == list(zip(keys, values))
    assert hg.to_dict() == dict(hg)


def test_getitem():
    hg = HighLevelGraph(
        {"a": {"a": 1, ("a", 0): 2, "b": 3}, "b": {"c": 4}}, {"a": set(), "b": set()}
    )
    # Key is a string and it exists in a layer with the same name
    assert hg["a"] == 1
    # Key is a tuple and the name exists in a layer with the same name
    assert hg["a", 0] == 2
    # Key is in the wrong layer, while the right layer does not contain it
    assert hg["b"] == 3
    # Key is in the wrong layer, while the right layer does not exist
    assert hg["c"] == 4

    for k in ("d", "", 1, ()):
        with pytest.raises(KeyError):
            hg[k]

    class Unhashable:
        __hash__ = None

    for k in (Unhashable(), (Unhashable(),)):
        with pytest.raises(TypeError):
            hg[k]


def test_copy():
    h1 = HighLevelGraph(
        {"a": {"a": "b"}, "b": {"b": 1}},
        {"a": {"b"}, "b": set()},
    )
    h1.get_all_dependencies()
    assert h1.key_dependencies
    h2 = h1.copy()
    for k in ("layers", "dependencies", "key_dependencies"):
        v1 = getattr(h1, k)
        v2 = getattr(h2, k)
        assert v1 is not v2
        assert v1 == v2


def test_cull():
    a = {"x": 1, "y": (inc, "x")}
    hg = HighLevelGraph({"a": a}, {"a": set()})

    culled_by_x = hg.cull({"x"})
    assert dict(culled_by_x) == {"x": 1}

    # parameter is the raw output of __dask_keys__()
    culled_by_y = hg.cull([[["y"]]])
    assert dict(culled_by_y) == a


def test_cull_layers():
    hg = HighLevelGraph(
        {
            "a": {"a1": "d1", "a2": "e1"},
            "b": {"b": "d", "dontcull_b": 1},
            "c": {"dontcull_c": 1},
            "d": {"d": 1, "dontcull_d": 1},
            "e": {"e": 1, "dontcull_e": 1},
        },
        {"a": {"d", "e"}, "b": {"d"}, "c": set(), "d": set(), "e": set()},
    )

    # Deep-copy layers before calling method to test they aren't modified in place
    expect = HighLevelGraph(
        {k: dict(v) for k, v in hg.layers.items() if k != "c"},
        {k: set(v) for k, v in hg.dependencies.items() if k != "c"},
    )

    culled = hg.cull_layers(["a", "b"])

    assert culled.layers == expect.layers
    assert culled.dependencies == expect.dependencies
    for k in culled.layers:
        assert culled.layers[k] is hg.layers[k]
        assert culled.dependencies[k] is hg.dependencies[k]


def test_repr_html_hlg_layers():
    pytest.importorskip("jinja2")
    hg = HighLevelGraph(
        {"a": {"a": 1, ("a", 0): 2, "b": 3}, "b": {"c": 4}},
        {"a": set(), "b": set()},
    )
    assert xml.etree.ElementTree.fromstring(hg._repr_html_()) is not None
    for layer in hg.layers.values():
        assert xml.etree.ElementTree.fromstring(layer._repr_html_()) is not None


def annot_map_fn(key):
    return key[1:]


@pytest.mark.parametrize(
    "annotation",
    [
        {"worker": "alice"},
        {"block_id": annot_map_fn},
    ],
)
def test_single_annotation(annotation):
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    with dask.annotate(**annotation):
        A = da.ones((10, 10), chunks=(5, 5))

    alayer = A.__dask_graph__().layers[A.name]
    assert alayer.annotations == annotation
    assert not dask.get_annotations()


def test_multiple_annotations():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    with dask.annotate(block_id=annot_map_fn):
        with dask.annotate(resources={"GPU": 1}):
            A = da.ones((10, 10), chunks=(5, 5))

        B = A + 1

    C = B + 1

    assert not dask.get_annotations()

    alayer = A.__dask_graph__().layers[A.name]
    blayer = B.__dask_graph__().layers[B.name]
    clayer = C.__dask_graph__().layers[C.name]
    assert alayer.annotations == {"resources": {"GPU": 1}, "block_id": annot_map_fn}
    assert blayer.annotations == {"block_id": annot_map_fn}
    assert clayer.annotations is None


def test_annotation_cleared_on_error():
    with dask.annotate(x=1):
        with pytest.raises(ZeroDivisionError):
            with dask.annotate(x=2):
                assert dask.get_annotations() == {"x": 2}
                1 / 0
        assert dask.get_annotations() == {"x": 1}
    assert not dask.get_annotations()


def test_materializedlayer_cull_preserves_annotations():
    layer = MaterializedLayer(
        {"a": 42, "b": 3.14},
        annotations={"foo": "bar"},
    )

    culled_layer, _ = layer.cull({"a"}, [])
    assert len(culled_layer) == 1
    assert culled_layer.annotations == {"foo": "bar"}


def test_annotations_leak():
    """Annotations shouldn't leak between threads.
    See https://github.com/dask/dask/issues/10340."""
    b1 = threading.Barrier(2)
    b2 = threading.Barrier(2)

    def f(n):
        with dask.annotate(foo=n):
            b1.wait()
            out = dask.get_annotations()
            b2.wait()
            return out

    with ThreadPoolExecutor(2) as ex:
        f1 = ex.submit(f, 1)
        f2 = ex.submit(f, 2)
        result = [f1.result(), f2.result()]
    assert result == [{"foo": 1}, {"foo": 2}]


@pytest.mark.parametrize("flat", [True, False])
def test_blockwise_cull(flat):
    np = pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    if flat:
        # Simple "flat" mapping between input and
        # output indices
        x = da.from_array(np.arange(40).reshape((4, 10)), (2, 4)) + 100
    else:
        # Complex mapping between input and output
        # indices (outer product and transpose)
        x = da.from_array(np.arange(10).reshape((10,)), (4,))
        y = da.from_array(np.arange(10).reshape((10,)), (4,))
        x = da.outer(x, y).transpose()

    # Check that blockwise culling results in correct
    # output keys and that full graph is not materialized
    dsk = x.__dask_graph__()
    select = (1, 1)  # Select a single chunk
    keys = {(x._name, *select)}
    dsk_cull = dsk.cull(keys)
    for name, layer in dsk_cull.layers.items():
        old_name = name.rsplit("-", 1)[0]
        if not isinstance(layer, dask.blockwise.Blockwise):
            # The original layer shouldn't be Blockwise if the new one isn't
            assert not isinstance(dsk.layers[old_name], dask.blockwise.Blockwise)
            continue
        assert isinstance(dsk.layers[old_name], dask.blockwise.Blockwise)
        assert not layer.is_materialized()
        out_keys = layer.get_output_keys()
        assert out_keys == {(layer.output, *select)}
        assert not layer.is_materialized()


def test_len_does_not_materialize():
    from dask._task_spec import Task

    a = {"x": 1}
    b = Blockwise(
        output="b",
        output_indices=tuple("ij"),
        task=Task("b", lambda: "1"),
        indices=(),
        numblocks={},
        new_axes={"i": (1, 1, 1), "j": (1, 1)},
    )
    assert len(b) == len(b.get_output_keys())

    layers = {"a": a, "b": b}
    dependencies = {"a": set(), "b": {"a"}}
    hg = HighLevelGraph(layers, dependencies)

    assert hg.layers["a"].is_materialized()
    assert not hg.layers["b"].is_materialized()

    assert len(hg) == len(a) + len(b) == 7

    assert not hg.layers["b"].is_materialized()


def test_node_tooltips_exist():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    pytest.importorskip("graphviz")

    a = da.ones((1000, 1000), chunks=(100, 100))
    b = a + a.T
    c = b.sum(axis=1)

    hg = c.dask
    g = to_graphviz(hg)

    for layer in g.body:
        if "label" in layer:
            assert "tooltip" in layer
            start = layer.find('tooltip="') + len('tooltip="')
            end = layer.find('"', start)
            tooltip = layer[start:end]
            assert len(tooltip) > 0


def test_tokenize_hlg():
    import dask.bag as db

    a = db.from_sequence(list(range(10)), npartitions=2).max()
    b = db.from_sequence(list(range(10)), npartitions=2).max()
    c = db.from_sequence(list(range(10)), npartitions=3).max()
    assert tokenize(a.dask) == tokenize(b.dask)
    assert tokenize(a.dask) != tokenize(c.dask)


def test_culling_changes_layer_names():
    """Test that culling changes the layer names in the graph."""
    dsk = MaterializedLayer({"x": 1, "y": (inc, "x"), "z": (dec, "x")})
    hg = HighLevelGraph({"a": dsk}, {"a": set()})

    # This construction of Delayed objects is similar to array.to_delayed
    y = Delayed("y", hg, layer="a")
    z = Delayed("z", hg, layer="a")

    from dask.delayed import optimize

    # This is roughly what dask.optimize is doing
    # The exception is that collection_to_dsk takes the Delayed object and
    # extracts the layer as a low level graph in which case we're losing the
    # layer name, i.e. all internal manipulations of the graphs with and without
    # optimization could be affected but triggering this with user APIs is
    # almost impossible
    yopt_hlg = optimize(y._dask, ["y"])
    assert isinstance(yopt_hlg, HighLevelGraph)
    layers = yopt_hlg.layers
    # This test doesn't care about there being a single layer but after
    # optimization we no longer want to hardcode the delayed layer key
    assert len(layers) == 1
    layer_key = next(iter(layers))
    yopt = Delayed("y", yopt_hlg, layer=layer_key)

    zopt_hlg = optimize(z._dask, ["z"])
    assert isinstance(zopt_hlg, HighLevelGraph)
    layers = zopt_hlg.layers
    # This test doesn't care about there being a single layer but after
    # optimization we no longer want to hardcode the delayed layer key
    assert len(layers) == 1
    layer_key = next(iter(layers))
    zopt = Delayed("z", zopt_hlg, layer=layer_key)

    # This internally merged the HLG graphs i.e. if layer names are not
    # deduplicated we will loose keys.
    expr_opt = collections_to_expr([yopt, zopt]).optimize()
    dsk_out = expr_opt.__dask_graph__()
    assert set(dsk_out) == set("xyz")
