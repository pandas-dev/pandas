from __future__ import annotations

import pickle

import pytest

from dask._expr import HLGExpr, _ExprSequence
from dask._task_spec import (
    DataNode,
    DependenciesMapping,
    Task,
    TaskRef,
    fuse_linear_task_spec,
    resolve_aliases,
)
from dask.blockwise import blockwise, optimize_blockwise
from dask.core import flatten, reverse_dict
from dask.highlevelgraph import HighLevelGraph, MaterializedLayer
from dask.tokenize import tokenize
from dask.utils import ensure_dict


def test_tokenize():
    # Ensure tokens are different for different high-level graphs The current
    # implementation actually ensures that no HLGExpr are tokenizing equally.
    # Technically, we do not need such a strong guarantee. but tokenizing a full
    # HLG reliably is tricky and we do not require the reproducibility for
    # HLGExpr since they do not undergo the same kind of optimization as the
    # rest of the graph.
    from dask.highlevelgraph import HighLevelGraph

    dsk = HighLevelGraph.from_collections("x", {"foo": None})
    dsk2 = HighLevelGraph.from_collections("x", {"bar": None})
    dsk3 = HighLevelGraph.from_collections("y", {"foo": None})
    assert tokenize(HLGExpr(dsk)) != tokenize(HLGExpr(dsk2))
    assert tokenize(HLGExpr(dsk)) != tokenize(HLGExpr(dsk3))
    assert tokenize(HLGExpr(dsk2)) != tokenize(HLGExpr(dsk3))

    # Roundtrip preserves the tokens
    for expr in [HLGExpr(dsk), HLGExpr(dsk2), HLGExpr(dsk3)]:
        assert tokenize(pickle.loads(pickle.dumps(expr))) == tokenize(expr)


def optimizer(dsk, keys):
    dsk = ensure_dict(dsk)
    keys = list(flatten(keys))
    dsk = fuse_linear_task_spec(dsk, keys)
    return resolve_aliases(
        dsk,
        keys,
        reverse_dict(DependenciesMapping(dsk)),
    )


def optimizer2(dsk, keys):
    return optimizer(dsk, keys)


def func(*args, **kwargs):
    pass


def test_hlg_expr_sequence_finalize():
    hlgx = HighLevelGraph(
        xlayer := {"xlayer": MaterializedLayer({"x": DataNode("x", 1)})},
        dependencies=(xdeps := {"xlayer": set()}),
    )
    ylayer = {"ylayer": MaterializedLayer({"y": Task("y", func, TaskRef("x"))})}
    ylayer.update(xlayer)
    ydeps = {"ylayer": {"xlayer"}}
    ydeps.update(xdeps)
    hlgy = HighLevelGraph(ylayer, dependencies=ydeps)
    zlayer = {"zlayer": MaterializedLayer({"z": Task("z", func, TaskRef("x"))})}
    zlayer.update(xlayer)
    zdeps = {"zlayer": {"xlayer"}}

    zdeps.update(xdeps)
    hlgz = HighLevelGraph(zlayer, dependencies=zdeps)
    hlgexprx = HLGExpr(
        hlgx,
        low_level_optimizer=optimizer,
        output_keys=["x"],
    )
    hlgexpry = HLGExpr(
        hlgy,
        low_level_optimizer=optimizer,
        output_keys=["y"],
    )
    hlgexprz = HLGExpr(
        hlgz,
        low_level_optimizer=optimizer,
        output_keys=["z"],
    )
    dskx = hlgexprx.finalize_compute().optimize().__dask_graph__()
    assert isinstance(dskx, dict)
    assert len(dskx) == 1
    assert "x" in dskx
    assert dskx["x"] is hlgy.layers["xlayer"]["x"]

    dsky = hlgexpry.finalize_compute().optimize().__dask_graph__()
    assert isinstance(dsky, dict)
    # Linear low level fusion
    assert len(dsky) == 1
    assert "y" in dsky
    assert dsky["y"] != hlgy.layers["ylayer"]["y"]

    expryz_opt = _ExprSequence(hlgexprz, hlgexpry).finalize_compute().optimize()
    keys_yz = expryz_opt.__dask_keys__()
    assert len(keys_yz) == 2

    dskyz = expryz_opt.__dask_graph__()
    assert isinstance(dskyz, dict)
    expected = {}
    expected.update(next(iter(hlgx.layers.values())).mapping)
    expected.update(next(iter(hlgy.layers.values())).mapping)
    expected.update(next(iter(hlgz.layers.values())).mapping)
    # This is building the graph properly without fusing anything
    assert dskyz == expected

    hlgexprz_different_optimizer = HLGExpr(
        hlgz,
        low_level_optimizer=optimizer2,
        output_keys=["z"],
    )

    dskyz2 = (
        _ExprSequence(hlgexprz_different_optimizer, hlgexpry)
        .finalize_compute()
        .optimize()
        .__dask_graph__()
    )
    # both are fusing x
    assert "x" not in dskyz2
    assert len(dskyz2) == 2


def test_hlg_expr_sequence_nested_keys():
    xlayer = {"xlayer": MaterializedLayer({"x": DataNode("x", 1)})}
    xdeps = {"xlayer": set()}
    ylayer = {"ylayer": MaterializedLayer({"y": Task("y", func, TaskRef("x"))})}
    ylayer.update(xlayer)
    ydeps = {"ylayer": {"xlayer"}}
    ydeps.update(xdeps)
    hlgy = HighLevelGraph(ylayer, dependencies=ydeps)
    zlayer = {"zlayer": MaterializedLayer({"z": Task("z", func, TaskRef("x"))})}
    zlayer.update(xlayer)
    zdeps = {"zlayer": {"xlayer"}}

    zdeps.update(xdeps)
    hlgz = HighLevelGraph(zlayer, dependencies=zdeps)
    hlgexpry = HLGExpr(
        hlgy,
        low_level_optimizer=optimizer,
        output_keys=[["y"], ["x"]],
    )
    hlgexprz = HLGExpr(
        hlgz,
        low_level_optimizer=optimizer,
        output_keys=[["z"]],
    )
    expr = _ExprSequence(hlgexprz, hlgexpry)
    expected = [[["z"]], [["y"], ["x"]]]
    assert expr.__dask_keys__() == expected
    assert expr.optimize().__dask_keys__() == expected

    # Now with a different grouping / optimizer pass. These are handled
    # separately internalyl and we want to make sure the sequence is putting it
    # back together properly
    hlgexprz = HLGExpr(
        hlgz,
        low_level_optimizer=optimizer2,
        output_keys=[["z"]],
    )
    expr = _ExprSequence(hlgexprz, hlgexpry)
    assert expr.__dask_keys__() == expected
    assert expr.optimize().__dask_keys__() == expected


def optimizer_with_annotations(dsk, keys):
    annots = {}
    for layer in dsk.layers.values():
        if layer.annotations:
            annots.update(layer.annotations)
    out = optimizer(dsk, keys)
    out = {f"{key}-optimized": val for key, val in out.items()}
    return HighLevelGraph(
        layers={"fused": MaterializedLayer(out, annotations=annots)},
        dependencies={"fused": set()},
    )


def test_hlg_sequence_uses_annotations_of_optimized_dsk():
    """
    There was an issue where generating annotations would fall back to the HLG
    objects that were not optimized in combination.
    """
    xlayer = {
        "xlayer": MaterializedLayer({"x": DataNode("x", 1)}, annotations={"foo": 1})
    }
    xdeps = {"xlayer": set()}
    ylayer = {"ylayer": MaterializedLayer({"y": Task("y", func, TaskRef("x"))})}
    ylayer.update(xlayer)
    ydeps = {"ylayer": {"xlayer"}}
    ydeps.update(xdeps)
    hlgy = HighLevelGraph(ylayer, dependencies=ydeps)
    zlayer = {"zlayer": MaterializedLayer({"z": Task("z", func, TaskRef("x"))})}
    zlayer.update(xlayer)
    zdeps = {"zlayer": {"xlayer"}}

    zdeps.update(xdeps)
    hlgz = HighLevelGraph(zlayer, dependencies=zdeps)
    hlgexpry = HLGExpr(
        hlgy,
        output_keys=[["y"], ["x"]],
    )
    hlgexprz = HLGExpr(
        hlgz,
        output_keys=[["z"]],
    )
    expr = _ExprSequence(hlgexprz, hlgexpry)

    expected = {"foo": {"x": 1}}
    assert hlgexprz.__dask_annotations__() == expected
    assert hlgexprz.optimize().__dask_annotations__() == expected
    assert expr.__dask_annotations__() == expected
    assert expr.optimize().__dask_annotations__() == expected

    def only_optimize_with_z(dsk, keys):
        if "z" not in list(flatten(keys)):
            raise RuntimeError("Don't optimize me")
        return optimizer_with_annotations(dsk, keys)

    hlgexpry = HLGExpr(
        hlgy,
        low_level_optimizer=only_optimize_with_z,
        output_keys=[["y"], ["x"]],
    )
    hlgexprz = HLGExpr(
        hlgz,
        low_level_optimizer=only_optimize_with_z,
        output_keys=[["z"]],
    )
    expr = _ExprSequence(hlgexprz, hlgexpry)

    expected = {"foo": {"z-optimized": 1}}
    assert hlgexprz.__dask_annotations__() == expected
    assert hlgexprz.optimize().__dask_annotations__() == expected
    expected = {"foo": {f"{key}-optimized": 1 for key in list("xyz")}}
    with pytest.raises(RuntimeError, match="Don't optimize me"):
        expr.__dask_annotations__()
    assert expr.optimize().__dask_annotations__() == expected


def test_hlg_blockwise_fusion():
    b1 = blockwise(
        lambda x: x,
        "out1",
        "i",
        "x",
        "i",
        numblocks={"x": (2,)},
    )
    b2 = blockwise(
        lambda x: x,
        "out2",
        "i",
        "out1",
        "i",
        numblocks={"out1": (2,)},
    )

    layers = {"b1": b1, "b2": b2}
    deps = {
        "b1": set(),
        "b2": {"b1"},
    }
    hlg = HighLevelGraph(layers, deps)
    # Sanity check that the graph is correct
    assert len(dict(hlg)) == 4
    hlgexpr = HLGExpr(hlg, low_level_optimizer=optimize_blockwise)
    assert len(hlgexpr.__dask_graph__()) == 2
    assert len(hlgexpr.optimize().__dask_graph__()) == 2
    assert len(hlgexpr.finalize_compute().optimize().__dask_graph__()) == 2
