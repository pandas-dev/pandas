from __future__ import annotations

import random
import time
from operator import add

import pytest

import dask
import dask.bag as db
from dask import delayed
from dask.base import clone_key
from dask.blockwise import Blockwise
from dask.graph_manipulation import bind, checkpoint, chunks, clone, wait_on
from dask.highlevelgraph import HighLevelGraph
from dask.tests.test_base import Tuple
from dask.utils_test import import_or_none

da = import_or_none("dask.array")
dd = import_or_none("dask.dataframe")
pd = import_or_none("pandas")
zarr = import_or_none("zarr")


class NodeCounter:
    def __init__(self):
        self.n = 0

    def __dask_tokenize__(self):
        return type(self), self.n

    def f(self, x):
        time.sleep(random.random() / 100)
        self.n += 1
        return x


def assert_no_common_keys(a, b, omit=None, *, layers: bool) -> None:
    dsk1 = a.__dask_graph__()
    dsk2 = b.__dask_graph__()

    if omit is not None:
        dsko = omit.__dask_graph__()
        assert not (dsk1.keys() - dsko.keys()) & dsk2.keys()
        assert not dsko.keys() - dsk1.keys()
        assert not dsko.keys() - dsk2.keys()
        if layers:
            assert not (dsk1.layers.keys() - dsko.layers.keys()) & dsk2.layers.keys()
            assert (
                not (dsk1.dependencies.keys() - dsko.dependencies.keys())
                & dsk2.dependencies.keys()
            )
            assert not dsko.layers.keys() - dsk1.layers.keys()
            assert not dsko.layers.keys() - dsk2.layers.keys()
            assert not dsko.dependencies.keys() - dsk1.dependencies.keys()
            assert not dsko.dependencies.keys() - dsk2.dependencies.keys()

    else:
        assert not dsk1.keys() & dsk2.keys()
        if layers:
            assert not dsk1.layers.keys() & dsk2.layers.keys()
            assert not dsk1.dependencies.keys() & dsk2.dependencies.keys()


def assert_did_not_materialize(cloned, orig):
    """Test that all layers of the original collection exist in the cloned collection
    too and that Blockwise layers have not been materialized
    """
    olayers = orig.__dask_graph__().layers
    clayers = cloned.__dask_graph__().layers
    for k, v in olayers.items():
        try:
            cv = clayers[k]
        except KeyError:
            cv = clayers[clone_key(k, 0)]
        if isinstance(v, Blockwise):
            assert not v.is_materialized()
            assert not cv.is_materialized()


# Generic valid keys
h1 = (1.2, "foo", ())
h2 = "h2"


def collections_with_node_counters():
    cnt = NodeCounter()
    df = pd.DataFrame({"x": list(range(10))})

    # Test two samples of all collections where applicable, one with multiple chunks
    # and one with a single chunk
    colls = [
        # dask.delayed
        delayed(cnt.f)("Hello 1"),  # 1 chunk
        # dask.array
        da.ones((10, 10), chunks=5).map_blocks(cnt.f),  # 4 chunks
        da.ones((1,), chunks=-1).map_blocks(cnt.f),  # 1 chunk
        # dask.bag
        db.from_sequence([1, 2], npartitions=2).map(cnt.f),  # 2 chunks
        db.from_sequence([1], npartitions=1).map(cnt.f),  # 1 chunk
        db.Item.from_delayed(delayed(cnt.f)("Hello 2")),  # 1 chunk
        # dask.dataframe
        dd.from_pandas(df, npartitions=2).map_partitions(cnt.f),  # 2 chunks
        dd.from_pandas(df, npartitions=1).map_partitions(cnt.f),  # 1 chunk
        dd.from_pandas(df["x"], npartitions=2).map_partitions(cnt.f),  # 2 chunks
        dd.from_pandas(df["x"], npartitions=1).map_partitions(cnt.f),  # 1 chunk
    ]
    cnt.n = 0
    return colls, cnt


def demo_tuples(layers: bool) -> tuple[Tuple, Tuple, NodeCounter]:
    cnt = NodeCounter()
    # Collections have multiple names
    dsk1 = HighLevelGraph(
        {"a": {("a", h1): (cnt.f, 1), ("a", h2): (cnt.f, 2)}, "b": {"b": (cnt.f, 3)}},
        {"a": set(), "b": set()},
    )
    dsk2 = HighLevelGraph(
        {"c": {"c": (cnt.f, 4)}, "d": {"d": (cnt.f, 5)}},
        {"c": set(), "d": set()},
    )
    if not layers:
        dsk1 = dsk1.to_dict()  # type: ignore
        dsk2 = dsk2.to_dict()  # type: ignore

    return Tuple(dsk1, list(dsk1)), Tuple(dsk2, list(dsk2)), cnt


@pytest.mark.parametrize("layers", [False, True])
def test_checkpoint(layers):
    t1, t2, cnt = demo_tuples(layers)
    cp = checkpoint(t1, {"x": [t2]})
    assert cp.compute(scheduler="sync") is None
    assert cnt.n == 5


@pytest.mark.skipif("not da or not dd")
def test_checkpoint_collections():
    colls, cnt = collections_with_node_counters()
    cp = checkpoint(*colls)
    cp.compute(scheduler="sync")
    assert cnt.n == 16


@pytest.mark.parametrize("layers", [False, True])
def test_wait_on_one(layers):
    t1, _, cnt = demo_tuples(layers)
    t1w = wait_on(t1)
    assert t1w.compute(scheduler="sync") == (1, 2, 3)
    assert cnt.n == 3


@pytest.mark.parametrize("layers", [False, True])
def test_wait_on_many(layers):
    t1, t2, cnt = demo_tuples(layers)
    out = wait_on(t1, {"x": [t2]})
    assert dask.compute(*out, scheduler="sync") == ((1, 2, 3), {"x": [(4, 5)]})
    assert cnt.n == 5


@pytest.mark.skipif("not da or not dd")
def test_wait_on_collections():
    dd = pytest.importorskip("dask.dataframe")
    if dd._dask_expr_enabled():
        pytest.skip("hlg doesn't make sense")
    colls, cnt = collections_with_node_counters()

    # Create a delayed that depends on a single one among all collections
    @delayed
    def f(x):
        pass

    colls2 = wait_on(*colls)
    f(colls2[0]).compute()
    assert cnt.n == 16

    # dask.delayed
    assert colls2[0].compute() == colls[0].compute()
    # dask.array
    da.utils.assert_eq(colls2[1], colls[1])
    da.utils.assert_eq(colls2[2], colls[2])
    # dask.bag
    db.utils.assert_eq(colls2[3], colls[3])
    db.utils.assert_eq(colls2[4], colls[4])
    db.utils.assert_eq(colls2[5], colls[5])
    # dask.dataframe
    dd.utils.assert_eq(colls2[6], colls[6])
    dd.utils.assert_eq(colls2[7], colls[7])
    dd.utils.assert_eq(colls2[8], colls[8])
    dd.utils.assert_eq(colls2[9], colls[9])


@pytest.mark.parametrize("layers", [False, True])
def test_clone(layers):
    dsk1 = {("a", h1): 1, ("a", h2): 2}
    dsk2 = {"b": (add, ("a", h1), ("a", h2))}
    dsk3 = {"c": 1, "d": 1}  # Multiple names
    if layers:
        dsk1 = HighLevelGraph.from_collections("a", dsk1)
        dsk2 = HighLevelGraph(
            {"a": dsk1, "b": dsk2}, dependencies={"a": set(), "b": {"a"}}
        )
        dsk3 = HighLevelGraph.from_collections("c", dsk3)
    else:
        dsk2.update(dsk1)

    t1 = Tuple(dsk1, [("a", h1), ("a", h2)])
    t2 = Tuple(dsk2, ["b"])
    t3 = Tuple(dsk3, ["c"])

    c1 = clone(t2, seed=1, assume_layers=layers)
    c2 = clone(t2, seed=1, assume_layers=layers)
    c3 = clone(t2, seed=2, assume_layers=layers)
    c4 = clone(c1, seed=1, assume_layers=layers)  # Clone of a clone has different keys
    c5 = clone(t2, assume_layers=layers)  # Random seed
    c6 = clone(t2, assume_layers=layers)  # Random seed
    c7 = clone(t2, omit=t1, seed=1, assume_layers=layers)

    assert c1.__dask_graph__() == c2.__dask_graph__()
    assert_no_common_keys(c1, t2, layers=layers)
    assert_no_common_keys(c1, c3, layers=layers)
    assert_no_common_keys(c1, c4, layers=layers)
    assert_no_common_keys(c1, c5, layers=layers)
    assert_no_common_keys(c5, c6, layers=layers)
    assert_no_common_keys(c7, t2, omit=t1, layers=layers)
    assert dask.compute(t2, c1, c2, c3, c4, c5, c6, c7) == ((3,),) * 8

    # Clone nested; some of the collections in omit are unrelated
    out = clone({"x": [t2]}, omit={"y": [t1, t3]}, assume_layers=layers)
    assert dask.compute(out) == ({"x": [(3,)]},)
    c8 = out["x"][0]
    assert_no_common_keys(c8, t2, omit=t1, layers=layers)
    assert_no_common_keys(c8, t3, layers=layers)


@pytest.mark.skipif("not da")
@pytest.mark.parametrize(
    "literal",
    [
        1,
        (1,),
        [1],
        {1: 1},
        {1},
    ],
)
def test_blockwise_clone_with_literals(literal):
    """https://github.com/dask/dask/issues/8978

    clone() on the result of a dask.array.blockwise operation with a (iterable) literal
    argument
    """
    arr = da.ones(10, chunks=1)

    def noop(arr, lit):
        return arr

    blk = da.blockwise(noop, "x", arr, "x", literal, None)

    cln = clone(blk)

    assert_no_common_keys(blk, cln, layers=True)
    da.utils.assert_eq(blk, cln)


@pytest.mark.skipif("not da or not zarr")
def test_blockwise_clone_with_no_indices():
    """https://github.com/dask/dask/issues/9621

    clone() on a Blockwise layer on top of a dependency layer with no indices
    """
    blk = da.from_zarr(zarr.ones(10))
    # This use case leverages the current implementation details of from_array when the
    # input is neither a numpy.ndarray nor a list. If it were to change in the future,
    # please find another way to create a use case that satisfies these assertions.
    assert isinstance(blk.dask.layers[blk.name], Blockwise)
    assert any(isinstance(k, str) for k in blk.dask)

    cln = clone(blk)
    assert_no_common_keys(blk, cln, layers=True)
    da.utils.assert_eq(blk, cln)


@pytest.mark.parametrize("layers", [False, True])
def test_bind(layers):
    dsk1 = {("a-1", h1): 1, ("a-1", h2): 2}
    dsk2 = {"b-1": (add, ("a-1", h1), ("a-1", h2))}
    dsk3 = {"c-1": "b-1"}
    cnt = NodeCounter()
    dsk4 = {("d-1", h1): (cnt.f, 1), ("d-1", h2): (cnt.f, 2)}
    dsk4b = {"e": (cnt.f, 3)}

    if layers:
        dsk1 = HighLevelGraph({"a-1": dsk1}, {"a-1": set()})
        dsk2 = HighLevelGraph(
            {"a-1": dsk1, "b-1": dsk2}, {"a-1": set(), "b-1": {"a-1"}}
        )
        dsk3 = HighLevelGraph(
            {"a-1": dsk1, "b-1": dsk2, "c-1": dsk3},
            {"a-1": set(), "b-1": {"a-1"}, "c-1": {"b-1"}},
        )
        dsk4 = HighLevelGraph({"d-1": dsk4, "e": dsk4b}, {"d-1": set(), "e": set()})
    else:
        dsk2.update(dsk1)
        dsk3.update(dsk2)
        dsk4.update(dsk4b)

    # t1 = Tuple(dsk1, [("a", h1), ("a", h2)])
    t2 = Tuple(dsk2, ["b-1"])
    t3 = Tuple(dsk3, ["c-1"])
    t4 = Tuple(dsk4, [("d-1", h1), ("d-1", h2), "e"])  # Multiple names

    bound1 = bind(t3, t4, seed=1, assume_layers=layers)
    cloned_a_name = clone_key("a-1", seed=1)
    assert bound1.__dask_graph__()[cloned_a_name, h1][0] is chunks.bind
    assert bound1.__dask_graph__()[cloned_a_name, h2][0] is chunks.bind
    assert bound1.compute() == (3,)
    assert cnt.n == 3

    bound2 = bind(t3, t4, omit=t2, seed=1, assume_layers=layers)
    cloned_c_name = clone_key("c-1", seed=1)
    assert bound2.__dask_graph__()[cloned_c_name][0] is chunks.bind
    assert bound2.compute() == (3,)
    assert cnt.n == 6

    bound3 = bind(t4, t3, seed=1, assume_layers=layers)
    cloned_d_name = clone_key("d-1", seed=1)
    cloned_e_name = clone_key("e", seed=1)
    assert bound3.__dask_graph__()[cloned_d_name, h1][0] is chunks.bind
    assert bound3.__dask_graph__()[cloned_d_name, h2][0] is chunks.bind
    assert bound3.__dask_graph__()[cloned_e_name][0] is chunks.bind
    assert bound3.compute() == (1, 2, 3)
    assert cnt.n == 9


@pytest.mark.skipif("not da or not dd")
@pytest.mark.parametrize("func", [bind, clone])
def test_bind_clone_collections(func):
    if dd._dask_expr_enabled():
        pytest.skip("not supported")

    @delayed
    def double(x):
        return x * 2

    # dask.delayed
    d1 = double(2)
    d2 = double(d1)
    # dask.array
    a1 = da.ones((10, 10), chunks=5)
    a2 = a1 + 1
    a3 = a2.T
    # dask.bag
    b1 = db.from_sequence([1, 2], npartitions=2)
    # b1's tasks are not callable, so we need an extra step to properly test bind
    b2 = b1.map(lambda x: x * 2)
    b3 = b2.map(lambda x: x + 1)
    b4 = b3.min()
    # dask.dataframe
    df = pd.DataFrame({"x": list(range(10))})
    ddf1 = dd.from_pandas(df, npartitions=2)
    # ddf1's tasks are not callable, so we need an extra step to properly test bind
    ddf2 = ddf1.map_partitions(lambda x: x * 2)
    ddf3 = ddf2.map_partitions(lambda x: x + 1)
    ddf4 = ddf3["x"]  # dd.Series
    ddf5 = ddf4.min()  # dd.Scalar

    cnt = NodeCounter()
    if func is bind:
        parent = da.ones((10, 10), chunks=5).map_blocks(cnt.f)
        cnt.n = 0
        d2c, a3c, b3c, b4c, ddf3c, ddf4c, ddf5c = bind(
            children=(d2, a3, b3, b4, ddf3, ddf4, ddf5),
            parents=parent,
            omit=(d1, a1, b2, ddf2),
            seed=0,
        )
    else:
        d2c, a3c, b3c, b4c, ddf3c, ddf4c, ddf5c = clone(
            d2,
            a3,
            b3,
            b4,
            ddf3,
            ddf4,
            ddf5,
            omit=(d1, a1, b2, ddf2),
            seed=0,
        )

    assert_did_not_materialize(d2c, d2)
    assert_did_not_materialize(a3c, a3)
    assert_did_not_materialize(b3c, b3)
    assert_did_not_materialize(b4c, b4)
    assert_did_not_materialize(ddf3c, ddf3)
    assert_did_not_materialize(ddf4c, ddf4)
    assert_did_not_materialize(ddf5c, ddf5)

    assert_no_common_keys(d2c, d2, omit=d1, layers=True)
    assert_no_common_keys(a3c, a3, omit=a1, layers=True)
    assert_no_common_keys(b3c, b3, omit=b2, layers=True)
    assert_no_common_keys(ddf3c, ddf3, omit=ddf2, layers=True)
    assert_no_common_keys(ddf4c, ddf4, omit=ddf2, layers=True)
    assert_no_common_keys(ddf5c, ddf5, omit=ddf2, layers=True)

    assert d2.compute() == d2c.compute()
    assert cnt.n == 4 or func is clone
    da.utils.assert_eq(a3c, a3)
    assert cnt.n == 8 or func is clone
    db.utils.assert_eq(b3c, b3)
    assert cnt.n == 12 or func is clone
    db.utils.assert_eq(b4c, b4)
    assert cnt.n == 16 or func is clone
    dd.utils.assert_eq(ddf3c, ddf3)
    assert cnt.n == 24 or func is clone  # dd.utils.assert_eq calls compute() twice
    dd.utils.assert_eq(ddf4c, ddf4)
    assert cnt.n == 32 or func is clone  # dd.utils.assert_eq calls compute() twice
    dd.utils.assert_eq(ddf5c, ddf5)
    assert cnt.n == 36 or func is clone


@pytest.mark.parametrize(
    "split_every,nkeys",
    [
        (2, 299),
        (3, 250),
        (8, 215),
        (None, 215),  # default is 8
        (8.1, 215),
        (1e9, 201),
        (False, 201),
    ],
)
def test_split_every(split_every, nkeys):
    dsk = {("a", i): i for i in range(100)}
    t1 = Tuple(dsk, list(dsk))
    c = checkpoint(t1, split_every=split_every)
    assert len(c.__dask_graph__()) == nkeys
    assert c.compute(scheduler="sync") is None

    t2 = wait_on(t1, split_every=split_every)
    assert len(t2.__dask_graph__()) == nkeys + 100
    assert t2.compute(scheduler="sync") == tuple(range(100))

    dsk3 = {"b": 1, "c": 2}
    t3 = Tuple(dsk3, list(dsk3))
    t4 = bind(t3, t1, split_every=split_every, assume_layers=False)
    assert len(t4.__dask_graph__()) == nkeys + 2
    assert t4.compute(scheduler="sync") == (1, 2)


def test_split_every_invalid():
    t = Tuple({"a": 1, "b": 2}, ["a", "b"])
    with pytest.raises(ValueError):
        checkpoint(t, split_every=1)
    with pytest.raises(ValueError):
        checkpoint(t, split_every=1.9)
    with pytest.raises(ValueError):
        checkpoint(t, split_every=0)  # Not to be confused with False or None
    with pytest.raises(ValueError):
        checkpoint(t, split_every=-2)
    with pytest.raises(TypeError):
        checkpoint(t, split_every={0: 2})  # This is legal for dask.array but not here
