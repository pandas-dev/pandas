from __future__ import annotations

import dataclasses
import inspect
import os
import subprocess
import sys
import time
from collections import OrderedDict
from concurrent.futures import Executor
from operator import add, mul
from typing import NamedTuple

import pytest
from tlz import merge, partial

import dask
import dask.bag as db
from dask.base import (
    DaskMethodsMixin,
    clone_key,
    compute,
    compute_as_if_collection,
    get_collection_names,
    get_name_from_key,
    get_scheduler,
    is_dask_collection,
    named_schedulers,
    optimize,
    persist,
    replace_name_in_key,
    tokenize,
    unpack_collections,
    visualize,
)
from dask.core import validate_key
from dask.delayed import Delayed, delayed
from dask.diagnostics import Profiler
from dask.highlevelgraph import HighLevelGraph
from dask.utils import key_split, tmpdir, tmpfile
from dask.utils_test import dec, import_or_none, inc

da = import_or_none("dask.array")
dd = import_or_none("dask.dataframe")
np = import_or_none("numpy")
pd = import_or_none("pandas")

# Arbitrary dask keys
h1 = (1.2, "foo", (3,))
h2 = "h2"


def test_is_dask_collection():
    class DummyCollection:
        def __init__(self, dsk):
            self.dask = dsk

        def __dask_graph__(self):
            return self.dask

    x = delayed(1) + 2
    assert is_dask_collection(x)
    assert not is_dask_collection(2)
    assert is_dask_collection(DummyCollection({}))
    assert not is_dask_collection(DummyCollection)


def test_is_dask_collection_dask_expr():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    df = pd.Series([1, 2, 3])
    ddf = dd.from_pandas(df)
    assert not is_dask_collection(df)
    assert is_dask_collection(ddf)


def test_is_dask_collection_dask_expr_does_not_materialize():
    pytest.importorskip("pandas")
    dx = pytest.importorskip("dask.dataframe.dask_expr")

    class DoNotMaterialize(dask._expr.Expr):
        @property
        def _meta(self):
            return 0

        def __dask_keys__(self):
            assert False, "must not reach"

        def __dask_graph__(self):
            assert False, "must not reach"

        def optimize(self, fuse=False):
            assert False, "must not reach"

    coll = dx._collection.new_collection(DoNotMaterialize())

    with pytest.raises(AssertionError, match="must not reach"):
        coll.__dask_keys__()
    with pytest.raises(AssertionError, match="must not reach"):
        coll.__dask_graph__()

    assert is_dask_collection(coll)


def test_unpack_collections():
    @dataclasses.dataclass
    class ADataClass:
        a: int  # type: ignore[annotation-unchecked]

    class ANamedTuple(NamedTuple):
        a: int  # type: ignore[annotation-unchecked]

    a = delayed(1) + 5
    b = a + 1
    c = a + 2

    def build(a, b, c, iterator):
        t = (
            a,
            b,  # Top-level collections
            {
                "a": a,  # dict
                a: b,  # collections as keys
                "b": [1, 2, [b]],  # list
                "c": 10,  # other builtins pass through unchanged
                "d": (c, 2),  # tuple
                "e": {a, 2, 3},  # set
                "f": OrderedDict([("a", a)]),
                "g": ADataClass(a=a),  # dataclass instance
                "h": (ADataClass, a),  # dataclass constructor
                "i": ANamedTuple(a=a),  # namedtuple instance
            },  # OrderedDict
            iterator,
        )  # Iterator

        return t

    args = build(a, b, c, (i for i in [a, b, c]))

    collections, repack = unpack_collections(*args)
    assert len(collections) == 3

    # Replace collections with `'~a'` strings
    result = repack(["~a", "~b", "~c"])
    sol = build("~a", "~b", "~c", ["~a", "~b", "~c"])
    assert result == sol

    # traverse=False
    collections, repack = unpack_collections(*args, traverse=False)
    assert len(collections) == 2  # just a and b
    assert repack(collections) == args

    # No collections
    collections, repack = unpack_collections(1, 2, {"a": 3})
    assert not collections
    assert repack(collections) == (1, 2, {"a": 3})

    # Result that looks like a task
    def fail(*args):
        raise ValueError("Shouldn't have been called")  # pragma: nocover

    collections, repack = unpack_collections(
        a, (fail, 1), [(fail, 2, 3)], traverse=False
    )
    repack(collections)  # Smoketest task literals
    repack([(fail, 1)])  # Smoketest results that look like tasks


def test_get_collection_names():
    class DummyCollection:
        def __init__(self, dsk, keys):
            self.dask = dsk
            self.keys = keys

        def __dask_graph__(self):
            return self.dask

        def __dask_keys__(self):
            return self.keys

    with pytest.raises(TypeError):
        get_collection_names(object())
    # Keys must either be a string or a tuple where the first element is a string
    with pytest.raises(TypeError):
        get_collection_names(DummyCollection({1: 2}, [1]))
    with pytest.raises(TypeError):
        get_collection_names(DummyCollection({(): 1}, [()]))
    with pytest.raises(TypeError):
        get_collection_names(DummyCollection({(1,): 1}, [(1,)]))

    assert get_collection_names(DummyCollection({}, [])) == set()

    # __dask_keys__() returns a nested list
    assert get_collection_names(
        DummyCollection(
            {("a-1", h1): 1, ("a-1", h2): 2, "b-2": 3, "c": 4},
            [[[("a-1", h1), ("a-1", h2), "b-2", "c"]]],
        )
    ) == {"a-1", "b-2", "c"}


def test_get_name_from_key():
    assert get_name_from_key("foo") == "foo"
    assert get_name_from_key("foo-123"), "foo-123"
    assert get_name_from_key(("foo-123", h1, h2)) == "foo-123"
    with pytest.raises(TypeError):
        get_name_from_key(1)
    with pytest.raises(TypeError):
        get_name_from_key(())
    with pytest.raises(TypeError):
        get_name_from_key((1,))


def test_replace_name_in_keys():
    assert replace_name_in_key("foo", {}) == "foo"
    assert replace_name_in_key("foo", {"bar": "baz"}) == "foo"
    assert replace_name_in_key("foo", {"foo": "bar", "baz": "asd"}) == "bar"
    assert replace_name_in_key("foo-123", {"foo-123": "bar-456"}) == "bar-456"
    assert replace_name_in_key(("foo-123", h1, h2), {"foo-123": "bar"}) == (
        "bar",
        h1,
        h2,
    )
    with pytest.raises(TypeError):
        replace_name_in_key(1, {})
    with pytest.raises(TypeError):
        replace_name_in_key((), {})
    with pytest.raises(TypeError):
        replace_name_in_key((1,), {})


class Tuple(DaskMethodsMixin):
    __slots__ = ("_dask", "_keys")
    __dask_scheduler__ = staticmethod(dask.threaded.get)
    __dask_optimize__ = None

    def __init__(self, dsk, keys):
        self._dask = dsk
        self._keys = keys

    def __add__(self, other):
        if not isinstance(other, Tuple):
            return NotImplemented  # pragma: nocover
        return Tuple(merge(self._dask, other._dask), self._keys + other._keys)

    def __dask_graph__(self):
        return self._dask

    def __dask_keys__(self):
        return self._keys

    def __dask_layers__(self):
        return tuple(get_collection_names(self))

    def __dask_tokenize__(self):
        return self._keys

    def __dask_postcompute__(self):
        return tuple, ()

    def __dask_postpersist__(self):
        return Tuple._rebuild, (self._keys,)

    @staticmethod
    def _rebuild(dsk, keys, *, rename=None):
        if rename:
            keys = [replace_name_in_key(key, rename) for key in keys]
        return Tuple(dsk, keys)


def test_custom_collection():
    dsk = {("x", h1): 1, ("x", h2): 2}
    dsk2 = {("y", h1): (add, ("x", h1), ("x", h2)), ("y", h2): (add, ("y", h1), 1)}
    dsk2.update(dsk)
    dsk3 = {"z": (add, ("y", h1), ("y", h2))}
    dsk3.update(dsk2)

    w = Tuple({}, [])  # A collection can have no keys at all
    x = Tuple(dsk, [("x", h1), ("x", h2)])
    y = Tuple(dsk2, [("y", h1), ("y", h2)])
    z = Tuple(dsk3, ["z"])
    # Collection with multiple names
    t = w + x + y + z

    # __slots__ defined on base mixin class propagates
    with pytest.raises(AttributeError):
        x.foo = 1

    # is_dask_collection
    assert is_dask_collection(w)
    assert is_dask_collection(x)
    assert is_dask_collection(y)
    assert is_dask_collection(z)
    assert is_dask_collection(t)

    # tokenize
    assert tokenize(w) == tokenize(w)
    assert tokenize(x) == tokenize(x)
    assert tokenize(y) == tokenize(y)
    assert tokenize(z) == tokenize(z)
    assert tokenize(t) == tokenize(t)
    # All tokens are unique
    assert len({tokenize(coll) for coll in (w, x, y, z, t)}) == 5

    # get_collection_names
    assert get_collection_names(w) == set()
    assert get_collection_names(x) == {"x"}
    assert get_collection_names(y) == {"y"}
    assert get_collection_names(z) == {"z"}
    assert get_collection_names(t) == {"x", "y", "z"}

    # compute
    assert w.compute() == ()
    assert x.compute() == (1, 2)
    assert y.compute() == (3, 4)
    assert z.compute() == (7,)
    assert dask.compute(w, [{"x": x}, y, z]) == ((), [{"x": (1, 2)}, (3, 4), (7,)])
    assert t.compute() == (1, 2, 3, 4, 7)

    # persist
    t2 = t.persist()
    assert isinstance(t2, Tuple)
    assert t2._keys == t._keys
    assert sorted(t2._dask.values()) == [1, 2, 3, 4, 7]
    assert t2.compute() == (1, 2, 3, 4, 7)

    w2, x2, y2, z2 = dask.persist(w, x, y, z)
    assert y2._keys == y._keys
    assert y2._dask == {("y", h1): 3, ("y", h2): 4}
    assert y2.compute() == (3, 4)

    t3 = x2 + y2 + z2
    assert t3.compute() == (1, 2, 3, 4, 7)

    # __dask_postpersist__ with name change
    rebuild, args = w.__dask_postpersist__()
    w3 = rebuild({}, *args, rename={"w": "w3"})
    assert w3.compute() == ()

    rebuild, args = x.__dask_postpersist__()
    x3 = rebuild({("x3", h1): 10, ("x3", h2): 20}, *args, rename={"x": "x3"})
    assert x3.compute() == (10, 20)

    rebuild, args = z.__dask_postpersist__()
    z3 = rebuild({"z3": 70}, *args, rename={"z": "z3"})
    assert z3.compute() == (70,)


def test_compute_no_opt():
    # Bag does `fuse` by default. Test that with `optimize_graph=False` that
    # doesn't get called. We check this by using a callback to track the keys
    # that are computed.
    from dask.callbacks import Callback

    b = db.from_sequence(range(100), npartitions=4)
    add1 = partial(add, 1)
    mul2 = partial(mul, 2)
    o = b.map(add1).map(mul2)
    # Check that with the kwarg, the optimization doesn't happen
    keys = []
    with Callback(pretask=lambda key, *args: keys.append(key)):
        o.compute(scheduler="single-threaded", optimize_graph=False)
    assert len([k for k in keys if "mul" in k[0]]) == 4
    assert len([k for k in keys if "add" in k[0]]) == 4
    # Check that without the kwarg, the optimization does happen
    keys = []
    with Callback(pretask=lambda key, *args: keys.append(key)):
        o.compute(scheduler="single-threaded")
    # Names of fused tasks have been merged, and the original key is an alias.
    # Otherwise, the lengths below would be 4 and 0.
    assert len([k for k in keys if "mul" in k[0]]) == 8
    assert len([k for k in keys if "add" in k[0]]) == 4
    assert len([k for k in keys if "add-from_sequence-mul" in k[0]]) == 4


@pytest.mark.skipif("not da")
def test_compute_array():
    arr = np.arange(100).reshape((10, 10))
    darr = da.from_array(arr, chunks=(5, 5))
    darr1 = darr + 1
    darr2 = darr + 2
    out1, out2 = compute(darr1, darr2)
    assert np.allclose(out1, arr + 1)
    assert np.allclose(out2, arr + 2)


@pytest.mark.skipif("not da")
def test_persist_array():
    from dask.array.utils import assert_eq

    arr = np.arange(100).reshape((10, 10))
    x = da.from_array(arr, chunks=(5, 5))
    x = (x + 1) - x.mean(axis=0)
    y = x.persist()

    assert_eq(x, y)
    assert set(y.dask).issubset(x.dask)
    assert len(y.dask) == y.npartitions


@pytest.mark.skipif("not da")
def test_persist_array_rename():
    a = da.zeros(4, dtype=int, chunks=2)
    rebuild, args = a.__dask_postpersist__()
    dsk = {("b", 0): np.array([1, 2]), ("b", 1): np.array([3, 4])}
    b = rebuild(dsk, *args, rename={a.name: "b"})
    assert isinstance(b, da.Array)
    assert b.name == "b"
    assert b.__dask_keys__() == [("b", 0), ("b", 1)]
    da.utils.assert_eq(b, [1, 2, 3, 4])


@pytest.mark.skipif("not dd")
def test_compute_dataframe():
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [5, 5, 3, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf1 = ddf.a + 1
    ddf2 = ddf.a + ddf.b
    out1, out2 = compute(ddf1, ddf2)
    dd.utils.assert_eq(out1, df.a + 1)
    dd.utils.assert_eq(out2, df.a + df.b)


@pytest.mark.skipif("not dd")
def test_persist_dataframe():
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [5, 6, 7, 8]})
    ddf1 = dd.from_pandas(df, npartitions=2) * 2
    assert len(ddf1.__dask_graph__()) == 4
    ddf2 = ddf1.persist()
    assert isinstance(ddf2, dd.DataFrame)
    assert len(ddf2.__dask_graph__()) == 4
    dd.utils.assert_eq(ddf2, ddf1)


@pytest.mark.skipif("not dd")
def test_persist_series():
    ds = pd.Series([1, 2, 3, 4])
    dds1 = dd.from_pandas(ds, npartitions=2) * 2
    assert len(dds1.__dask_graph__()) == 4
    dds2 = dds1.persist()
    assert isinstance(dds2, dd.Series)
    assert len(dds2.__dask_graph__()) == 4
    dd.utils.assert_eq(dds2, dds1)


@pytest.mark.skipif("not dd")
def test_persist_scalar():
    import dask.dataframe as dd

    ds = pd.Series([1, 2, 3, 4])
    dds1 = dd.from_pandas(ds, npartitions=2).min()
    assert len(dds1.__dask_graph__()) == 5
    dds2 = dds1.persist()
    assert isinstance(dds2, dd.Scalar)
    assert len(dds2.__dask_graph__()) == 2
    dd.utils.assert_eq(dds2, dds1)


@pytest.mark.skipif("not dd or not da")
def test_compute_array_dataframe():
    arr = np.arange(100).reshape((10, 10))
    darr = da.from_array(arr, chunks=(5, 5)) + 1
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [5, 5, 3, 3]})
    ddf = dd.from_pandas(df, npartitions=2).a + 2
    with pytest.warns(UserWarning, match="mixed.*materialize"):
        arr_out, df_out = compute(darr, ddf)
    assert np.allclose(arr_out, arr + 1)
    dd.utils.assert_eq(df_out, df.a + 2)


@pytest.mark.skipif("not dd")
def test_compute_dataframe_valid_unicode_in_bytes():
    df = pd.DataFrame(data=np.random.random((3, 1)), columns=["ö".encode()])
    dd.from_pandas(df, npartitions=4)


@pytest.mark.skipif("not dd")
def test_compute_dataframe_invalid_unicode():
    # see https://github.com/dask/dask/issues/2713
    df = pd.DataFrame(data=np.random.random((3, 1)), columns=["\ud83d"])
    dd.from_pandas(df, npartitions=4)


@pytest.mark.skipif("not da")
def test_compute_array_bag():
    x = da.arange(5, chunks=2)
    b = db.from_sequence([1, 2, 3])

    pytest.raises(ValueError, lambda: compute(x, b))

    xx, bb = compute(x, b, scheduler="single-threaded")
    assert np.allclose(xx, np.arange(5))
    assert bb == [1, 2, 3]


@pytest.mark.skipif("not da")
def test_compute_with_literal():
    x = da.arange(5, chunks=2)
    y = 10

    xx, yy = compute(x, y)
    assert (xx == x.compute()).all()
    assert yy == y

    assert compute(5) == (5,)


def test_compute_nested():
    a = delayed(1) + 5
    b = a + 1
    c = a + 2
    assert compute({"a": a, "b": [1, 2, b]}, (c, 2)) == (
        {"a": 6, "b": [1, 2, 7]},
        (8, 2),
    )

    res = compute([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1] == 8


@pytest.mark.skipif("not da")
@pytest.mark.skipif(
    bool(sys.flags.optimize), reason="graphviz exception with Python -OO flag"
)
@pytest.mark.xfail(
    sys.platform == "win32",
    reason="graphviz/pango on conda-forge currently broken for windows",
    strict=False,
)
def test_visualize():
    pytest.importorskip("graphviz")
    pytest.importorskip("ipycytoscape")
    with tmpdir() as d:
        x = da.arange(5, chunks=2)
        x.visualize(filename=os.path.join(d, "mydask"))
        assert os.path.exists(os.path.join(d, "mydask.png"))

        x.visualize(filename=os.path.join(d, "mydask.pdf"))
        assert os.path.exists(os.path.join(d, "mydask.pdf"))

        visualize(x, 1, 2, filename=os.path.join(d, "mydask.png"))
        assert os.path.exists(os.path.join(d, "mydask.png"))

        dsk = {"a": 1, "b": (add, "a", 2), "c": (mul, "a", 1)}
        visualize(x, dsk, filename=os.path.join(d, "mydask.png"))
        assert os.path.exists(os.path.join(d, "mydask.png"))

        x = Tuple(dsk, ["a", "b", "c"])
        visualize(x, filename=os.path.join(d, "mydask.png"))
        assert os.path.exists(os.path.join(d, "mydask.png"))

        x = Tuple(dsk, ["a", "b", "c"])
        visualize(x, filename=os.path.join(d, "cyt"), engine="cytoscape")
        assert os.path.exists(os.path.join(d, "cyt.html"))

        visualize(x, filename=os.path.join(d, "cyt2.html"), engine="ipycytoscape")
        assert os.path.exists(os.path.join(d, "cyt2.html"))

        with dask.config.set(visualization__engine="cytoscape"):
            visualize(x, filename=os.path.join(d, "cyt3.html"))
            assert os.path.exists(os.path.join(d, "cyt3.html"))

        with pytest.raises(ValueError, match="not-real"):
            visualize(x, engine="not-real")

        # To see if visualize() works when the filename parameter is set to None
        # If the function raises an error, the test will fail
        x.visualize(filename=None)


@pytest.mark.skipif("not da")
@pytest.mark.skipif(
    bool(sys.flags.optimize), reason="graphviz exception with Python -OO flag"
)
def test_visualize_highlevelgraph():
    graphviz = pytest.importorskip("graphviz")
    with tmpdir() as d:
        x = da.arange(5, chunks=2)
        viz = x.dask.visualize(filename=os.path.join(d, "mydask.png"))
        # check visualization will automatically render in the jupyter notebook
        assert isinstance(viz, graphviz.Digraph)


@pytest.mark.skipif("not da")
@pytest.mark.skipif(
    bool(sys.flags.optimize), reason="graphviz exception with Python -OO flag"
)
def test_visualize_order():
    pytest.importorskip("graphviz")
    pytest.importorskip("matplotlib.pyplot")
    x = da.arange(5, chunks=2)
    with tmpfile(extension="dot") as fn:
        x.visualize(color="order", filename=fn, cmap="RdBu")
        with open(fn) as f:
            text = f.read()
        assert 'color="#' in text


def inc_to_dec(dsk, keys):
    dsk = dict(dsk)
    for key in dsk:
        if dsk[key][0] == inc:
            dsk[key] = (dec,) + dsk[key][1:]
    return dsk


def test_optimize():
    x = dask.delayed(inc)(1)
    y = dask.delayed(inc)(x)
    z = x + y

    x2, y2, z2, constant = optimize(x, y, z, 1)
    assert constant == 1

    # Same graphs for each
    dsk = dict(x2.dask)
    assert dict(y2.dask) == dsk
    assert dict(z2.dask) == dsk

    # Computationally equivalent
    assert dask.compute(x2, y2, z2) == dask.compute(x, y, z)


def test_optimize_nested():
    a = dask.delayed(inc)(1)
    b = dask.delayed(inc)(a)
    c = a + b

    result = optimize({"a": a, "b": [1, 2, b]}, (c, 2))

    a2 = result[0]["a"]
    b2 = result[0]["b"][2]
    c2 = result[1][0]

    assert isinstance(a2, Delayed)
    assert isinstance(b2, Delayed)
    assert isinstance(c2, Delayed)
    assert dict(a2.dask) == dict(b2.dask) == dict(c2.dask)
    assert compute(*result) == ({"a": 2, "b": [1, 2, 3]}, (5, 2))

    res = optimize([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1].compute() == 5


def test_default_imports():
    """
    Startup time: `import dask` should not import too many modules.
    """
    code = """if 1:
        import dask
        import sys

        print(sorted(sys.modules))
        """

    out = subprocess.check_output([sys.executable, "-c", code])
    modules = set(eval(out.decode()))
    assert "dask" in modules
    blacklist = [
        "dask.array",
        "dask.dataframe",
        "numpy",
        "pandas",
        "partd",
        "s3fs",
        "distributed",
    ]
    for mod in blacklist:
        assert mod not in modules


def test_persist_literals():
    assert persist(1, 2, 3) == (1, 2, 3)


def test_persist_nested():
    a = delayed(1) + 5
    b = a + 1
    c = a + 2
    result = persist({"a": a, "b": [1, 2, b]}, (c, 2))
    assert isinstance(result[0]["a"], Delayed)
    assert isinstance(result[0]["b"][2], Delayed)
    assert isinstance(result[1][0], Delayed)
    assert compute(*result) == ({"a": 6, "b": [1, 2, 7]}, (8, 2))

    res = persist([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1].compute() == 8


def test_persist_delayed():
    x1 = delayed(1)
    x2 = delayed(inc)(x1)
    x3 = delayed(inc)(x2)
    (xx,) = persist(x3)
    assert isinstance(xx, Delayed)
    assert xx.key == x3.key
    assert len(xx.dask) == 1

    assert x3.compute() == xx.compute()


@pytest.mark.parametrize("key", ["a", ("a-123", h1)])
def test_persist_delayed_custom_key(key):
    d = Delayed(key, {key: "b", "b": 1})
    assert d.compute() == 1
    dp = d.persist()
    assert dp.compute() == 1
    assert dp.key == key
    assert dict(dp.dask) == {key: 1}


@pytest.mark.parametrize(
    "key,rename,new_key",
    [
        ("a", {}, "a"),
        ("a", {"c": "d"}, "a"),
        ("a", {"a": "b"}, "b"),
        (("a-123", h1), {"a-123": "b-123"}, ("b-123", h1)),
    ],
)
def test_persist_delayed_rename(key, rename, new_key):
    d = Delayed(key, {key: 1})
    assert d.compute() == 1
    rebuild, args = d.__dask_postpersist__()
    dp = rebuild({new_key: 2}, *args, rename=rename)
    assert dp.compute() == 2
    assert dp.key == new_key
    assert dict(dp.dask) == {new_key: 2}


def test_persist_delayedleaf():
    x = delayed(1)
    (xx,) = persist(x)
    assert isinstance(xx, Delayed)
    assert xx.compute() == 1


def test_persist_delayedattr():
    class C:
        x = 1

    x = delayed(C).x
    (xx,) = persist(x)
    assert isinstance(xx, Delayed)
    assert xx.compute() == 1


@pytest.mark.skipif("not da")
def test_persist_array_bag():
    x = da.arange(5, chunks=2) + 1
    b = db.from_sequence([1, 2, 3]).map(inc)

    with pytest.raises(ValueError):
        persist(x, b)

    xx, bb = persist(x, b, scheduler="single-threaded")

    assert isinstance(xx, da.Array)
    assert isinstance(bb, db.Bag)

    assert xx.name == x.name
    assert bb.name == b.name
    assert len(xx.dask) == xx.npartitions < len(x.dask)
    assert len(bb.dask) == bb.npartitions < len(b.dask)

    assert np.allclose(x, xx)
    assert list(b) == list(bb)


def test_persist_bag():
    a = db.from_sequence([1, 2, 3], npartitions=2).map(lambda x: x * 2)
    assert len(a.__dask_graph__()) == 4
    b = a.persist(scheduler="sync")
    assert isinstance(b, db.Bag)
    assert len(b.__dask_graph__()) == 2
    db.utils.assert_eq(a, b)


def test_persist_item():
    a = db.from_sequence([1, 2, 3], npartitions=2).map(lambda x: x * 2).min()
    assert len(a.__dask_graph__()) == 7
    b = a.persist(scheduler="sync")
    assert isinstance(b, db.Item)
    assert len(b.__dask_graph__()) == 1
    db.utils.assert_eq(a, b)


def test_persist_bag_rename():
    a = db.from_sequence([1, 2, 3], npartitions=2)
    rebuild, args = a.__dask_postpersist__()
    dsk = {("b", 0): [4], ("b", 1): [5, 6]}
    b = rebuild(dsk, *args, rename={a.name: "b"})
    assert isinstance(b, db.Bag)
    assert b.name == "b"
    assert b.__dask_keys__() == [("b", 0), ("b", 1)]
    db.utils.assert_eq(b, [4, 5, 6])


def test_persist_item_change_name():
    a = db.from_sequence([1, 2, 3]).min()
    rebuild, args = a.__dask_postpersist__()
    b = rebuild({"x": 4}, *args, rename={a.name: "x"})
    assert isinstance(b, db.Item)
    assert b.__dask_keys__() == ["x"]
    db.utils.assert_eq(b, 4)


def test_optimize_globals():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")

    x = da.ones(10, chunks=(5,))

    def optimize_double(dsk, keys):
        return {
            k: (mul, 2, v) if not key_split(k).startswith("finalize") else v
            for k, v in dsk.items()
        }

    from dask.array.utils import assert_eq

    assert_eq(x + 1, np.ones(10) + 1)

    with dask.config.set(array_optimize=optimize_double):
        assert_eq(x + 1, (np.ones(10) * 2 + 1) * 2, check_chunks=False)

    assert_eq(x + 1, np.ones(10) + 1)

    b = db.range(10, npartitions=2)

    with dask.config.set(array_optimize=optimize_double):
        xx, bb = dask.compute(x + 1, b.map(inc), scheduler="single-threaded")
        assert_eq(xx, (np.ones(10) * 2 + 1) * 2)


def test_optimize_None():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")

    x = da.ones(10, chunks=(5,))
    y = x[:9][1:8][::2] + 1  # normally these slices would be fused

    def my_get(dsk, keys):
        # but they aren't. +1 for the finalize task
        assert len(dsk.__dask_graph__()) == len(y.dask) + 1
        return dask.get(dsk, keys)

    with dask.config.set(array_optimize=None, scheduler=my_get):
        y.compute()


def test_scheduler_keyword():
    def schedule(dsk, keys, **kwargs):
        return [123]

    named_schedulers["foo"] = schedule

    x = delayed(inc)(1)

    try:
        assert x.compute() == 2
        assert x.compute(scheduler="foo") == 123

        with dask.config.set(scheduler="foo"):
            assert x.compute() == 123
        assert x.compute() == 2

        with dask.config.set(scheduler="foo"):
            assert x.compute(scheduler="threads") == 2
    finally:
        del named_schedulers["foo"]


def test_raise_get_keyword():
    x = delayed(inc)(1)

    with pytest.raises(TypeError) as info:
        x.compute(get=dask.get)

    assert "scheduler=" in str(info.value)


class MyExecutor(Executor):
    _max_workers = None


def test_get_scheduler():
    assert get_scheduler() is None
    assert get_scheduler(scheduler=dask.local.get_sync) is dask.local.get_sync
    assert get_scheduler(scheduler="threads") is dask.threaded.get
    assert get_scheduler(scheduler="sync") is dask.local.get_sync
    assert callable(get_scheduler(scheduler=dask.local.synchronous_executor))
    assert callable(get_scheduler(scheduler=MyExecutor()))
    with dask.config.set(scheduler="threads"):
        assert get_scheduler() is dask.threaded.get
    assert get_scheduler() is None


def test_callable_scheduler():
    called = [False]

    def get(dsk, keys, *args, **kwargs):
        called[0] = True
        return dask.get(dsk, keys)

    assert delayed(lambda: 1)().compute(scheduler=get) == 1
    assert called[0]


@pytest.mark.flaky(reruns=10, reruns_delay=5)
@pytest.mark.slow
@pytest.mark.parametrize("scheduler", ["threads", "processes"])
def test_num_workers_config(scheduler):
    # Regression test for issue #4082

    f = delayed(pure=False)(time.sleep)
    # Be generous with the initial sleep times, as process have been observed
    # to take >0.5s to spin up
    num_workers = 3
    a = [f(1.0) for i in range(num_workers)]
    with dask.config.set(num_workers=num_workers, chunksize=1), Profiler() as prof:
        compute(*a, scheduler=scheduler)

    workers = {i.worker_id for i in prof.results}

    assert len(workers) == num_workers


def test_clone_key():
    for key, seed in [("x", 123), (("x", 1), 456), (("sum-1-2-3", h1, 1), 123)]:
        validate_key(clone_key(key, seed))
        assert clone_key(key, seed) != key
        assert clone_key(key, seed) == clone_key(key, seed)
        assert clone_key(key, seed) != clone_key(key, seed + 1)
        assert key_split(clone_key(key, seed)) == key_split(key)

    with pytest.raises(TypeError):
        clone_key(1, 2)


def test_compute_as_if_collection_low_level_task_graph():
    # See https://github.com/dask/dask/pull/7969
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    x = da.arange(10)

    # Boolean flag to ensure MyDaskArray.__dask_optimize__ is called
    optimized = False

    class MyDaskArray(da.Array):
        """Dask Array subclass with validation logic in __dask_optimize__"""

        @classmethod
        def __dask_optimize__(cls, dsk, keys, **kwargs):
            # Ensure `compute_as_if_collection` don't convert to a low-level task graph
            assert type(dsk) is HighLevelGraph
            nonlocal optimized
            optimized = True
            return super().__dask_optimize__(dsk, keys, **kwargs)

    result = compute_as_if_collection(
        MyDaskArray, x.__dask_graph__(), x.__dask_keys__()
    )[0]
    assert optimized
    da.utils.assert_eq(x, result)


# A function designed to be run in a subprocess with dask._compatibility.EMSCRIPTEN
# patched. This allows for checking for different default schedulers depending on the
# platform. One might prefer patching `sys.platform` for a more direct test, but that
# causes problems in other libraries.
def check_default_scheduler(module, collection, expected, emscripten):
    from contextlib import nullcontext
    from unittest import mock

    from dask.local import get_sync

    if emscripten:
        ctx = mock.patch("dask.base.named_schedulers", {"sync": get_sync})
    else:
        ctx = nullcontext()
    with ctx:
        import importlib

        if expected == "sync":
            from dask.local import get_sync as get
        elif expected == "threads":
            from dask.threaded import get
        elif expected == "processes":
            from dask.multiprocessing import get

        mod = importlib.import_module(module)

        assert getattr(mod, collection).__dask_scheduler__ == get


@pytest.mark.parametrize(
    "params",
    (
        "'dask.array', 'Array', 'sync', True",
        "'dask.array', 'Array', 'threads', False",
        "'dask.bag', 'Bag', 'sync', True",
        "'dask.bag', 'Bag', 'processes', False",
    ),
)
def test_emscripten_default_scheduler(params):
    pytest.importorskip("numpy")
    pytest.importorskip("dask.array")
    pytest.importorskip("pandas")
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                inspect.getsource(check_default_scheduler)
                + f"check_default_scheduler({params})\n"
            ),
        ]
    )
    proc.check_returncode()
