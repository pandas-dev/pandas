from __future__ import annotations

import contextlib
import copy
import pathlib
import re
import xml.etree.ElementTree

import pytest

from dask._task_spec import DataNode, List, Task, TaskRef

np = pytest.importorskip("numpy")

import itertools
import math
import operator
import os
import time
import warnings
from operator import add, sub
from threading import Lock

from packaging.version import Version
from tlz import concat, merge
from tlz.curried import identity

import dask
import dask.array as da
from dask.array.chunk import getitem
from dask.array.core import (
    Array,
    BlockView,
    blockdims_from_blockshape,
    broadcast_chunks,
    broadcast_shapes,
    broadcast_to,
    common_blockdim,
    concatenate,
    concatenate3,
    concatenate_axes,
    dotmany,
    from_array,
    from_delayed,
    from_func,
    getter,
    graph_from_arraylike,
    load_store_chunk,
    normalize_chunks,
    normalize_chunks_cached,
    optimize,
    stack,
    store,
)
from dask.array.numpy_compat import NUMPY_GE_200, NUMPY_GE_210
from dask.array.reshape import _not_implemented_message
from dask.array.utils import assert_eq, same_keys
from dask.base import collections_to_dsk, compute_as_if_collection, tokenize
from dask.blockwise import (
    _make_blockwise_graph,
    broadcast_dimensions,
    optimize_blockwise,
)
from dask.delayed import Delayed, delayed
from dask.highlevelgraph import HighLevelGraph, MaterializedLayer
from dask.layers import Blockwise
from dask.utils import SerializableLock, key_split, tmpdir, tmpfile
from dask.utils_test import dec, hlg_layer_topological, inc

if da._array_expr_enabled():
    pytest.skip("parametrize using unsupported functions", allow_module_level=True)


def skip_if_no_sparray():
    try:
        import scipy
    except ImportError:
        skip = True
    else:
        skip = Version(scipy.__version__) < Version("1.11")

    return pytest.mark.skipif(skip, reason="scipy<1.11 has no sparray")


@pytest.mark.parametrize("inline_array", [True, False])
def test_graph_from_arraylike(inline_array):
    d = 2
    chunk = (2, 3)
    shape = tuple(d * n for n in chunk)
    arr = np.ones(shape)

    dsk = graph_from_arraylike(
        arr, chunk, shape=shape, name="X", inline_array=inline_array
    )

    assert isinstance(dsk, HighLevelGraph)
    if inline_array:
        assert len(dsk.layers) == 1
        assert isinstance(hlg_layer_topological(dsk, 0), Blockwise)
    else:
        assert len(dsk.layers) == 2
        assert isinstance(hlg_layer_topological(dsk, 0), MaterializedLayer)
        assert isinstance(hlg_layer_topological(dsk, 1), Blockwise)
    dsk = dict(dsk)

    # Somewhat odd membership check to avoid numpy elemwise __in__ overload
    assert any(arr is v for v in dsk.values()) is not inline_array


def test_top():
    t = Task("inc", inc, TaskRef("A"))
    assert _make_blockwise_graph(
        t, "z", "ij", "x", "ij", numblocks={"x": (2, 2)}, keys=("A",)
    ) == {
        ("z", 0, 0): Task(("z", 0, 0), inc, TaskRef(("x", 0, 0))),
        ("z", 0, 1): Task(("z", 0, 1), inc, TaskRef(("x", 0, 1))),
        ("z", 1, 0): Task(("z", 1, 0), inc, TaskRef(("x", 1, 0))),
        ("z", 1, 1): Task(("z", 1, 1), inc, TaskRef(("x", 1, 1))),
    }

    keys = list(map(str, range(2)))
    t = Task("add", add, *(TaskRef(k) for k in keys))
    assert _make_blockwise_graph(
        t,
        "z",
        "ij",
        "x",
        "ij",
        "y",
        "ij",
        numblocks={"x": (2, 2), "y": (2, 2)},
        keys=keys,
    ) == {
        ("z", 0, 0): Task(("z", 0, 0), add, TaskRef(("x", 0, 0)), TaskRef(("y", 0, 0))),
        ("z", 0, 1): Task(("z", 0, 1), add, TaskRef(("x", 0, 1)), TaskRef(("y", 0, 1))),
        ("z", 1, 0): Task(("z", 1, 0), add, TaskRef(("x", 1, 0)), TaskRef(("y", 1, 0))),
        ("z", 1, 1): Task(("z", 1, 1), add, TaskRef(("x", 1, 1)), TaskRef(("y", 1, 1))),
    }
    keys = list(map(str, range(2)))
    t = Task(
        "dotmany",
        dotmany,
        TaskRef(keys[0]),
        TaskRef(keys[1]),
    )
    out = _make_blockwise_graph(
        t,
        "z",
        "ik",
        "x",
        "ij",
        "y",
        "jk",
        numblocks={"x": (2, 2), "y": (2, 2)},
        keys=keys,
    )
    expected = {
        ("z", 0, 0): Task(
            ("z", 0, 0),
            dotmany,
            List(TaskRef(("x", 0, 0)), TaskRef(("x", 0, 1))),
            List(TaskRef(("y", 0, 0)), TaskRef(("y", 1, 0))),
        ),
        ("z", 0, 1): Task(
            ("z", 0, 1),
            dotmany,
            List(TaskRef(("x", 0, 0)), TaskRef(("x", 0, 1))),
            List(TaskRef(("y", 0, 1)), TaskRef(("y", 1, 1))),
        ),
        ("z", 1, 0): Task(
            ("z", 1, 0),
            dotmany,
            List(TaskRef(("x", 1, 0)), TaskRef(("x", 1, 1))),
            List(TaskRef(("y", 0, 0)), TaskRef(("y", 1, 0))),
        ),
        ("z", 1, 1): Task(
            ("z", 1, 1),
            dotmany,
            List(TaskRef(("x", 1, 0)), TaskRef(("x", 1, 1))),
            List(TaskRef(("y", 0, 1)), TaskRef(("y", 1, 1))),
        ),
    }
    assert out[("z", 0, 0)] == expected[("z", 0, 0)]
    assert out == expected

    t = Task("identity", identity, TaskRef("0"))
    assert _make_blockwise_graph(
        t, "z", "", "x", "ij", numblocks={"x": (2, 2)}, keys=("0",)
    ) == {
        ("z",): Task(
            ("z",),
            identity,
            List(
                List(TaskRef(("x", 0, 0)), TaskRef(("x", 0, 1))),
                List(TaskRef(("x", 1, 0)), TaskRef(("x", 1, 1))),
            ),
        )
    }


def test_top_supports_broadcasting_rules():

    t = Task("add", add, TaskRef("A"), TaskRef("B"))
    assert _make_blockwise_graph(
        t,
        "z",
        "ij",
        "x",
        "ij",
        "y",
        "ij",
        numblocks={"x": (1, 2), "y": (2, 1)},
        keys=("A", "B"),
    ) == {
        ("z", 0, 0): Task(("z", 0, 0), add, TaskRef(("x", 0, 0)), TaskRef(("y", 0, 0))),
        ("z", 0, 1): Task(("z", 0, 1), add, TaskRef(("x", 0, 1)), TaskRef(("y", 0, 0))),
        ("z", 1, 0): Task(("z", 1, 0), add, TaskRef(("x", 0, 0)), TaskRef(("y", 1, 0))),
        ("z", 1, 1): Task(("z", 1, 1), add, TaskRef(("x", 0, 1)), TaskRef(("y", 1, 0))),
    }


def test_top_literals():
    t = Task("add", add, TaskRef("A"), TaskRef("B"))
    assert _make_blockwise_graph(
        t, "z", "ij", "x", "ij", 123, None, numblocks={"x": (2, 2)}, keys=("A", "B")
    ) == {
        ("z", 0, 0): Task(("z", 0, 0), add, TaskRef(("x", 0, 0)), DataNode(None, 123)),
        ("z", 0, 1): Task(("z", 0, 1), add, TaskRef(("x", 0, 1)), DataNode(None, 123)),
        ("z", 1, 0): Task(("z", 1, 0), add, TaskRef(("x", 1, 0)), DataNode(None, 123)),
        ("z", 1, 1): Task(("z", 1, 1), add, TaskRef(("x", 1, 1)), DataNode(None, 123)),
    }


def test_blockwise_literals():
    x = da.ones((10, 10), chunks=(5, 5))
    z = da.blockwise(add, "ij", x, "ij", 100, None, dtype=x.dtype)
    assert_eq(z, x + 100)

    z = da.blockwise(
        lambda x, y, z: x * y + z, "ij", 2, None, x, "ij", 100, None, dtype=x.dtype
    )
    assert_eq(z, 2 * x + 100)

    z = da.blockwise(getitem, "ij", x, "ij", slice(None), None, dtype=x.dtype)
    assert_eq(z, x)


def test_blockwise_1_in_shape_I():
    def test_f(a, b):
        assert 1 in b.shape

    p, k, N = 7, 2, 5
    arr = da.blockwise(
        test_f,
        "x",
        da.zeros((2 * p, 9, k * N), chunks=(p, 3, k)),
        "xzt",
        da.zeros((2 * p, 9, 1), chunks=(p, 3, -1)),
        "xzt",
        concatenate=True,
        dtype=float,
    )
    arr.compute()


def test_blockwise_1_in_shape_II():
    def test_f(a, b):
        assert 1 in b.shape

    p, k, N = 7, 2, 5
    da.blockwise(
        test_f,
        "x",
        da.zeros((2 * p, 9, k * N, 8), chunks=(p, 9, k, 4)),
        "xztu",
        da.zeros((2 * p, 9, 1, 8), chunks=(p, 9, -1, 4)),
        "xztu",
        concatenate=True,
        dtype=float,
    ).compute()


def test_blockwise_1_in_shape_III():
    def test_f(a, b):
        assert 1 in b.shape

    k, N = 2, 5
    da.blockwise(
        test_f,
        "x",
        da.zeros((k * N, 9, 8), chunks=(k, 3, 4)),
        "xtu",
        da.zeros((1, 9, 8), chunks=(-1, 3, 4)),
        "xtu",
        concatenate=True,
        dtype=float,
    ).compute()


def test_concatenate3_on_scalars():
    assert_eq(concatenate3([1, 2]), np.array([1, 2]))


def test_chunked_dot_product():
    x = np.arange(400).reshape((20, 20))
    o = np.ones((20, 20))

    getx = graph_from_arraylike(x, (5, 5), shape=(20, 20), name="x")
    geto = graph_from_arraylike(o, (5, 5), shape=(20, 20), name="o")

    keys = list(map(str, range(2)))
    t = Task(
        "dotmany",
        dotmany,
        TaskRef(keys[0]),
        TaskRef(keys[1]),
    )
    result = _make_blockwise_graph(
        t,
        "out",
        "ik",
        "x",
        "ij",
        "o",
        "jk",
        numblocks={"x": (4, 4), "o": (4, 4)},
        keys=keys,
    )

    dsk = merge(getx, geto, result)
    out = dask.get(dsk, [[("out", i, j) for j in range(4)] for i in range(4)])

    assert_eq(np.dot(x, o), concatenate3(out))


def test_chunked_transpose_plus_one():
    x = np.arange(400).reshape((20, 20))

    getx = graph_from_arraylike(x, (5, 5), shape=(20, 20), name="x")

    f = Task("f", lambda x: x.T + 1, TaskRef("x"))
    comp = _make_blockwise_graph(
        f, "out", "ij", "x", "ji", numblocks={"x": (4, 4)}, keys=("x",)
    )

    dsk = merge(getx, comp)
    out = dask.get(dsk, [[("out", i, j) for j in range(4)] for i in range(4)])

    assert_eq(concatenate3(out), x.T + 1)


def test_broadcast_dimensions_works_with_singleton_dimensions():
    argpairs = [("x", "i")]
    numblocks = {"x": ((1,),)}
    assert broadcast_dimensions(argpairs, numblocks) == {"i": (1,)}


def test_broadcast_dimensions():
    argpairs = [("x", "ij"), ("y", "ij")]
    d = {"x": ("Hello", 1), "y": (1, (2, 3))}
    assert broadcast_dimensions(argpairs, d) == {"i": "Hello", "j": (2, 3)}


def test_Array():
    arr = object()  # arraylike is unimportant since we never compute
    shape = (1000, 1000)
    chunks = (100, 100)
    name = "x"
    dsk = graph_from_arraylike(arr, chunks, shape, name)
    a = Array(dsk, name, chunks, shape=shape, dtype="f8")

    assert a.numblocks == (10, 10)

    assert a.__dask_keys__() == [[("x", i, j) for j in range(10)] for i in range(10)]

    assert a.chunks == ((100,) * 10, (100,) * 10)

    assert a.shape == shape

    assert len(a) == shape[0]

    with pytest.raises(ValueError):
        Array(dsk, name, chunks, shape=shape)
    with pytest.raises(TypeError):
        Array(dsk, name, chunks, shape=shape, dtype="f8", meta=np.empty(0, 0))


def test_uneven_chunks():
    a = Array({}, "x", chunks=(3, 3), shape=(10, 10), dtype="f8")
    assert a.chunks == ((3, 3, 3, 1), (3, 3, 3, 1))


def test_numblocks_suppoorts_singleton_block_dims():
    arr = object()  # arraylike is unimportant since we never compute
    shape = (100, 10)
    chunks = (10, 10)
    name = "x"
    dsk = graph_from_arraylike(arr, chunks, shape, name)
    a = Array(dsk, name, chunks, shape=shape, dtype="f8")

    assert set(concat(a.__dask_keys__())) == {("x", i, 0) for i in range(10)}


def test_keys():
    dsk = {("x", i, j): () for i in range(5) for j in range(6)}
    dx = Array(dsk, "x", chunks=(10, 10), shape=(50, 60), dtype="f8")
    assert dx.__dask_keys__() == [[(dx.name, i, j) for j in range(6)] for i in range(5)]
    # Cache works
    assert dx.__dask_keys__() is dx.__dask_keys__()
    # Test mutating names clears key cache
    dx.dask = {("y", i, j): () for i in range(5) for j in range(6)}
    dx._name = "y"
    new_keys = [[(dx.name, i, j) for j in range(6)] for i in range(5)]
    assert dx.__dask_keys__() == new_keys
    assert np.array_equal(dx._key_array, np.array(new_keys, dtype="object"))
    d = Array({}, "x", (), shape=(), dtype="f8")
    assert d.__dask_keys__() == [("x",)]


def test_Array_computation():
    a = Array({("x", 0, 0): np.eye(3)}, "x", shape=(3, 3), chunks=(3, 3), dtype="f8")
    assert_eq(np.array(a), np.eye(3))
    assert isinstance(a.compute(), np.ndarray)
    assert float(a[0, 0]) == 1


@pytest.mark.parametrize("asarray", [np.asarray, np.asanyarray, np.array])
def test_numpy_asarray_dtype(asarray):
    """Test dtype= parameter of np.*array()"""
    x = da.arange(10, dtype="i2")
    nx = asarray(x, dtype="i4")
    assert nx.dtype == np.int32


@pytest.mark.parametrize(
    "asarray",
    [
        np.array,
        pytest.param(
            np.asarray,
            marks=pytest.mark.skipif(not NUMPY_GE_200, reason="no copy kwarg"),
        ),
        pytest.param(
            np.asanyarray,
            marks=pytest.mark.skipif(not NUMPY_GE_210, reason="no copy kwarg"),
        ),
    ],
)
@pytest.mark.parametrize("chunks", [5, 10])
def test_numpy_asarray_copy_true(asarray, chunks):
    """Test np.*array(x, copy=True)"""
    x = da.asarray(np.arange(10), chunks=chunks)

    nx = asarray(x, copy=True)
    nx[0] = 42
    # Did not write back to the buffer held in the Dask graph
    assert x[0].compute() == 0


@pytest.mark.parametrize(
    "asarray",
    [
        pytest.param(
            np.array,
            marks=pytest.mark.skipif(
                not NUMPY_GE_200, reason="copy kwarg not forwarded to dask"
            ),
        ),
        pytest.param(
            np.asarray,
            marks=pytest.mark.skipif(not NUMPY_GE_200, reason="no copy kwarg"),
        ),
        pytest.param(
            np.asanyarray,
            marks=pytest.mark.skipif(not NUMPY_GE_210, reason="no copy kwarg"),
        ),
    ],
)
@pytest.mark.parametrize("chunks", [5, 10])
def test_numpy_asarray_copy_false(asarray, chunks):
    """Test that np.*array(x, copy=False) is forbidden"""
    x = da.asarray(np.arange(10), chunks=chunks)
    # Loudly crash if numpy demands for a view of dask's memory
    with pytest.warns(FutureWarning, match="Can't acquire a memory view"):
        y = asarray(x, copy=False)
    assert_eq(y, np.arange(10))


@pytest.mark.parametrize(
    "asarray",
    [
        pytest.param(
            np.array,
            marks=pytest.mark.skipif(not NUMPY_GE_200, reason="copy=None invalid"),
        ),
        pytest.param(
            np.asarray,
            marks=pytest.mark.skipif(not NUMPY_GE_200, reason="no copy kwarg"),
        ),
        pytest.param(
            np.asanyarray,
            marks=pytest.mark.skipif(not NUMPY_GE_210, reason="no copy kwarg"),
        ),
    ],
)
@pytest.mark.parametrize("chunks", [5, 10])
def test_numpy_asarray_copy_none(asarray, chunks):
    """Test np.*array(x, copy=None)"""
    x = da.asarray(np.arange(10), chunks=chunks)

    nx = asarray(x, copy=None)
    nx[0] = 42
    # Did not write back to the buffer held in the Dask graph
    assert x[0].compute() == 0


@pytest.mark.parametrize("asarray", [np.asarray, np.asanyarray, np.array])
@pytest.mark.parametrize("chunks", [5, 10])
def test_numpy_asarray_copy_default(asarray, chunks):
    """Test that np.*array() never returns an object that shares
    a buffer with the dask graph or a process-local Worker
    """
    x = da.asarray(np.arange(10), chunks=chunks)

    # array() defaults to copy=True.
    # asarray() and asanyarray() default to copy=None.
    # For Dask, it doesn't make a difference.
    nx = asarray(x)
    nx[0] = 42
    # Did not write back to the buffer held in the Dask graph
    assert x[0].compute() == 0


def test_array_interface_deprecated_kwargs():
    x = da.ones(10)
    with pytest.warns(FutureWarning, match="ignored"):
        x.__array__(something="foo")


@pytest.mark.parametrize("chunks", [5, 10])
def test_compute_copy(chunks):
    """Test that compute() never returns an object that shares
    a buffer with the dask graph or a process-local Worker
    """
    x = da.asarray(np.arange(10), chunks=chunks)

    nx = x.compute()
    nx[0] = 42
    # Did not write back to the buffer held in the Dask graph
    assert x[0].compute() == 0

    (nx,) = dask.compute(x)
    nx[0] = 42
    assert x[0].compute() == 0


def test_Array_numpy_gufunc_call__array_ufunc__01():
    x = da.random.default_rng().normal(size=(3, 10, 10), chunks=(2, 10, 10))
    nx = x.compute()
    ny = np.linalg._umath_linalg.inv(nx)
    y = np.linalg._umath_linalg.inv(x)
    assert_eq(ny, y)


def test_Array_numpy_gufunc_call__array_ufunc__02():
    x = da.random.default_rng().normal(size=(3, 10, 10), chunks=(2, 10, 10))
    nx = x.compute()
    nw, nv = np.linalg._umath_linalg.eig(nx)
    w, v = np.linalg._umath_linalg.eig(x)
    assert_eq(nw, w)
    assert_eq(nv, v)


def test_stack():
    a, b, c = (
        Array(
            graph_from_arraylike(object(), chunks=(2, 3), shape=(4, 6), name=name),
            name,
            chunks=(2, 3),
            dtype="f8",
            shape=(4, 6),
        )
        for name in "ABC"
    )

    s = stack([a, b, c], axis=0)

    colon = slice(None, None, None)

    assert s.shape == (3, 4, 6)
    assert s.chunks == ((1, 1, 1), (2, 2), (3, 3))
    assert s.chunksize == (1, 2, 3)
    assert s.dask[(s.name, 0, 1, 0)] == (getitem, ("A", 1, 0), (None, colon, colon))
    assert s.dask[(s.name, 2, 1, 0)] == (getitem, ("C", 1, 0), (None, colon, colon))
    assert same_keys(s, stack([a, b, c], axis=0))

    s2 = stack([a, b, c], axis=1)
    assert s2.shape == (4, 3, 6)
    assert s2.chunks == ((2, 2), (1, 1, 1), (3, 3))
    assert s2.chunksize == (2, 1, 3)
    assert s2.dask[(s2.name, 0, 1, 0)] == (getitem, ("B", 0, 0), (colon, None, colon))
    assert s2.dask[(s2.name, 1, 1, 0)] == (getitem, ("B", 1, 0), (colon, None, colon))
    assert same_keys(s2, stack([a, b, c], axis=1))

    s2 = stack([a, b, c], axis=2)
    assert s2.shape == (4, 6, 3)
    assert s2.chunks == ((2, 2), (3, 3), (1, 1, 1))
    assert s2.chunksize == (2, 3, 1)
    assert s2.dask[(s2.name, 0, 1, 0)] == (getitem, ("A", 0, 1), (colon, colon, None))
    assert s2.dask[(s2.name, 1, 1, 2)] == (getitem, ("C", 1, 1), (colon, colon, None))
    assert same_keys(s2, stack([a, b, c], axis=2))

    pytest.raises(ValueError, lambda: stack([]))
    pytest.raises(ValueError, lambda: stack([a, b, c], axis=3))

    assert set(b.dask.keys()).issubset(s2.dask.keys())

    assert stack([a, b, c], axis=-1).chunks == stack([a, b, c], axis=2).chunks


def test_stack_zero_size():
    x = np.empty((2, 0, 3))
    y = da.from_array(x, chunks=1)

    result_np = np.concatenate([x, x])
    result_da = da.concatenate([y, y])

    assert_eq(result_np, result_da)


def test_short_stack():
    x = np.array([1])
    d = da.from_array(x, chunks=(1,))
    s = da.stack([d])
    assert s.shape == (1, 1)
    chunks = compute_as_if_collection(Array, s.dask, s.__dask_keys__())
    assert chunks[0][0].shape == (1, 1)


def test_stack_scalars():
    d = da.arange(4, chunks=2)

    s = da.stack([d.mean(), d.sum()])

    assert s.compute().tolist() == [np.arange(4).mean(), np.arange(4).sum()]


def test_stack_promote_type():
    i = np.arange(10, dtype="i4")
    f = np.arange(10, dtype="f4")
    di = da.from_array(i, chunks=5)
    df = da.from_array(f, chunks=5)
    res = da.stack([di, df])
    assert_eq(res, np.stack([i, f]))


def test_stack_rechunk():
    rng = da.random.default_rng()
    x = rng.random(10, chunks=5)
    y = rng.random(10, chunks=4)

    z = da.stack([x, y], axis=0)
    assert z.shape == (2, 10)
    assert z.chunks == ((1, 1), (4, 1, 3, 2))

    assert_eq(z, np.stack([x.compute(), y.compute()], axis=0))


def test_stack_unknown_chunksizes():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    a_df = pd.DataFrame({"x": np.arange(12)})
    b_df = pd.DataFrame({"y": np.arange(12) * 10})

    a_ddf = dd.from_pandas(a_df, sort=False, npartitions=3)
    b_ddf = dd.from_pandas(b_df, sort=False, npartitions=3)

    a_x = a_ddf.values
    b_x = b_ddf.values

    assert np.isnan(a_x.shape[0])
    assert np.isnan(b_x.shape[0])

    with pytest.raises(ValueError) as exc_info:
        da.stack([a_x, b_x], axis=0)

    assert "shape" in str(exc_info.value)
    assert "nan" in str(exc_info.value)

    c_x = da.stack([a_x, b_x], axis=0, allow_unknown_chunksizes=True)

    assert_eq(c_x, np.stack([a_df.values, b_df.values], axis=0))

    with pytest.raises(ValueError) as exc_info:
        da.stack([a_x, b_x], axis=1)

    assert "shape" in str(exc_info.value)
    assert "nan" in str(exc_info.value)

    c_x = da.stack([a_x, b_x], axis=1, allow_unknown_chunksizes=True)

    assert_eq(c_x, np.stack([a_df.values, b_df.values], axis=1))

    m_df = pd.DataFrame({"m": np.arange(12) * 100})
    n_df = pd.DataFrame({"n": np.arange(12) * 1000})

    m_ddf = dd.from_pandas(m_df, sort=False, npartitions=3)
    n_ddf = dd.from_pandas(n_df, sort=False, npartitions=3)

    m_x = m_ddf.values
    n_x = n_ddf.values

    assert np.isnan(m_x.shape[0])
    assert np.isnan(n_x.shape[0])

    with pytest.raises(ValueError) as exc_info:
        da.stack([[a_x, b_x], [m_x, n_x]])

    assert "shape" in str(exc_info.value)
    assert "nan" in str(exc_info.value)

    c_x = da.stack([[a_x, b_x], [m_x, n_x]], allow_unknown_chunksizes=True)

    assert_eq(c_x, np.stack([[a_df.values, b_df.values], [m_df.values, n_df.values]]))


def test_concatenate():
    a, b, c = (
        Array(
            graph_from_arraylike(object(), chunks=(2, 3), shape=(4, 6), name=name),
            name,
            chunks=(2, 3),
            dtype="f8",
            shape=(4, 6),
        )
        for name in "ABC"
    )

    x = concatenate([a, b, c], axis=0)

    assert x.shape == (12, 6)
    assert x.chunks == ((2, 2, 2, 2, 2, 2), (3, 3))
    assert x.dask[(x.name, 0, 1)] == ("A", 0, 1)
    assert x.dask[(x.name, 5, 0)] == ("C", 1, 0)
    assert same_keys(x, concatenate([a, b, c], axis=0))

    y = concatenate([a, b, c], axis=1)

    assert y.shape == (4, 18)
    assert y.chunks == ((2, 2), (3, 3, 3, 3, 3, 3))
    assert y.dask[(y.name, 1, 0)] == ("A", 1, 0)
    assert y.dask[(y.name, 1, 5)] == ("C", 1, 1)
    assert same_keys(y, concatenate([a, b, c], axis=1))

    assert set(b.dask.keys()).issubset(y.dask.keys())

    z = concatenate([a], axis=0)

    assert z.shape == a.shape
    assert z.chunks == a.chunks
    assert z.dask == a.dask
    assert z is a

    assert (
        concatenate([a, b, c], axis=-1).chunks == concatenate([a, b, c], axis=1).chunks
    )

    pytest.raises(ValueError, lambda: concatenate([]))
    pytest.raises(ValueError, lambda: concatenate([a, b, c], axis=2))


@pytest.mark.parametrize(
    "dtypes", [((">f8", ">f8"), "float64"), (("<f4", "<f8"), "float64")]
)
def test_concatenate_types(dtypes):
    dts_in, dt_out = dtypes
    arrs = [np.zeros(4, dtype=dt) for dt in dts_in]
    darrs = [from_array(arr, chunks=(2,)) for arr in arrs]

    x = concatenate(darrs, axis=0)
    assert x.dtype == dt_out


def test_concatenate_unknown_axes():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    a_df = pd.DataFrame({"x": np.arange(12)})
    b_df = pd.DataFrame({"y": np.arange(12) * 10})

    a_ddf = dd.from_pandas(a_df, sort=False, npartitions=3)
    b_ddf = dd.from_pandas(b_df, sort=False, npartitions=3)

    a_x = a_ddf.values
    b_x = b_ddf.values

    assert np.isnan(a_x.shape[0])
    assert np.isnan(b_x.shape[0])

    da.concatenate([a_x, b_x], axis=0)  # works fine

    with pytest.raises(ValueError) as exc_info:
        da.concatenate([a_x, b_x], axis=1)  # unknown chunks

    assert "nan" in str(exc_info.value)
    assert "allow_unknown_chunksize" in str(exc_info.value)

    c_x = da.concatenate(
        [a_x, b_x], axis=1, allow_unknown_chunksizes=True
    )  # unknown chunks

    assert_eq(c_x, np.concatenate([a_df.values, b_df.values], axis=1))


def test_concatenate_flatten():
    x = np.array([1, 2])
    y = np.array([[3, 4], [5, 6]])

    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2, 1))

    assert_eq(np.concatenate([x, y], axis=None), da.concatenate([a, b], axis=None))


def test_concatenate_rechunk():
    rng = da.random.default_rng()
    x = rng.random((6, 6), chunks=(3, 3))
    y = rng.random((6, 6), chunks=(2, 2))

    z = da.concatenate([x, y], axis=0)
    assert z.shape == (12, 6)
    assert z.chunks == ((3, 3, 2, 2, 2), (2, 1, 1, 2))
    assert_eq(z, np.concatenate([x.compute(), y.compute()], axis=0))

    z = da.concatenate([x, y], axis=1)
    assert z.shape == (6, 12)
    assert z.chunks == ((2, 1, 1, 2), (3, 3, 2, 2, 2))
    assert_eq(z, np.concatenate([x.compute(), y.compute()], axis=1))


def test_concatenate_fixlen_strings():
    x = np.array(["a", "b", "c"])
    y = np.array(["aa", "bb", "cc"])

    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))

    assert_eq(np.concatenate([x, y]), da.concatenate([a, b]))


def test_concatenate_zero_size():
    x = np.random.default_rng().random(10)
    y = da.from_array(x, chunks=3)
    result_np = np.concatenate([x, x[:0]])
    result_da = da.concatenate([y, y[:0]])
    assert_eq(result_np, result_da)
    assert result_da is y

    # dtype of a size 0 arrays can affect the output dtype
    result_np = np.concatenate([np.zeros(0, dtype=float), np.zeros(1, dtype=int)])
    result_da = da.concatenate([da.zeros(0, dtype=float), da.zeros(1, dtype=int)])

    assert_eq(result_np, result_da)

    # All empty arrays case
    result_np = np.concatenate([np.zeros(0), np.zeros(0)])
    result_da = da.concatenate([da.zeros(0), da.zeros(0)])

    assert_eq(result_np, result_da)


def test_block_simple_row_wise():
    a1 = np.ones((2, 2))
    a2 = 2 * a1

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([a1, a2])
    result = da.block([d1, d2])

    assert_eq(expected, result)

    expected = np.block([a1, a2[:, :0]])
    result = da.block([d1, d2[:, :0]])

    assert result is d1
    assert_eq(expected, result)


def test_block_simple_column_wise():
    a1 = np.ones((2, 2))
    a2 = 2 * a1

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1], [a2]])
    result = da.block([[d1], [d2]])

    assert_eq(expected, result)


def test_block_with_1d_arrays_row_wise():
    # # # 1-D vectors are treated as row arrays
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([a1, a2])
    result = da.block([d1, d2])

    assert_eq(expected, result)

    expected = np.block([a1, a2[:0]])
    result = da.block([d1, d2[:0]])

    assert result is d1
    assert_eq(expected, result)


def test_block_with_1d_arrays_multiple_rows():
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1, a2], [a1, a2]])
    result = da.block([[d1, d2], [d1, d2]])

    assert_eq(expected, result)


def test_block_with_1d_arrays_column_wise():
    # # # 1-D vectors are treated as row arrays
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1], [a2]])
    result = da.block([[d1], [d2]])

    assert_eq(expected, result)


def test_block_mixed_1d_and_2d():
    a1 = np.ones((2, 2))
    a2 = np.array([2, 2])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[d1], [d2]])
    result = da.block([[a1], [a2]])

    assert_eq(expected, result)


def test_block_complicated():
    # a bit more complicated
    a1 = np.array([[1, 1, 1]])
    a2 = np.array([[2, 2, 2]])
    a3 = np.array([[3, 3, 3, 3, 3, 3]])
    a4 = np.array([4, 4, 4, 4, 4, 4])
    a5 = np.array(5)
    a6 = np.array([6, 6, 6, 6, 6])
    a7 = np.zeros((2, 6))

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)
    d3 = da.asarray(a3)
    d4 = da.asarray(a4)
    d5 = da.asarray(a5)
    d6 = da.asarray(a6)
    d7 = da.asarray(a7)

    expected = np.block([[a1, a2], [a3], [a4], [a5, a6], [a7]])
    result = da.block([[d1, d2], [d3], [d4], [d5, d6], [d7]])

    assert_eq(expected, result)


def test_block_nested():
    a1 = np.array([1, 1, 1])
    a2 = np.array([[2, 2, 2], [2, 2, 2], [2, 2, 2]])
    a3 = np.array([3, 3, 3])
    a4 = np.array([4, 4, 4])
    a5 = np.array(5)
    a6 = np.array([6, 6, 6, 6, 6])
    a7 = np.zeros((2, 6))

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)
    d3 = da.asarray(a3)
    d4 = da.asarray(a4)
    d5 = da.asarray(a5)
    d6 = da.asarray(a6)
    d7 = da.asarray(a7)

    expected = np.block([[np.block([[a1], [a3], [a4]]), a2], [a5, a6], [a7]])
    result = da.block([[da.block([[d1], [d3], [d4]]), d2], [d5, d6], [d7]])

    assert_eq(expected, result)


def test_block_3d():
    a000 = np.ones((2, 2, 2), int) * 1

    a100 = np.ones((3, 2, 2), int) * 2
    a010 = np.ones((2, 3, 2), int) * 3
    a001 = np.ones((2, 2, 3), int) * 4

    a011 = np.ones((2, 3, 3), int) * 5
    a101 = np.ones((3, 2, 3), int) * 6
    a110 = np.ones((3, 3, 2), int) * 7

    a111 = np.ones((3, 3, 3), int) * 8

    d000 = da.asarray(a000)

    d100 = da.asarray(a100)
    d010 = da.asarray(a010)
    d001 = da.asarray(a001)

    d011 = da.asarray(a011)
    d101 = da.asarray(a101)
    d110 = da.asarray(a110)

    d111 = da.asarray(a111)

    expected = np.block([[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]])
    result = da.block([[[d000, d001], [d010, d011]], [[d100, d101], [d110, d111]]])

    assert_eq(expected, result)

    expected = np.block(
        [
            [[a000, a001[:, :, :0]], [a010[:, :0, :], a011[:, :0, :0]]],
            [[a100[:0, :, :], a101[:0, :, :0]], [a110[:0, :0, :], a111[:0, :0, :0]]],
        ]
    )
    result = da.block(
        [
            [[d000, d001[:, :, :0]], [d010[:, :0, :], d011[:, :0, :0]]],
            [[d100[:0, :, :], d101[:0, :, :0]], [d110[:0, :0, :], d111[:0, :0, :0]]],
        ]
    )

    assert result is d000
    assert_eq(expected, result)


def test_block_with_mismatched_shape():
    a = np.array([0, 0])
    b = np.eye(2)

    for arrays in [[a, b], [b, a]]:
        with pytest.raises(ValueError):
            da.block(arrays)


def test_block_no_lists():
    assert_eq(da.block(1), np.block(1))
    assert_eq(da.block(np.eye(3)), np.block(np.eye(3)))


def test_block_invalid_nesting():
    for arrays in [
        [1, [2]],
        [1, []],
        [[1], 2],
        [[], 2],
        [[[1], [2]], [[3, 4]], [5]],  # missing brackets
    ]:
        with pytest.raises(ValueError) as e:
            da.block(arrays)
        e.match(r"depths are mismatched")


def test_block_empty_lists():
    for arrays in [[], [[]], [[1], []]]:
        with pytest.raises(ValueError) as e:
            da.block(arrays)
        e.match(r"empty")


def test_block_tuple():
    for arrays in [([1, 2], [3, 4]), [(1, 2), (3, 4)]]:
        with pytest.raises(TypeError) as e:
            da.block(arrays)
        e.match(r"tuple")


def test_broadcast_shapes():
    assert () == broadcast_shapes()
    assert (2, 5) == broadcast_shapes((2, 5))
    assert (0, 5) == broadcast_shapes((0, 1), (1, 5))
    assert np.allclose(
        (2, np.nan), broadcast_shapes((1, np.nan), (2, 1)), equal_nan=True
    )
    assert np.allclose(
        (2, np.nan), broadcast_shapes((2, 1), (1, np.nan)), equal_nan=True
    )
    assert (3, 4, 5) == broadcast_shapes((3, 4, 5), (4, 1), ())
    assert (3, 4) == broadcast_shapes((3, 1), (1, 4), (4,))
    assert (5, 6, 7, 3, 4) == broadcast_shapes((3, 1), (), (5, 6, 7, 1, 4))

    assert all(
        isinstance(i, int) for i in itertools.chain(*broadcast_shapes([(1, 2), (3, 1)]))
    )

    pytest.raises(ValueError, lambda: broadcast_shapes((3,), (3, 4)))
    pytest.raises(ValueError, lambda: broadcast_shapes((2, 3), (2, 3, 1)))
    pytest.raises(ValueError, lambda: broadcast_shapes((2, 3), (1, np.nan)))


def test_elemwise_on_scalars():
    nx = np.arange(10, dtype=np.int64)
    ny = np.arange(10, dtype=np.int32)
    nz = nx.sum() * ny

    dx = from_array(nx, chunks=(5,))
    dy = from_array(ny, chunks=(5,))
    dz = dx.sum() * dy

    if NUMPY_GE_200:
        assert_eq(dz, nz)
    else:
        # Dask 0-d arrays do not behave like numpy scalars for type promotion
        assert_eq(dz, nz, check_dtype=False)
        assert nz.dtype == np.int32
        assert dz.dtype == np.int64
        assert dz.compute().dtype == np.int64


def test_elemwise_with_ndarrays():
    x = np.arange(3)
    y = np.arange(12).reshape(4, 3)
    a = from_array(x, chunks=(3,))
    b = from_array(y, chunks=(2, 3))

    assert_eq(x + a, 2 * x)
    assert_eq(a + x, 2 * x)

    assert_eq(x + b, x + y)
    assert_eq(b + x, x + y)
    assert_eq(a + y, x + y)
    assert_eq(y + a, x + y)
    # Error on shape mismatch
    pytest.raises(ValueError, lambda: a + y.T)
    pytest.raises(ValueError, lambda: a + np.arange(2))


def test_elemwise_differently_chunked():
    x = np.arange(3)
    y = np.arange(12).reshape(4, 3)
    a = from_array(x, chunks=(3,))
    b = from_array(y, chunks=(2, 2))

    assert_eq(a + b, x + y)
    assert_eq(b + a, x + y)


@pytest.mark.filterwarnings("ignore:overflow encountered in cast")  # numpy >=2.0
def test_elemwise_dtype():
    values = [
        da.from_array(np.ones(5, np.float32), chunks=3),
        da.from_array(np.ones(5, np.int16), chunks=3),
        da.from_array(np.ones(5, np.int64), chunks=3),
        da.from_array(np.ones((), np.float64), chunks=()) * 1e200,
        np.ones(5, np.float32),
        1,
        1.0,
        1e200,
        np.int64(1),
        np.ones((), np.int64),
    ]
    for x in values:
        for y in values:
            assert da.maximum(x, y).dtype == da.result_type(x, y)


def test_operators():
    x = np.arange(10)
    y = np.arange(10).reshape((10, 1))
    a = from_array(x, chunks=(5,))
    b = from_array(y, chunks=(5, 1))

    c = a + 1
    assert_eq(c, x + 1)

    c = a + b
    assert_eq(c, x + x.reshape((10, 1)))

    expr = (3 / a * b) ** 2 > 5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # divide by zero
        assert_eq(expr, (3 / x * y) ** 2 > 5)

    c = da.exp(a)
    assert_eq(c, np.exp(x))

    assert_eq(abs(-a), a)
    assert_eq(a, +x)


def test_binary_operator_delegation(mocker):
    # Binary operator delegation in Dask should follow Numpy's:
    # https://numpy.org/neps/nep-0013-ufunc-overrides.html#behavior-in-combination-with-python-s-binary-operations

    x = from_array([1, 2, 3])

    # Mock various types of `other` objects
    ufunc_none = mocker.Mock()
    ufunc_none.__array_ufunc__ = None
    ufunc_none.__radd__ = mocker.Mock()

    ufunc_high_priority = mocker.Mock()
    ufunc_high_priority.__array_priority__ = x.__array_priority__ + 1
    ufunc_high_priority.__radd__ = mocker.Mock()

    ufunc_low_priority = mocker.Mock()
    ufunc_low_priority.__array_priority__ = x.__array_priority__ - 1
    ufunc_low_priority.__radd__ = mocker.Mock()

    ufunc_no_priority = mocker.Mock()
    ufunc_no_priority.__radd__ = mocker.Mock()

    # If `other.__array_ufunc__ is None`, delegates back to Python
    # and therefore call reflected operator on `other`
    x + ufunc_none
    ufunc_none.__radd__.assert_called_once()

    # If the `__array_ufunc__` attribute is absent on other and
    # `other.__array_priority__ > self.__array_priority__`, also delegates back
    # to Python and therefore call reflected operator on `other`
    x + ufunc_high_priority
    ufunc_high_priority.__radd__.assert_called_once()

    # If `other.__array_priority__ <= self.__array_priority__`, does not
    # delegate (here it raises an error)
    with pytest.raises(TypeError):
        x + ufunc_low_priority
    ufunc_low_priority.__radd__.assert_not_called()

    # If `other.__array_priority__` is absent, does not delegate (raises)
    with pytest.raises(TypeError):
        x + ufunc_no_priority
    ufunc_no_priority.__radd__.assert_not_called()


@pytest.mark.filterwarnings("ignore:overflow encountered in cast")  # numpy >=2.0
def test_operator_dtype_promotion():
    x = np.arange(10, dtype=np.float32)
    y = np.array([1])
    a = from_array(x, chunks=(5,))

    assert_eq(x + 1, a + 1)  # still float32
    assert_eq(x + 1e50, a + 1e50)  # now float64
    assert_eq(x + y, a + y)  # also float64


def test_field_access():
    x = np.array([(1, 1.0), (2, 2.0)], dtype=[("a", "i4"), ("b", "f4")])
    y = from_array(x, chunks=(1,))
    assert_eq(y["a"], x["a"])
    assert_eq(y[["b", "a"]], x[["b", "a"]])
    assert same_keys(y[["b", "a"]], y[["b", "a"]])


def test_field_access_with_shape():
    dtype = [("col1", ("f4", (3, 2))), ("col2", ("f4", 3))]
    data = np.ones((100, 50), dtype=dtype)
    x = da.from_array(data, 10)
    assert_eq(x["col1"], data["col1"])
    assert_eq(x[["col1"]], data[["col1"]])
    assert_eq(x["col2"], data["col2"])
    assert_eq(x[["col1", "col2"]], data[["col1", "col2"]])


def test_matmul():
    rng = np.random.default_rng()
    x = rng.random((5, 5))
    y = rng.random((5, 2))
    a = from_array(x, chunks=(1, 5))
    b = from_array(y, chunks=(5, 1))
    assert_eq(operator.matmul(a, b), a.dot(b))
    assert_eq(operator.matmul(a, b), operator.matmul(x, y))
    assert_eq(operator.matmul(a, y), operator.matmul(x, b))
    list_vec = list(range(1, 6))
    assert_eq(operator.matmul(list_vec, b), operator.matmul(list_vec, y))
    assert_eq(operator.matmul(x, list_vec), operator.matmul(a, list_vec))
    z = rng.random((5, 5, 5))
    c = from_array(z, chunks=(1, 5, 1))
    assert_eq(operator.matmul(a, z), operator.matmul(x, c))
    assert_eq(operator.matmul(z, a), operator.matmul(c, x))


def test_matmul_array_ufunc():
    # regression test for https://github.com/dask/dask/issues/4353
    rng = np.random.default_rng()
    x = rng.random((5, 5))
    y = rng.random((5, 2))
    a = from_array(x, chunks=(1, 5))
    b = from_array(y, chunks=(5, 1))
    result = b.__array_ufunc__(np.matmul, "__call__", a, b)
    assert_eq(result, x.dot(y))


def test_T():
    x = np.arange(400).reshape((20, 20))
    a = from_array(x, chunks=(5, 5))

    assert_eq(x.T, a.T)


def test_broadcast_to():
    x = np.random.default_rng().integers(10, size=(5, 1, 6))
    a = from_array(x, chunks=(3, 1, 3))

    for shape in [a.shape, (5, 0, 6), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        xb = np.broadcast_to(x, shape)
        ab = broadcast_to(a, shape)

        assert_eq(xb, ab)

        if a.shape == ab.shape:
            assert a is ab

    pytest.raises(ValueError, lambda: broadcast_to(a, (2, 1, 6)))
    pytest.raises(ValueError, lambda: broadcast_to(a, (3,)))


def test_broadcast_to_array():
    x = np.random.default_rng().integers(10, size=(5, 1, 6))

    for shape in [(5, 0, 6), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        a = np.broadcast_to(x, shape)
        d = broadcast_to(x, shape)

        assert_eq(a, d)


def test_broadcast_to_scalar():
    x = 5

    for shape in [tuple(), (0,), (2, 3), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        a = np.broadcast_to(x, shape)
        d = broadcast_to(x, shape)

        assert_eq(a, d)


def test_broadcast_to_chunks():
    x = np.random.default_rng().integers(10, size=(5, 1, 6))
    a = from_array(x, chunks=(3, 1, 3))

    for shape, chunks, expected_chunks in [
        ((5, 3, 6), (3, -1, 3), ((3, 2), (3,), (3, 3))),
        ((5, 3, 6), (3, 1, 3), ((3, 2), (1, 1, 1), (3, 3))),
        ((2, 5, 3, 6), (1, 3, 1, 3), ((1, 1), (3, 2), (1, 1, 1), (3, 3))),
    ]:
        xb = np.broadcast_to(x, shape)
        ab = broadcast_to(a, shape, chunks=chunks)
        assert_eq(xb, ab)
        assert ab.chunks == expected_chunks

    with pytest.raises(ValueError):
        broadcast_to(a, a.shape, chunks=((2, 3), (1,), (3, 3)))
    with pytest.raises(ValueError):
        broadcast_to(a, a.shape, chunks=((3, 2), (3,), (3, 3)))
    with pytest.raises(ValueError):
        broadcast_to(a, (5, 2, 6), chunks=((3, 2), (3,), (3, 3)))


def test_broadcast_arrays():
    assert np.broadcast_arrays() == da.broadcast_arrays()

    a = np.arange(4)
    d_a = da.from_array(a, chunks=tuple(s // 2 for s in a.shape))

    a_0 = np.arange(4)[None, :]
    a_1 = np.arange(4)[:, None]

    d_a_0 = d_a[None, :]
    d_a_1 = d_a[:, None]

    a_r = np.broadcast_arrays(a_0, a_1)
    d_r = da.broadcast_arrays(d_a_0, d_a_1)

    assert isinstance(d_r, (list, tuple))
    assert len(a_r) == len(d_r)

    for e_a_r, e_d_r in zip(a_r, d_r):
        assert_eq(e_a_r, e_d_r)


def test_broadcast_arrays_uneven_chunks():
    x = da.ones(30, chunks=(3,))
    y = da.ones(30, chunks=(5,))
    z = np.broadcast_arrays(x, y)

    assert_eq(z, z)

    x = da.ones((1, 30), chunks=(1, 3))
    y = da.ones(30, chunks=(5,))
    z = np.broadcast_arrays(x, y)

    assert_eq(z, z)


@pytest.mark.parametrize(
    "u_shape, v_shape",
    [
        [tuple(), (2, 3)],
        [(1,), (2, 3)],
        [(1, 1), (2, 3)],
        [(0, 3), (1, 3)],
        [(2, 0), (2, 1)],
        [(1, 0), (2, 1)],
        [(0, 1), (1, 3)],
    ],
)
def test_broadcast_operator(u_shape, v_shape):
    rng = np.random.default_rng()
    u = rng.random(u_shape)
    v = rng.random(v_shape)

    d_u = from_array(u, chunks=1)
    d_v = from_array(v, chunks=1)

    w = u * v
    d_w = d_u * d_v

    assert_eq(w, d_w)


@pytest.mark.parametrize(
    "original_shape,new_shape,chunks",
    [
        ((10,), (10,), (3, 3, 4)),
        ((10,), (10, 1, 1), 5),
        ((10,), (1, 10), 5),
        ((24,), (2, 3, 4), 12),
        ((1, 24), (2, 3, 4), 12),
        ((2, 3, 4), (24,), (1, 3, 4)),
        ((2, 3, 4), (24,), 4),
        ((2, 3, 4), (24, 1), 4),
        ((2, 3, 4), (1, 24), 4),
        ((4, 4, 1), (4, 4), 2),
        ((4, 4), (4, 4, 1), 2),
        ((1, 4, 4), (4, 4), 2),
        ((1, 4, 4), (4, 4, 1), 2),
        ((1, 4, 4), (1, 1, 4, 4), 2),
        ((4, 4), (1, 4, 4, 1), 2),
        ((4, 4), (1, 4, 4), 2),
        ((2, 3), (2, 3), (1, 2)),
        ((2, 3), (3, 2), 3),
        ((4, 2, 3), (4, 6), 4),
        ((3, 4, 5, 6), (3, 4, 5, 6), (2, 3, 4, 5)),
        ((), (1,), 1),
        ((1,), (), 1),
        ((24,), (3, 8), 24),
        ((24,), (4, 6), 6),
        ((24,), (4, 3, 2), 6),
        ((24,), (4, 6, 1), 6),
        ((24,), (4, 6), (6, 12, 6)),
        ((64, 4), (8, 8, 4), (16, 2)),
        ((4, 64), (4, 8, 4, 2), (2, 16)),
        ((4, 8, 4, 2), (2, 1, 2, 32, 2), (2, 4, 2, 2)),
        ((4, 1, 4), (4, 4), (2, 1, 2)),
        ((0, 10), (0, 5, 2), (5, 5)),
        ((5, 0, 2), (0, 10), (5, 2, 2)),
        ((0,), (2, 0, 2), (4,)),
        ((2, 0, 2), (0,), (4, 4, 4)),
    ],
)
def test_reshape(original_shape, new_shape, chunks):
    x = np.random.default_rng().integers(10, size=original_shape)
    a = from_array(x, chunks=chunks)

    xr = x.reshape(new_shape)
    ar = a.reshape(new_shape)

    if a.shape == new_shape:
        assert a is ar

    assert_eq(xr, ar)


def test_reshape_exceptions():
    x = np.random.default_rng().integers(10, size=(5,))
    a = from_array(x, chunks=(2,))
    with pytest.raises(ValueError):
        da.reshape(a, (100,))


def test_reshape_splat():
    x = da.ones((5, 5), chunks=(2, 2))
    assert_eq(x.reshape((25,)), x.reshape(25))


def test_reshape_not_implemented_error():
    a = da.ones((4, 5, 6), chunks=(2, 2, 3))
    for new_shape in [(2, 10, 6), (5, 4, 6), (6, 5, 4)]:
        with pytest.raises(
            NotImplementedError, match=re.escape(_not_implemented_message)
        ):
            a.reshape(new_shape)


def test_reshape_unknown_dimensions():
    for original_shape in [(24,), (2, 12), (2, 3, 4)]:
        for new_shape in [(-1,), (2, -1), (-1, 3, 4)]:
            x = np.random.default_rng().integers(10, size=original_shape)
            a = from_array(x, 24)
            assert_eq(x.reshape(new_shape), a.reshape(new_shape))

    pytest.raises(ValueError, lambda: da.reshape(a, (-1, -1)))


def test_full():
    d = da.full((3, 4), 2, chunks=((2, 1), (2, 2)))
    assert d.chunks == ((2, 1), (2, 2))
    assert_eq(d, np.full((3, 4), 2))


def test_map_blocks():
    x = np.arange(400).reshape((20, 20))
    d = from_array(x, chunks=(7, 7))

    e = d.map_blocks(inc, dtype=d.dtype)

    assert d.chunks == e.chunks
    assert_eq(e, x + 1)

    e = d.map_blocks(inc, name="increment")
    assert e.name.startswith("increment-")

    assert d.map_blocks(inc, name="foo").name != d.map_blocks(dec, name="foo").name

    d = from_array(x, chunks=(10, 10))
    e = d.map_blocks(lambda x: x[::2, ::2], chunks=(5, 5), dtype=d.dtype)

    assert e.chunks == ((5, 5), (5, 5))
    assert_eq(e, x[::2, ::2])

    d = from_array(x, chunks=(8, 8))
    e = d.map_blocks(
        lambda x: x[::2, ::2], chunks=((4, 4, 2), (4, 4, 2)), dtype=d.dtype
    )

    assert_eq(e, x[::2, ::2])


def test_map_blocks2():
    x = np.arange(10, dtype="i8")
    d = from_array(x, chunks=(2,))

    def func(block, block_id=None, c=0):
        return np.ones_like(block) * sum(block_id) + c

    out = d.map_blocks(func, dtype="i8")
    expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4], dtype="i8")

    assert_eq(out, expected)
    assert same_keys(d.map_blocks(func, dtype="i8"), out)

    out = d.map_blocks(func, dtype="i8", c=1)
    expected = expected + 1

    assert_eq(out, expected)
    assert same_keys(d.map_blocks(func, dtype="i8", c=1), out)


def test_map_blocks_block_info():
    x = da.arange(50, chunks=10)

    def func(a, b, c, block_info=None):
        for idx in [0, 2, None]:  # positions in args
            assert block_info[idx]["shape"] == (50,)
            assert block_info[idx]["num-chunks"] == (5,)
            start, stop = block_info[idx]["array-location"][0]
            assert stop - start == 10
            assert 0 <= start <= 40
            assert 10 <= stop <= 50

            assert 0 <= block_info[idx]["chunk-location"][0] <= 4
        assert block_info[None]["chunk-shape"] == (10,)
        assert block_info[None]["dtype"] == x.dtype

        return a + b + c

    z = da.map_blocks(func, x, 100, x + 1, dtype=x.dtype)
    assert_eq(z, x + x + 1 + 100)


def test_map_blocks_block_info_with_new_axis():
    # https://github.com/dask/dask/issues/4298
    values = da.from_array(np.array(["a", "a", "b", "c"]), 2)

    def func(x, block_info=None):
        assert block_info.keys() == {0, None}
        assert block_info[0]["shape"] == (4,)
        assert block_info[0]["num-chunks"] == (2,)
        assert block_info[None]["shape"] == (4, 3)
        assert block_info[None]["num-chunks"] == (2, 1)
        assert block_info[None]["chunk-shape"] == (2, 3)
        assert block_info[None]["dtype"] == np.dtype("f8")

        assert block_info[0]["chunk-location"] in {(0,), (1,)}

        if block_info[0]["chunk-location"] == (0,):
            assert block_info[0]["array-location"] == [(0, 2)]
            assert block_info[None]["chunk-location"] == (0, 0)
            assert block_info[None]["array-location"] == [(0, 2), (0, 3)]
        elif block_info[0]["chunk-location"] == (1,):
            assert block_info[0]["array-location"] == [(2, 4)]
            assert block_info[None]["chunk-location"] == (1, 0)
            assert block_info[None]["array-location"] == [(2, 4), (0, 3)]

        return np.ones((len(x), 3))

    z = values.map_blocks(func, chunks=((2, 2), 3), new_axis=1, dtype="f8")
    assert_eq(z, np.ones((4, 3), dtype="f8"))


def test_map_blocks_block_info_with_drop_axis():
    # https://github.com/dask/dask/issues/4584
    values = da.from_array(
        np.array(
            [[1, 2, 4], [8, 16, 32], [64, 128, 256], [1024, 2048, 4096]], dtype="u4"
        ),
        (2, 1),
    )

    def func(x, block_info=None):
        assert block_info.keys() == {0, None}
        assert block_info[0]["shape"] == (4, 3)
        # drop_axis concatenates along the dropped dimension, hence not (2, 3)
        assert block_info[0]["num-chunks"] == (2, 1)
        assert block_info[None]["shape"] == (4,)
        assert block_info[None]["num-chunks"] == (2,)
        assert block_info[None]["chunk-shape"] == (2,)
        assert block_info[None]["dtype"] == np.dtype("u4")

        assert block_info[0]["chunk-location"] in {(0, 0), (1, 0)}

        if block_info[0]["chunk-location"] == (0, 0):
            assert block_info[0]["array-location"] == [(0, 2), (0, 3)]
            assert block_info[None]["chunk-location"] == (0,)
            assert block_info[None]["array-location"] == [(0, 2)]
        elif block_info[0]["chunk-location"] == (1, 0):
            assert block_info[0]["array-location"] == [(2, 4), (0, 3)]
            assert block_info[None]["chunk-location"] == (1,)
            assert block_info[None]["array-location"] == [(2, 4)]

        return np.sum(x, axis=1, dtype="u4")

    z = values.map_blocks(func, drop_axis=1, dtype="u4")
    assert_eq(z, np.array([7, 56, 448, 7168], dtype="u4"))


def test_map_blocks_block_info_with_broadcast():
    expected0 = [
        {
            "shape": (3, 4),
            "num-chunks": (1, 2),
            "array-location": [(0, 3), (0, 2)],
            "chunk-location": (0, 0),
        },
        {
            "shape": (3, 4),
            "num-chunks": (1, 2),
            "array-location": [(0, 3), (2, 4)],
            "chunk-location": (0, 1),
        },
    ]
    expected1 = [
        {
            "shape": (6, 2),
            "num-chunks": (2, 1),
            "array-location": [(0, 3), (0, 2)],
            "chunk-location": (0, 0),
        },
        {
            "shape": (6, 2),
            "num-chunks": (2, 1),
            "array-location": [(3, 6), (0, 2)],
            "chunk-location": (1, 0),
        },
    ]
    expected2 = [
        {
            "shape": (4,),
            "num-chunks": (2,),
            "array-location": [(0, 2)],
            "chunk-location": (0,),
        },
        {
            "shape": (4,),
            "num-chunks": (2,),
            "array-location": [(2, 4)],
            "chunk-location": (1,),
        },
    ]
    expected = [
        {
            0: expected0[0],
            1: expected1[0],
            2: expected2[0],
            None: {
                "shape": (6, 4),
                "num-chunks": (2, 2),
                "dtype": np.float64,
                "chunk-shape": (3, 2),
                "array-location": [(0, 3), (0, 2)],
                "chunk-location": (0, 0),
            },
        },
        {
            0: expected0[1],
            1: expected1[0],
            2: expected2[1],
            None: {
                "shape": (6, 4),
                "num-chunks": (2, 2),
                "dtype": np.float64,
                "chunk-shape": (3, 2),
                "array-location": [(0, 3), (2, 4)],
                "chunk-location": (0, 1),
            },
        },
        {
            0: expected0[0],
            1: expected1[1],
            2: expected2[0],
            None: {
                "shape": (6, 4),
                "num-chunks": (2, 2),
                "dtype": np.float64,
                "chunk-shape": (3, 2),
                "array-location": [(3, 6), (0, 2)],
                "chunk-location": (1, 0),
            },
        },
        {
            0: expected0[1],
            1: expected1[1],
            2: expected2[1],
            None: {
                "shape": (6, 4),
                "num-chunks": (2, 2),
                "dtype": np.float64,
                "chunk-shape": (3, 2),
                "array-location": [(3, 6), (2, 4)],
                "chunk-location": (1, 1),
            },
        },
    ]

    def func(x, y, z, block_info=None):
        for info in expected:
            if block_info[None]["chunk-location"] == info[None]["chunk-location"]:
                assert block_info == info
                break
        else:
            assert False
        return x + y + z

    a = da.ones((3, 4), chunks=(3, 2))
    b = da.ones((6, 2), chunks=(3, 2))
    c = da.ones((4,), chunks=(2,))
    d = da.map_blocks(func, a, b, c, chunks=((3, 3), (2, 2)), dtype=a.dtype)
    assert d.chunks == ((3, 3), (2, 2))
    assert_eq(d, 3 * np.ones((6, 4)))


def test_map_blocks_with_constants():
    d = da.arange(10, chunks=3)
    e = d.map_blocks(add, 100, dtype=d.dtype)

    assert_eq(e, np.arange(10) + 100)

    assert_eq(da.map_blocks(sub, d, 10, dtype=d.dtype), np.arange(10) - 10)
    assert_eq(da.map_blocks(sub, 10, d, dtype=d.dtype), 10 - np.arange(10))


def test_map_blocks_with_kwargs():
    d = da.arange(10, chunks=5)

    result = d.map_blocks(np.max, axis=0, keepdims=True, dtype=d.dtype, chunks=(1,))

    assert_eq(result, np.array([4, 9]))


def test_map_blocks_infer_chunks_broadcast():
    dx = da.from_array([[1, 2, 3, 4]], chunks=((1,), (2, 2)))
    dy = da.from_array([[10, 20], [30, 40]], chunks=((1, 1), (2,)))
    result = da.map_blocks(lambda x, y: x + y, dx, dy)
    assert result.chunks == ((1, 1), (2, 2))
    assert_eq(result, np.array([[11, 22, 13, 24], [31, 42, 33, 44]]))


def test_map_blocks_with_chunks():
    dx = da.ones((5, 3), chunks=(2, 2))
    dy = da.ones((5, 3), chunks=(2, 2))
    dz = da.map_blocks(np.add, dx, dy, chunks=dx.chunks)
    assert_eq(dz, np.ones((5, 3)) * 2)


def test_map_blocks_dtype_inference():
    x = np.arange(50).reshape((5, 10))
    y = np.arange(10)
    dx = da.from_array(x, chunks=5)
    dy = da.from_array(y, chunks=5)

    def foo(x, *args, **kwargs):
        cast = kwargs.pop("cast", "i8")
        return (x + sum(args)).astype(cast)

    assert_eq(dx.map_blocks(foo, dy, 1), foo(dx, dy, 1))
    assert_eq(dx.map_blocks(foo, dy, 1, cast="f8"), foo(dx, dy, 1, cast="f8"))
    assert_eq(
        dx.map_blocks(foo, dy, 1, cast="f8", dtype="f8"),
        foo(dx, dy, 1, cast="f8", dtype="f8"),
    )

    def foo(x):
        raise RuntimeError("Woops")

    with pytest.raises(ValueError) as e:
        dx.map_blocks(foo)
    msg = str(e.value)
    assert "dtype" in msg


def test_map_blocks_infer_newaxis():
    x = da.ones((5, 3), chunks=(2, 2))
    y = da.map_blocks(lambda x: x[None], x, chunks=((1,), (2, 2, 1), (2, 1)))
    assert_eq(y, da.ones((1, 5, 3)))


def test_map_blocks_no_array_args():
    def func(dtype, block_info=None):
        loc = block_info[None]["array-location"]
        return np.arange(loc[0][0], loc[0][1], dtype=dtype)

    x = da.map_blocks(func, np.float32, chunks=((5, 3),), dtype=np.float32)
    assert x.chunks == ((5, 3),)
    assert_eq(x, np.arange(8, dtype=np.float32))


def test_map_blocks_unique_name_enforce_dim():
    def func(some_3d, block_info=None):
        return some_3d

    input_arr = da.zeros((3, 4, 5), chunks=((3,), (4,), (5,)), dtype=np.float32)
    x = da.map_blocks(func, input_arr, enforce_ndim=True, dtype=np.float32)
    y = da.map_blocks(func, input_arr, enforce_ndim=False, dtype=np.float32)
    assert x._name != y._name


def test_map_blocks_unique_name_chunks_dtype():
    def func(block_info=None):
        loc = block_info[None]["array-location"]
        dtype = block_info[None]["dtype"]
        return np.arange(loc[0][0], loc[0][1], dtype=dtype)

    x = da.map_blocks(func, chunks=((5, 3),), dtype=np.float32)
    assert x.chunks == ((5, 3),)
    assert_eq(x, np.arange(8, dtype=np.float32))

    y = da.map_blocks(func, chunks=((2, 2, 1, 3),), dtype=np.float32)
    assert y.chunks == ((2, 2, 1, 3),)
    assert_eq(y, np.arange(8, dtype=np.float32))
    assert x.name != y.name

    z = da.map_blocks(func, chunks=((5, 3),), dtype=np.float64)
    assert z.chunks == ((5, 3),)
    assert_eq(z, np.arange(8, dtype=np.float64))
    assert x.name != z.name
    assert y.name != z.name


def test_map_blocks_unique_name_drop_axis():
    def func(some_3d, block_info=None):
        if not block_info:
            return some_3d
        dtype = block_info[None]["dtype"]
        return np.zeros(block_info[None]["shape"], dtype=dtype)

    input_arr = da.zeros((3, 4, 5), chunks=((3,), (4,), (5,)), dtype=np.float32)
    x = da.map_blocks(func, input_arr, drop_axis=[0], dtype=np.float32)
    assert x.chunks == ((4,), (5,))
    assert_eq(x, np.zeros((4, 5), dtype=np.float32))

    y = da.map_blocks(func, input_arr, drop_axis=[2], dtype=np.float32)
    assert y.chunks == ((3,), (4,))
    assert_eq(y, np.zeros((3, 4), dtype=np.float32))
    assert x.name != y.name


def test_map_blocks_unique_name_new_axis():
    def func(some_2d, block_info=None):
        if not block_info:
            return some_2d
        dtype = block_info[None]["dtype"]
        return np.zeros(block_info[None]["shape"], dtype=dtype)

    input_arr = da.zeros((3, 4), chunks=((3,), (4,)), dtype=np.float32)
    x = da.map_blocks(func, input_arr, new_axis=[0], dtype=np.float32)
    assert x.chunks == ((1,), (3,), (4,))
    assert_eq(x, np.zeros((1, 3, 4), dtype=np.float32))

    y = da.map_blocks(func, input_arr, new_axis=[2], dtype=np.float32)
    assert y.chunks == ((3,), (4,), (1,))
    assert_eq(y, np.zeros((3, 4, 1), dtype=np.float32))
    assert x.name != y.name


@pytest.mark.parametrize("func", [lambda x, y: x + y, lambda x, y, block_info: x + y])
def test_map_blocks_optimize_blockwise(func):
    # Check that map_blocks layers can merge with elementwise layers
    base = [da.full((1,), i, chunks=1) for i in range(4)]
    a = base[0] + base[1]
    b = da.map_blocks(func, a, base[2], dtype=np.int8)
    c = b + base[3]
    dsk = c.__dask_graph__()
    optimized = optimize_blockwise(dsk)

    # Everything should be fused into a single layer.
    # If the lambda includes block_info, there will be two layers.
    assert len(optimized.layers) == max(len(dsk.layers) - 7, 1)


def test_repr():
    d = da.ones((4, 4), chunks=(2, 2))
    assert key_split(d.name) in repr(d)
    assert str(d.shape) in repr(d)
    assert str(d.dtype) in repr(d)
    d = da.ones((4000, 4), chunks=(4, 2))
    assert len(str(d)) < 1000


def test_repr_meta():
    d = da.ones((4, 4), chunks=(2, 2))
    assert "chunktype=numpy.ndarray" in repr(d)

    # Test non-numpy meta
    sparse = pytest.importorskip("sparse")
    s = d.map_blocks(sparse.COO)
    assert "chunktype=sparse.COO" in repr(s)


def test_repr_html_array_highlevelgraph():
    pytest.importorskip("jinja2")
    x = da.ones((9, 9), chunks=(3, 3)).T[0:4, 0:4]
    hg = x.dask
    assert xml.etree.ElementTree.fromstring(hg._repr_html_()) is not None
    for layer in hg.layers.values():
        assert xml.etree.ElementTree.fromstring(layer._repr_html_()) is not None


def test_slicing_with_ellipsis():
    x = np.arange(256).reshape((4, 4, 4, 4))
    d = da.from_array(x, chunks=((2, 2, 2, 2)))

    assert_eq(d[..., 1], x[..., 1])
    assert_eq(d[0, ..., 1], x[0, ..., 1])


def test_slicing_with_ndarray():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=((4, 4)))

    assert_eq(d[np.arange(8)], x)
    assert_eq(d[np.ones(8, dtype=bool)], x)
    assert_eq(d[np.array([1])], x[[1]])
    assert_eq(d[np.array([True, False, True] + [False] * 5)], x[[0, 2]])


def test_slicing_flexible_type():
    a = np.array([["a", "b"], ["c", "d"]])
    b = da.from_array(a, 2)

    assert_eq(a[:, 0], b[:, 0])


def test_slicing_with_object_dtype():
    # https://github.com/dask/dask/issues/6892
    d = da.from_array(np.array(["a", "b"], dtype=object), chunks=(1,))
    assert d.dtype == d[(0,)].dtype


def test_dtype():
    d = da.ones((4, 4), chunks=(2, 2))

    assert d.dtype == d.compute().dtype
    assert (d * 1.0).dtype == (d + 1.0).compute().dtype
    assert d.sum().dtype == d.sum().compute().dtype  # no shape


def test_blockdims_from_blockshape():
    assert blockdims_from_blockshape((10, 10), (4, 3)) == ((4, 4, 2), (3, 3, 3, 1))
    pytest.raises(TypeError, lambda: blockdims_from_blockshape((10,), None))
    assert blockdims_from_blockshape((1e2, 3), [1e1, 3]) == ((10,) * 10, (3,))
    assert blockdims_from_blockshape((np.int8(10),), (5,)) == ((5, 5),)


def test_coerce():
    d0 = da.from_array(np.array(1), chunks=(1,))
    d1 = da.from_array(np.array([1]), chunks=(1,))
    with dask.config.set(scheduler="sync"):
        for d in d0, d1:
            assert bool(d) is True
            assert int(d) == 1
            assert float(d) == 1.0
            assert complex(d) == complex(1)

    a2 = np.arange(2)
    d2 = da.from_array(a2, chunks=(2,))
    for func in (int, float, complex):
        pytest.raises(TypeError, lambda func=func: func(d2))


def test_bool():
    arr = np.arange(100).reshape((10, 10))
    darr = da.from_array(arr, chunks=(10, 10))
    with pytest.raises(ValueError):
        bool(darr)
        bool(darr == darr)


def test_store_kwargs():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    called = [False]

    def get_func(*args, **kwargs):
        assert kwargs.pop("foo") == "test kwarg"
        r = dask.get(*args, **kwargs)
        called[0] = True
        return r

    called[0] = False
    at = np.zeros(shape=(10, 10))
    store([a], [at], scheduler=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = np.zeros(shape=(10, 10))
    a.store(at, scheduler=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = np.zeros(shape=(10, 10))
    store([a], [at], scheduler=get_func, return_stored=True, foo="test kwarg")
    assert called[0]


def test_store_delayed_target():
    from dask.delayed import delayed

    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    # empty buffers to be used as targets
    targs = {}

    def make_target(key):
        a = np.empty((4, 4))
        targs[key] = a
        return a

    # delayed calls to these targets
    atd = delayed(make_target)("at")
    btd = delayed(make_target)("bt")

    # test not keeping result
    st = store([a, b], [atd, btd])

    at = targs["at"]
    bt = targs["bt"]

    assert st is None
    assert_eq(at, a)
    assert_eq(bt, b)

    # test keeping result
    for st_compute in [False, True]:
        targs.clear()

        st = store([a, b], [atd, btd], return_stored=True, compute=st_compute)
        if st_compute:
            assert all(not any(dask.core.get_deps(e.dask)[0].values()) for e in st)

        st = dask.compute(*st)

        at = targs["at"]
        bt = targs["bt"]

        assert st is not None
        assert isinstance(st, tuple)
        assert all([isinstance(v, np.ndarray) for v in st])
        assert_eq(at, a)
        assert_eq(bt, b)
        assert_eq(st[0], a)
        assert_eq(st[1], b)

        pytest.raises(ValueError, lambda at=at, bt=bt: store([a], [at, bt]))
        pytest.raises(ValueError, lambda at=at: store(at, at))
        pytest.raises(ValueError, lambda at=at, bt=bt: store([at, bt], [at, bt]))


def test_store():
    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.empty(shape=(4, 4))
    bt = np.empty(shape=(4, 4))

    st = store([a, b], [at, bt])
    assert st is None
    assert (at == 2).all()
    assert (bt == 3).all()

    pytest.raises(ValueError, lambda: store([a], [at, bt]))
    pytest.raises(ValueError, lambda: store(at, at))
    pytest.raises(ValueError, lambda: store([at, bt], [at, bt]))


def test_store_regions():
    d = da.ones((4, 4, 4), dtype=int, chunks=(2, 2, 2))
    a, b = d + 1, d + 2
    a = a[:, 1:, :].astype(float)

    region = (slice(None, None, 2), slice(None), [1, 2, 4, 5])

    # Single region:
    at = np.zeros(shape=(8, 3, 6))
    bt = np.zeros(shape=(8, 4, 6))
    v = store([a, b], [at, bt], regions=region, compute=False)
    assert isinstance(v, Delayed)
    assert (at == 0).all() and (bt[region] == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at[region] == 2).all() and (bt[region] == 3).all()
    assert not (bt == 3).all() and not (bt == 0).all()
    assert not (at == 2).all() and not (at == 0).all()

    # Multiple regions:
    at = np.zeros(shape=(8, 3, 6))
    bt = np.zeros(shape=(8, 4, 6))
    v = store([a, b], [at, bt], regions=[region, region], compute=False)
    assert isinstance(v, Delayed)
    assert (at == 0).all() and (bt[region] == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at[region] == 2).all() and (bt[region] == 3).all()
    assert not (bt == 3).all() and not (bt == 0).all()
    assert not (at == 2).all() and not (at == 0).all()

    # Single region (keep result):
    for st_compute in [False, True]:
        at = np.zeros(shape=(8, 3, 6))
        bt = np.zeros(shape=(8, 4, 6))
        v = store(
            [a, b], [at, bt], regions=region, compute=st_compute, return_stored=True
        )
        assert isinstance(v, tuple)
        assert all([isinstance(e, da.Array) for e in v])
        if st_compute:
            assert all(not any(dask.core.get_deps(e.dask)[0].values()) for e in v)
        else:
            assert (at == 0).all() and (bt[region] == 0).all()

        ar, br = v
        assert ar.dtype == a.dtype
        assert br.dtype == b.dtype
        assert ar.shape == a.shape
        assert br.shape == b.shape
        assert ar.chunks == a.chunks
        assert br.chunks == b.chunks

        ar, br = da.compute(ar, br)
        assert (at[region] == 2).all() and (bt[region] == 3).all()
        assert not (bt == 3).all() and not (bt == 0).all()
        assert not (at == 2).all() and not (at == 0).all()
        assert (br == 3).all()
        assert (ar == 2).all()

    # Multiple regions (keep result):
    for st_compute in [False, True]:
        at = np.zeros(shape=(8, 3, 6))
        bt = np.zeros(shape=(8, 4, 6))
        v = store(
            [a, b],
            [at, bt],
            regions=[region, region],
            compute=st_compute,
            return_stored=True,
        )
        assert isinstance(v, tuple)
        assert all([isinstance(e, da.Array) for e in v])
        if st_compute:
            assert all(not any(dask.core.get_deps(e.dask)[0].values()) for e in v)
        else:
            assert (at == 0).all() and (bt[region] == 0).all()

        ar, br = v
        assert ar.dtype == a.dtype
        assert br.dtype == b.dtype
        assert ar.shape == a.shape
        assert br.shape == b.shape
        assert ar.chunks == a.chunks
        assert br.chunks == b.chunks

        ar, br = da.compute(ar, br)
        assert (at[region] == 2).all() and (bt[region] == 3).all()
        assert not (bt == 3).all() and not (bt == 0).all()
        assert not (at == 2).all() and not (at == 0).all()
        assert (br == 3).all()
        assert (ar == 2).all()


def test_store_compute_false():
    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.zeros(shape=(4, 4))
    bt = np.zeros(shape=(4, 4))

    v = store([a, b], [at, bt], compute=False)
    assert isinstance(v, Delayed)

    # You need a well-formed HighLevelgraph for e.g. dask.graph_manipulation.bind
    for layer in v.__dask_layers__():
        assert layer in v.dask.layers

    assert (at == 0).all() and (bt == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at == 2).all() and (bt == 3).all()

    at = np.zeros(shape=(4, 4))
    bt = np.zeros(shape=(4, 4))
    dat, dbt = store([a, b], [at, bt], compute=False, return_stored=True)
    assert isinstance(dat, Array) and isinstance(dbt, Array)
    assert (at == 0).all() and (bt == 0).all()
    assert (dat.compute() == at).all() and (dbt.compute() == bt).all()
    assert (at == 2).all() and (bt == 3).all()

    at = np.zeros(shape=(4, 4))
    bt = np.zeros(shape=(4, 4))
    dat, dbt = store(
        [a, b], [at, bt], compute=False, return_stored=True, load_stored=False
    )
    assert isinstance(dat, Array) and isinstance(dbt, Array)
    assert (at == 0).all() and (bt == 0).all()
    dask.compute(dat, dbt)
    assert (at == 2).all() and (bt == 3).all()


def test_store_nocompute_regions():
    x = da.ones(10, chunks=1)
    y = np.zeros((2, 10))
    d1 = da.store(x, y, regions=(0,), compute=False)
    d2 = da.store(x, y, regions=(1,), compute=False)
    assert d1.key != d2.key


class ThreadSafetyError(Exception):
    pass


class NonthreadSafeStore:
    def __init__(self):
        self.in_use = False

    def __setitem__(self, key, value):
        if self.in_use:
            raise ThreadSafetyError()
        self.in_use = True
        time.sleep(0.001)
        self.in_use = False


class ThreadSafeStore:
    def __init__(self):
        self.concurrent_uses = 0
        self.max_concurrent_uses = 0

    def __setitem__(self, key, value):
        self.concurrent_uses += 1
        self.max_concurrent_uses = max(self.concurrent_uses, self.max_concurrent_uses)
        time.sleep(0.01)
        self.concurrent_uses -= 1


class CounterLock:
    def __init__(self, *args, **kwargs):
        self.lock = Lock(*args, **kwargs)

        self.acquire_count = 0
        self.release_count = 0

    def acquire(self, *args, **kwargs):
        self.acquire_count += 1
        return self.lock.acquire(*args, **kwargs)

    def release(self, *args, **kwargs):
        self.release_count += 1
        return self.lock.release(*args, **kwargs)


def test_store_locks():
    _Lock = type(Lock())
    d = da.ones((10, 10), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.zeros(shape=(10, 10))
    bt = np.zeros(shape=(10, 10))

    lock = Lock()
    v = store([a, b], [at, bt], compute=False, lock=lock)
    assert isinstance(v, Delayed)
    dsk = v.dask
    locks = {
        vv
        for v in dsk.values()
        for vv in (v.args if isinstance(v, Task) else v)
        if isinstance(vv, _Lock)
    }
    assert locks == {lock}

    # Ensure same lock applies over multiple stores
    at = NonthreadSafeStore()
    v = store([a, b], [at, at], lock=lock, scheduler="threads", num_workers=10)
    assert v is None

    # Don't assume thread safety by default
    at = NonthreadSafeStore()
    assert store(a, at, scheduler="threads", num_workers=10) is None
    assert a.store(at, scheduler="threads", num_workers=10) is None

    # Ensure locks can be removed
    at = ThreadSafeStore()
    for i in range(10):
        st = a.store(at, lock=False, scheduler="threads", num_workers=10)
        assert st is None
        if at.max_concurrent_uses > 1:
            break
        if i == 9:
            assert False

    # Verify number of lock calls
    nchunks = sum(math.prod(map(len, e.chunks)) for e in (a, b))
    for c in (False, True):
        at = np.zeros(shape=(10, 10))
        bt = np.zeros(shape=(10, 10))
        lock = CounterLock()

        v = store([a, b], [at, bt], lock=lock, compute=c, return_stored=True)
        assert all(isinstance(e, Array) for e in v)

        da.compute(v)

        # When `return_stored=True` and `compute=False`,
        # the lock should be acquired only once for store and load steps
        # as they are fused together into one step.
        assert lock.acquire_count == lock.release_count
        if c:
            assert lock.acquire_count == 2 * nchunks
        else:
            assert lock.acquire_count == nchunks


def test_store_method_return():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    for compute in [False, True]:
        for return_stored in [False, True]:
            at = np.zeros(shape=(10, 10))
            r = a.store(
                at, scheduler="threads", compute=compute, return_stored=return_stored
            )

            if return_stored:
                assert isinstance(r, Array)
            elif compute:
                assert r is None
            else:
                assert isinstance(r, Delayed)


@pytest.mark.xfail(reason="can't lock with multiprocessing")
def test_store_multiprocessing_lock():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    at = np.zeros(shape=(10, 10))
    st = a.store(at, scheduler="processes", num_workers=10)
    assert st is None


@pytest.mark.parametrize("return_stored", [False, True])
@pytest.mark.parametrize("delayed_target", [False, True])
def test_store_deterministic_keys(return_stored, delayed_target):
    a = da.ones((10, 10), chunks=(2, 2))
    at = np.zeros(shape=(10, 10))
    if delayed_target:
        at = delayed(at)
    st1 = a.store(at, return_stored=return_stored, compute=False)
    st2 = a.store(at, return_stored=return_stored, compute=False)
    assert st1.dask.keys() == st2.dask.keys()


def test_to_hdf5():
    h5py = pytest.importorskip("h5py")
    x = da.ones((4, 4), chunks=(2, 2))
    y = da.ones(4, chunks=2, dtype="i4")

    with tmpfile(".hdf5") as fn:
        x.to_hdf5(fn, "/x")
        with h5py.File(fn, mode="r+") as f:
            d = f["/x"]

            assert_eq(d[:], x)
            assert d.chunks == (2, 2)

    with tmpfile(".hdf5") as fn:
        x.to_hdf5(fn, "/x", chunks=None)
        with h5py.File(fn, mode="r+") as f:
            d = f["/x"]

            assert_eq(d[:], x)
            assert d.chunks is None

    with tmpfile(".hdf5") as fn:
        x.to_hdf5(fn, "/x", chunks=(1, 1))
        with h5py.File(fn, mode="r+") as f:
            d = f["/x"]

            assert_eq(d[:], x)
            assert d.chunks == (1, 1)

    with tmpfile(".hdf5") as fn:
        da.to_hdf5(fn, {"/x": x, "/y": y})

        with h5py.File(fn, mode="r+") as f:
            assert_eq(f["/x"][:], x)
            assert f["/x"].chunks == (2, 2)
            assert_eq(f["/y"][:], y)
            assert f["/y"].chunks == (2,)


def test_to_dask_dataframe():
    pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")
    a = da.ones((4,), chunks=(2,))
    d = a.to_dask_dataframe()
    assert isinstance(d, dd.Series)

    a = da.ones((4, 4), chunks=(2, 2))
    d = a.to_dask_dataframe()
    assert isinstance(d, dd.DataFrame)


def test_np_array_with_zero_dimensions():
    d = da.ones((4, 4), chunks=(2, 2))
    assert_eq(np.array(d.sum()), np.array(d.compute().sum()))


def test_dtype_complex():
    x = np.arange(24).reshape((4, 6)).astype("f4")
    y = np.arange(24).reshape((4, 6)).astype("i8")
    z = np.arange(24).reshape((4, 6)).astype("i2")

    a = da.from_array(x, chunks=(2, 3))
    b = da.from_array(y, chunks=(2, 3))
    c = da.from_array(z, chunks=(2, 3))

    def assert_eq(a, b):
        return isinstance(a, np.dtype) and isinstance(b, np.dtype) and str(a) == str(b)

    assert_eq(a.dtype, x.dtype)
    assert_eq(b.dtype, y.dtype)

    assert_eq((a + 1).dtype, (x + 1).dtype)
    assert_eq((a + b).dtype, (x + y).dtype)
    assert_eq(a.T.dtype, x.T.dtype)
    assert_eq(a[:3].dtype, x[:3].dtype)
    assert_eq((a.dot(b.T)).dtype, (x.dot(y.T)).dtype)

    assert_eq(stack([a, b]).dtype, np.vstack([x, y]).dtype)
    assert_eq(concatenate([a, b]).dtype, np.concatenate([x, y]).dtype)

    assert_eq(b.std().dtype, y.std().dtype)
    assert_eq(c.sum().dtype, z.sum().dtype)
    assert_eq(a.min().dtype, a.min().dtype)
    assert_eq(b.std().dtype, b.std().dtype)
    assert_eq(a.argmin(axis=0).dtype, a.argmin(axis=0).dtype)

    assert_eq(da.sin(c).dtype, np.sin(z).dtype)
    assert_eq(da.exp(b).dtype, np.exp(y).dtype)
    assert_eq(da.floor(a).dtype, np.floor(x).dtype)
    assert_eq(da.isnan(b).dtype, np.isnan(y).dtype)
    with contextlib.suppress(ImportError):
        assert da.isnull(b).dtype == "bool"
        assert da.notnull(b).dtype == "bool"

    x = np.array([("a", 1)], dtype=[("text", "S1"), ("numbers", "i4")])
    d = da.from_array(x, chunks=(1,))

    assert_eq(d["text"].dtype, x["text"].dtype)
    assert_eq(d[["numbers", "text"]].dtype, x[["numbers", "text"]].dtype)


def test_astype():
    x = np.ones((5, 5), dtype="f8")
    d = da.from_array(x, chunks=(2, 2))

    assert d.astype("i8").dtype == "i8"
    assert_eq(d.astype("i8"), x.astype("i8"))
    assert same_keys(d.astype("i8"), d.astype("i8"))

    with pytest.raises(TypeError):
        d.astype("i8", casting="safe")

    with pytest.raises(TypeError):
        d.astype("i8", not_a_real_kwarg="foo")

    # smoketest with kwargs
    assert_eq(d.astype("i8", copy=False), x.astype("i8", copy=False))

    # Check it's a noop
    assert d.astype("f8") is d


def test_astype_gh1151():
    a = np.arange(5).astype(np.int32)
    b = da.from_array(a, (1,))
    assert_eq(a.astype(np.int16), b.astype(np.int16))


def test_astype_gh9318():
    # `order`` kwarg in `astype` should not cause an error
    a = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]], order="C")
    b = da.from_array(a, chunks=(2, 2))
    result_a = a.astype(float, order="F")
    result_b = b.astype(float, order="F")  # if no error at this line, pytest passes
    assert_eq(result_a, result_b)  # won't check the order matches, but checks results


@pytest.mark.xfail(reason="Github issue https://github.com/dask/dask/issues/9316")
def test_astype_gh9316():
    # Issue https://github.com/dask/dask/issues/9316
    # Can be combined with test_astype_gh9318 above when XFAIL marker is removed
    a = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]], order="C")
    b = da.from_array(a, chunks=(2, 2))
    result_a = a.astype(float, order="F")
    result_b = b.astype(float, order="F")
    result_b = result_b.compute()
    assert result_a.flags.c_contiguous == result_b.flags.c_contiguous
    assert result_a.flags.f_contiguous == result_b.flags.f_contiguous


def test_arithmetic():
    x = np.arange(5).astype("f4") + 2
    y = np.arange(5).astype("i8") + 2
    z = np.arange(5).astype("i4") + 2
    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))
    c = da.from_array(z, chunks=(2,))
    assert_eq(a + b, x + y)
    assert_eq(a * b, x * y)
    assert_eq(a - b, x - y)
    assert_eq(a / b, x / y)
    assert_eq(b & b, y & y)
    assert_eq(b | b, y | y)
    assert_eq(b ^ b, y ^ y)
    assert_eq(a // b, x // y)
    assert_eq(a**b, x**y)
    assert_eq(a % b, x % y)
    assert_eq(a > b, x > y)
    assert_eq(a < b, x < y)
    assert_eq(a >= b, x >= y)
    assert_eq(a <= b, x <= y)
    assert_eq(a == b, x == y)
    assert_eq(a != b, x != y)

    assert_eq(a + 2, x + 2)
    assert_eq(a * 2, x * 2)
    assert_eq(a - 2, x - 2)
    assert_eq(a / 2, x / 2)
    assert_eq(b & True, y & True)
    assert_eq(b | True, y | True)
    assert_eq(b ^ True, y ^ True)
    assert_eq(a // 2, x // 2)
    assert_eq(a**2, x**2)
    assert_eq(a % 2, x % 2)
    assert_eq(a > 2, x > 2)
    assert_eq(a < 2, x < 2)
    assert_eq(a >= 2, x >= 2)
    assert_eq(a <= 2, x <= 2)
    assert_eq(a == 2, x == 2)
    assert_eq(a != 2, x != 2)

    assert_eq(2 + b, 2 + y)
    assert_eq(2 * b, 2 * y)
    assert_eq(2 - b, 2 - y)
    assert_eq(2 / b, 2 / y)
    assert_eq(True & b, True & y)
    assert_eq(True | b, True | y)
    assert_eq(True ^ b, True ^ y)
    assert_eq(2 // b, 2 // y)
    assert_eq(2**b, 2**y)
    assert_eq(2 % b, 2 % y)
    assert_eq(2 > b, 2 > y)
    assert_eq(2 < b, 2 < y)
    assert_eq(2 >= b, 2 >= y)
    assert_eq(2 <= b, 2 <= y)
    assert_eq(2 == b, 2 == y)
    assert_eq(2 != b, 2 != y)

    assert_eq(-a, -x)
    assert_eq(abs(a), abs(x))
    assert_eq(~(a == b), ~(x == y))
    assert_eq(~(a == b), ~(x == y))

    assert_eq(da.logaddexp(a, b), np.logaddexp(x, y))
    assert_eq(da.logaddexp2(a, b), np.logaddexp2(x, y))
    assert_eq(da.exp(b), np.exp(y))
    assert_eq(da.log(a), np.log(x))
    assert_eq(da.log10(a), np.log10(x))
    assert_eq(da.log1p(a), np.log1p(x))
    assert_eq(da.expm1(b), np.expm1(y))
    assert_eq(da.sqrt(a), np.sqrt(x))
    assert_eq(da.square(a), np.square(x))

    assert_eq(da.sin(a), np.sin(x))
    assert_eq(da.cos(b), np.cos(y))
    assert_eq(da.tan(a), np.tan(x))
    assert_eq(da.arcsin(b / 10), np.arcsin(y / 10))
    assert_eq(da.arccos(b / 10), np.arccos(y / 10))
    assert_eq(da.arctan(b / 10), np.arctan(y / 10))
    assert_eq(da.arctan2(b * 10, a), np.arctan2(y * 10, x))
    assert_eq(da.hypot(b, a), np.hypot(y, x))
    assert_eq(da.sinh(a), np.sinh(x))
    assert_eq(da.cosh(b), np.cosh(y))
    assert_eq(da.tanh(a), np.tanh(x))
    assert_eq(da.arcsinh(b * 10), np.arcsinh(y * 10))
    assert_eq(da.arccosh(b * 10), np.arccosh(y * 10))
    assert_eq(da.arctanh(b / 10), np.arctanh(y / 10))
    assert_eq(da.deg2rad(a), np.deg2rad(x))
    assert_eq(da.rad2deg(a), np.rad2deg(x))

    assert_eq(da.logical_and(a < 1, b < 4), np.logical_and(x < 1, y < 4))
    assert_eq(da.logical_or(a < 1, b < 4), np.logical_or(x < 1, y < 4))
    assert_eq(da.logical_xor(a < 1, b < 4), np.logical_xor(x < 1, y < 4))
    assert_eq(da.logical_not(a < 1), np.logical_not(x < 1))
    assert_eq(da.maximum(a, 5 - a), np.maximum(a, 5 - a))
    assert_eq(da.minimum(a, 5 - a), np.minimum(a, 5 - a))
    assert_eq(da.fmax(a, 5 - a), np.fmax(a, 5 - a))
    assert_eq(da.fmin(a, 5 - a), np.fmin(a, 5 - a))

    assert_eq(da.isreal(a + 1j * b), np.isreal(x + 1j * y))
    assert_eq(da.iscomplex(a + 1j * b), np.iscomplex(x + 1j * y))
    assert_eq(da.isfinite(a), np.isfinite(x))
    assert_eq(da.isinf(a), np.isinf(x))
    assert_eq(da.isnan(a), np.isnan(x))
    assert_eq(da.signbit(a - 3), np.signbit(x - 3))
    assert_eq(da.copysign(a - 3, b), np.copysign(x - 3, y))
    assert_eq(da.nextafter(a - 3, b), np.nextafter(x - 3, y))
    assert_eq(da.ldexp(c, c), np.ldexp(z, z))
    assert_eq(da.fmod(a * 12, b), np.fmod(x * 12, y))
    assert_eq(da.floor(a * 0.5), np.floor(x * 0.5))
    assert_eq(da.ceil(a), np.ceil(x))
    assert_eq(da.trunc(a / 2), np.trunc(x / 2))

    assert_eq(da.degrees(b), np.degrees(y))
    assert_eq(da.radians(a), np.radians(x))

    assert_eq(da.rint(a + 0.3), np.rint(x + 0.3))
    assert_eq(da.fix(a - 2.5), np.fix(x - 2.5))

    assert_eq(da.angle(a + 1j), np.angle(x + 1j))
    assert_eq(da.real(a + 1j), np.real(x + 1j))
    assert_eq((a + 1j).real, np.real(x + 1j))
    assert_eq(da.imag(a + 1j), np.imag(x + 1j))
    assert_eq((a + 1j).imag, np.imag(x + 1j))
    assert_eq(da.conj(a + 1j * b), np.conj(x + 1j * y))
    assert_eq((a + 1j * b).conj(), (x + 1j * y).conj())

    assert_eq(da.clip(b, 1, 4), np.clip(y, 1, 4))
    assert_eq(b.clip(1, 4), y.clip(1, 4))
    assert_eq(da.fabs(b), np.fabs(y))
    assert_eq(da.sign(b - 2), np.sign(y - 2))
    assert_eq(da.absolute(b - 2), np.absolute(y - 2))
    assert_eq(da.absolute(b - 2 + 1j), np.absolute(y - 2 + 1j))

    l1, l2 = da.frexp(a)
    r1, r2 = np.frexp(x)
    assert_eq(l1, r1)
    assert_eq(l2, r2)

    l1, l2 = da.modf(a)
    r1, r2 = np.modf(x)
    assert_eq(l1, r1)
    assert_eq(l2, r2)

    assert_eq(da.around(a, -1), np.around(x, -1))


def test_elemwise_consistent_names():
    a = da.from_array(np.arange(5, dtype="f4"), chunks=(2,))
    b = da.from_array(np.arange(5, dtype="f4"), chunks=(2,))
    assert same_keys(a + b, a + b)
    assert same_keys(a + 2, a + 2)
    assert same_keys(da.exp(a), da.exp(a))
    assert same_keys(da.exp(a, dtype="f8"), da.exp(a, dtype="f8"))
    assert same_keys(da.maximum(a, b), da.maximum(a, b))


def test_optimize():
    x = np.arange(5).astype("f4")
    a = da.from_array(x, chunks=(2,))
    expr = a[1:4] + 1
    result = optimize(expr.dask, expr.__dask_keys__())
    assert isinstance(result, dict)
    assert all(key in result for key in expr.__dask_keys__())


def test_slicing_with_non_ndarrays():
    class ARangeSlice:
        dtype = np.dtype("i8")
        ndim = 1

        def __init__(self, start, stop):
            self.start = start
            self.stop = stop

        def __array__(self):
            return np.arange(self.start, self.stop)

    class ARangeSlicable:
        dtype = np.dtype("i8")
        ndim = 1

        def __init__(self, n):
            self.n = n

        @property
        def shape(self):
            return (self.n,)

        def __getitem__(self, key):
            return ARangeSlice(key[0].start, key[0].stop)

    x = da.from_array(ARangeSlicable(10), chunks=(4,))

    assert_eq((x + 1).sum(), (np.arange(10, dtype=x.dtype) + 1).sum())


@pytest.mark.filterwarnings("ignore:the matrix subclass")
def test_getter():
    assert type(getter(np.matrix([[1]]), 0)) is np.ndarray
    assert type(getter(np.matrix([[1]]), 0, asarray=False)) is np.matrix
    assert_eq(getter([1, 2, 3, 4, 5], slice(1, 4)), np.array([2, 3, 4]))

    assert_eq(getter(np.arange(5), (None, slice(None, None))), np.arange(5)[None, :])


def test_size():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.size == np.array(x).size
    assert isinstance(x.size, int)


def test_nbytes():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.nbytes == np.array(x).nbytes


def test_itemsize():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.itemsize == 8


def test_Array_normalizes_dtype():
    x = da.ones((3,), chunks=(1,), dtype=int)
    assert isinstance(x.dtype, np.dtype)


@pytest.mark.parametrize("inline_array", [True, False])
def test_from_array_with_lock(inline_array):
    x = np.arange(10)

    class FussyLock(SerializableLock):
        def acquire(self, blocking=True, timeout=-1):
            if self.locked():
                raise RuntimeError("I am locked")
            return super().acquire(blocking, timeout)

    lock = FussyLock()
    d = da.from_array(x, chunks=5, lock=lock, inline_array=inline_array)

    lock.acquire()
    with pytest.raises(RuntimeError):
        d.compute()

    lock.release()
    assert_eq(d, x)

    lock = CounterLock()
    e = da.from_array(x, chunks=5, lock=lock, inline_array=inline_array)

    assert_eq(e, x)
    # Note: the specific counts for composite arithmetic operations can vary
    # significantly based on the complexity of the computation, whether we are inlining,
    # and optimization fusion settings. But for this simple comparison it seems pretty
    # stable.
    assert lock.release_count == 2
    assert lock.acquire_count == 2


class MyArray:
    def __init__(self, x):
        self.x = x
        self.dtype = x.dtype
        self.shape = x.shape
        self.ndim = len(x.shape)

    def __getitem__(self, i):
        return self.x[i]


@pytest.mark.parametrize(
    "x,chunks",
    [
        (np.arange(25).reshape((5, 5)), (5, 5)),
        (np.arange(25).reshape((5, 5)), -1),
        (np.array([[1]]), 1),
        (np.array(1), 1),
    ],
)
@pytest.mark.parametrize("inline_array", [True, False])
def test_from_array_tasks_always_call_getter(x, chunks, inline_array):
    dx = da.from_array(
        MyArray(x), chunks=chunks, asarray=False, inline_array=inline_array
    )
    assert_eq(x, dx)


@pytest.mark.parametrize(
    "x",
    [
        np.array([[1, 2], [3, 4]]),
        np.ma.array([[1, 2], [3, 4]], mask=[[True, False], [False, False]]),
        np.ma.array([1], mask=[True]),
        np.ma.array([1.5], mask=[True]),
        np.ma.array(1, mask=True),
        np.ma.array(1.5, mask=True),
    ],
)
def test_from_array_ndarray_onechunk(x):
    """ndarray with a single chunk produces a minimal single key dict"""
    dx = da.from_array(x, chunks=-1)
    assert_eq(x, dx)
    assert len(dx.dask) == 1
    assert not dx.dask[(dx.name,) + (0,) * dx.ndim] is x
    assert_eq(dx.dask[(dx.name,) + (0,) * dx.ndim], x)


def test_from_array_ndarray_getitem():
    """For ndarray, don't use getter / getter_nofancy; use the cleaner
    operator.getitem"""
    x = np.array([[1, 2], [3, 4]])
    dx = da.from_array(x, chunks=(1, 2))
    assert_eq(x, dx)
    assert (dx.dask[dx.name, 0, 0] == np.array([[1, 2]])).all()


@pytest.mark.parametrize("x", [[1, 2], (1, 2), memoryview(b"abc")])
def test_from_array_list(x):
    """Lists, tuples, and memoryviews are automatically converted to ndarray"""
    dx = da.from_array(x, chunks=-1)
    assert_eq(np.array(x), dx)
    assert isinstance(dx.dask[dx.name, 0], np.ndarray)

    dx = da.from_array(x, chunks=1)
    assert_eq(np.array(x), dx)
    assert dx.dask[dx.name, 0][0] == x[0]


# On MacOS Python 3.9, the order of the np.ScalarType tuple randomly changes across
# interpreter restarts, thus causing pytest-xdist failures; setting PYTHONHASHSEED does
# not help
@pytest.mark.parametrize(
    "type_", sorted((t for t in np.ScalarType if t is not memoryview), key=str)
)
def test_from_array_scalar(type_):
    """Python and numpy scalars are automatically converted to ndarray"""
    if type_ == np.datetime64:
        x = np.datetime64("2000-01-01")
    else:
        x = type_(1)

    dx = da.from_array(x, chunks=-1)
    assert_eq(np.array(x), dx)
    assert isinstance(
        dx.dask[dx.name,],
        np.ndarray,
    )


@pytest.mark.parametrize("asarray,cls", [(True, np.ndarray), (False, np.matrix)])
@pytest.mark.parametrize("inline_array", [True, False])
@pytest.mark.filterwarnings("ignore:the matrix subclass")
def test_from_array_no_asarray(asarray, cls, inline_array):
    def assert_chunks_are_of_type(x):
        chunks = compute_as_if_collection(Array, x.dask, x.__dask_keys__())
        # If it's a tuple of tuples we want to concat, but if it's a tuple
        # of 1d arrays, we just want to iterate directly
        for c in concat(chunks) if isinstance(chunks[0], tuple) else chunks:
            assert type(c) is cls

    x = np.matrix(np.arange(100).reshape((10, 10)))
    dx = da.from_array(x, chunks=(5, 5), asarray=asarray, inline_array=inline_array)
    assert_chunks_are_of_type(dx)
    assert_chunks_are_of_type(dx[0:5])
    assert_chunks_are_of_type(dx[0:5][:, 0])


@pytest.mark.parametrize("wrap", [True, False])
@pytest.mark.parametrize("inline_array", [True, False])
def test_from_array_getitem(wrap, inline_array):
    x = np.arange(10)
    called = False

    def my_getitem(a, ind):
        nonlocal called
        called = True
        return a[ind]

    xx = MyArray(x) if wrap else x
    y = da.from_array(xx, chunks=(5,), getitem=my_getitem, inline_array=inline_array)

    assert_eq(x, y)
    # If we have a raw numpy array we eagerly slice, so custom getters
    # are not called.
    assert called is wrap


def test_from_array_minus_one():
    x = np.arange(10)
    y = da.from_array(x, -1)
    assert y.chunks == ((10,),)
    assert_eq(x, y)


@pytest.mark.parametrize("chunks", [-1, 2])
def test_array_copy_noop(chunks):
    # Regression test for https://github.com/dask/dask/issues/9533
    # Which is a revert of the solution for https://github.com/dask/dask/issues/3751
    x = np.arange(10)
    y = da.from_array(x, chunks=chunks)
    y_c = y.copy()
    assert y.name == y_c.name


def test_from_array_dask_array():
    x = np.array([[1, 2], [3, 4]])
    dx = da.from_array(x, chunks=(1, 2))
    with pytest.raises(ValueError):
        da.from_array(dx)


def test_from_array_dask_collection_warns():
    class CustomCollection(np.ndarray):
        def __dask_graph__(self):
            return {"bar": 1}

    x = CustomCollection([1, 2, 3])
    with pytest.warns(UserWarning):
        da.from_array(x)

    # Ensure da.array warns too
    with pytest.warns(UserWarning):
        da.array(x)


def test_from_array_inline():
    class MyArray(np.ndarray):
        pass

    a = np.array([1, 2, 3]).view(MyArray)
    dsk = dict(da.from_array(a, name="my-array", inline_array=False).dask)
    assert not dsk["original-my-array"] is a
    assert_eq(dsk["original-my-array"], a)

    dsk = dict(da.from_array(a, name="my-array", inline_array=True).dask)
    assert "original-my-array" not in dsk


@pytest.mark.parametrize("asarray", [da.asarray, da.asanyarray])
def test_asarray(asarray):
    assert_eq(asarray([1, 2, 3]), np.asarray([1, 2, 3]))

    x = asarray([1, 2, 3])
    assert asarray(x) is x

    y = [x[0], 2, x[2]]
    assert_eq(asarray(y), x)


@pytest.mark.parametrize("asarray", [da.asarray, da.asanyarray])
def test_asarray_array_dtype(asarray):
    # test array input
    x = asarray([1, 2])
    assert_eq(asarray(x, dtype=da.float32), np.asarray(x, dtype=np.float32))

    # dask->dask
    x = asarray(x, dtype=da.float64)
    assert x.dtype == da.float64
    x = asarray(x, dtype=da.int32)
    assert x.dtype == da.int32
    x = asarray(x)
    assert x.dtype == da.int32
    # Test explicit null dtype. astype(None) converts to float!
    x = asarray(x, dtype=None)
    assert x.dtype == da.int32

    # non-dask->dask
    x = asarray(np.asarray([1, 2], dtype=np.int8))
    assert x.dtype == da.int8
    x = asarray(np.asarray([1, 2], dtype=np.int8), dtype=None)
    assert x.dtype == da.int8
    x = asarray(np.asarray([1, 2], dtype=np.int8), dtype=da.float64)
    assert x.dtype == da.float64


@pytest.mark.parametrize("asarray", [da.asarray, da.asanyarray])
def test_asarray_dask_dataframe(asarray):
    # https://github.com/dask/dask/issues/3885
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    s = dd.from_pandas(pd.Series([1, 2, 3, 4]), 2)
    result = asarray(s)
    expected = s.values
    assert_eq(result, expected)

    df = s.to_frame(name="s")
    result = asarray(df)
    expected = df.values
    assert_eq(result, expected)


@pytest.mark.parametrize("asarray", [da.asarray, da.asanyarray])
@pytest.mark.parametrize("inline_array", [True, False])
def test_asarray_h5py(asarray, inline_array):
    h5py = pytest.importorskip("h5py")

    with tmpfile(".hdf5") as fn:
        with h5py.File(fn, mode="a") as f:
            d = f.create_dataset("/x", shape=(2, 2), dtype=float)
            x = asarray(d, inline_array=inline_array)

            # Check for the array in the dsk
            dsk = dict(x.dask)
            assert (d in dsk.values()) is not inline_array
            assert not any(isinstance(v, np.ndarray) for v in dsk.values())


def test_asarray_chunks():
    with dask.config.set({"array.chunk-size": "100 B"}):
        x = np.ones(1000)
        d = da.asarray(x)
        assert d.npartitions > 1


@pytest.mark.filterwarnings("ignore:the matrix subclass")
def test_asanyarray():
    x = np.matrix([1, 2, 3])
    dx = da.asanyarray(x)
    assert dx.numblocks == (1, 1)
    chunks = compute_as_if_collection(Array, dx.dask, dx.__dask_keys__())
    assert isinstance(chunks[0][0], np.matrix)
    assert da.asanyarray(dx) is dx


def test_asanyarray_dataframe():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    df = pd.DataFrame({"x": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=2)

    x = np.asanyarray(df)
    dx = da.asanyarray(ddf)
    assert isinstance(dx, da.Array)

    assert_eq(x, dx)

    x = np.asanyarray(df.x)
    dx = da.asanyarray(ddf.x)
    assert isinstance(dx, da.Array)

    assert_eq(x, dx)


def test_asanyarray_datetime64():
    x = np.array(["2000-01-01"], dtype="datetime64")
    dx = da.asanyarray(x)
    assert isinstance(dx, da.Array)
    assert_eq(x, dx)


def test_from_func():
    x = np.arange(10)
    f = lambda n: n * x
    d = from_func(f, (10,), x.dtype, kwargs={"n": 2})

    assert d.shape == x.shape
    assert d.dtype == x.dtype
    assert_eq(d, 2 * x)
    assert same_keys(d, from_func(f, (10,), x.dtype, kwargs={"n": 2}))


def test_concatenate3_2():
    x = np.array([1, 2])
    assert_eq(concatenate3([x, x, x]), np.array([1, 2, 1, 2, 1, 2]))

    x = np.array([[1, 2]])
    assert (
        concatenate3([[x, x, x], [x, x, x]])
        == np.array([[1, 2, 1, 2, 1, 2], [1, 2, 1, 2, 1, 2]])
    ).all()

    assert (
        concatenate3([[x, x], [x, x], [x, x]])
        == np.array([[1, 2, 1, 2], [1, 2, 1, 2], [1, 2, 1, 2]])
    ).all()

    x = np.arange(12).reshape((2, 2, 3))
    assert_eq(
        concatenate3([[[x, x, x], [x, x, x]], [[x, x, x], [x, x, x]]]),
        np.array(
            [
                [
                    [0, 1, 2, 0, 1, 2, 0, 1, 2],
                    [3, 4, 5, 3, 4, 5, 3, 4, 5],
                    [0, 1, 2, 0, 1, 2, 0, 1, 2],
                    [3, 4, 5, 3, 4, 5, 3, 4, 5],
                ],
                [
                    [6, 7, 8, 6, 7, 8, 6, 7, 8],
                    [9, 10, 11, 9, 10, 11, 9, 10, 11],
                    [6, 7, 8, 6, 7, 8, 6, 7, 8],
                    [9, 10, 11, 9, 10, 11, 9, 10, 11],
                ],
                [
                    [0, 1, 2, 0, 1, 2, 0, 1, 2],
                    [3, 4, 5, 3, 4, 5, 3, 4, 5],
                    [0, 1, 2, 0, 1, 2, 0, 1, 2],
                    [3, 4, 5, 3, 4, 5, 3, 4, 5],
                ],
                [
                    [6, 7, 8, 6, 7, 8, 6, 7, 8],
                    [9, 10, 11, 9, 10, 11, 9, 10, 11],
                    [6, 7, 8, 6, 7, 8, 6, 7, 8],
                    [9, 10, 11, 9, 10, 11, 9, 10, 11],
                ],
            ]
        ),
    )


def test_map_blocks3():
    x = np.arange(10)
    y = np.arange(10) * 2

    d = da.from_array(x, chunks=5)
    e = da.from_array(y, chunks=5)

    assert_eq(
        da.core.map_blocks(lambda a, b: a + 2 * b, d, e, dtype=d.dtype), x + 2 * y
    )

    z = np.arange(100).reshape((10, 10))
    f = da.from_array(z, chunks=5)

    func = lambda a, b: a + 2 * b
    res = da.core.map_blocks(func, d, f, dtype=d.dtype)
    assert_eq(res, x + 2 * z)
    assert same_keys(da.core.map_blocks(func, d, f, dtype=d.dtype), res)

    assert_eq(da.map_blocks(func, f, d, dtype=d.dtype), z + 2 * x)


def test_from_array_with_missing_chunks():
    x = np.random.default_rng().standard_normal((2, 4, 3))
    d = da.from_array(x, chunks=(None, 2, None))
    assert d.chunks == da.from_array(x, chunks=(2, 2, 3)).chunks


@pytest.mark.parametrize("func", [normalize_chunks, normalize_chunks_cached])
def test_normalize_chunks(func):
    assert func(3, (4, 6)) == ((3, 1), (3, 3))
    assert func(((3, 3), (8,)), (6, 8)) == ((3, 3), (8,))
    assert func((4, 5), (9,)) == ((4, 5),)
    assert func((4, 5), (9, 9)) == ((4, 4, 1), (5, 4))
    assert func(-1, (5, 5)) == ((5,), (5,))
    assert func((3, -1), (5, 5)) == ((3, 2), (5,))
    assert func((3, None), (5, 5)) == ((3, 2), (5,))
    if func is normalize_chunks:
        assert func({0: 3}, (5, 5)) == ((3, 2), (5,))
        assert func([[2, 2], [3, 3]]) == ((2, 2), (3, 3))
    assert func(10, (30, 5)) == ((10, 10, 10), (5,))
    assert func((), (0, 0)) == ((0,), (0,))
    assert func(-1, (0, 3)) == ((0,), (3,))
    assert func(((float("nan"),),)) == ((np.nan,),)

    assert func("auto", shape=(20,), limit=5, dtype="uint8") == ((5, 5, 5, 5),)
    assert func(("auto", None), (5, 5), dtype=int) == ((5,), (5,))

    with pytest.raises(ValueError):
        func(((10,),), (11,))
    with pytest.raises(ValueError):
        func(((5,), (5,)), (5,))


def test_single_element_tuple():
    assert normalize_chunks(
        (100, "auto"), (500, 500_000), dtype=np.int64, previous_chunks=((1,), (500,))
    ) == (
        (100,) * 5,
        (
            167_500,
            167_500,
            165_000,
        ),
    )


def test_align_chunks_to_previous_chunks():
    chunks = normalize_chunks(
        "auto", shape=(2000,), previous_chunks=(512,), limit="600 B", dtype=np.uint8
    )
    assert chunks == ((512, 512, 512, 2000 - 512 * 3),)

    chunks = normalize_chunks(
        "auto", shape=(2000,), previous_chunks=(128,), limit="600 B", dtype=np.uint8
    )
    assert chunks == ((512, 512, 512, 2000 - 512 * 3),)

    chunks = normalize_chunks(
        "auto", shape=(2000,), previous_chunks=(512,), limit="1200 B", dtype=np.uint8
    )
    assert chunks == ((1024, 2000 - 1024),)

    chunks = normalize_chunks(
        "auto",
        shape=(3, 10211, 10376),
        previous_chunks=(1, 512, 512),
        limit="1MiB",
        dtype=np.float32,
    )
    assert chunks[0] == (1, 1, 1)
    assert all(c % 512 == 0 for c in chunks[1][:-1])
    assert all(c % 512 == 0 for c in chunks[2][:-1])

    chunks = normalize_chunks(
        "auto",
        shape=(48, 720, 1440),
        previous_chunks=((36, 12), (720,), (1440,)),
        limit=134217728,
        dtype=np.float32,
    )
    assert chunks == ((36, 12), (720,), (1440,))


def test_raise_on_no_chunks():
    x = da.ones(6, chunks=3)
    try:
        Array(x.dask, x.name, chunks=None, dtype=x.dtype, shape=None)
        assert False
    except ValueError as e:
        assert "dask" in str(e)
        assert ".org" in str(e)


def test_chunks_is_immutable():
    x = da.ones(6, chunks=3)
    try:
        x.chunks = 2
        assert False
    except TypeError as e:
        assert "rechunk(2)" in str(e)


def test_raise_on_bad_kwargs():
    x = da.ones(5, chunks=3)
    try:
        da.minimum(x, foo=None)
    except TypeError as e:
        assert "minimum" in str(e)
        assert "foo" in str(e)


def test_long_slice():
    x = np.arange(10000)
    d = da.from_array(x, chunks=1)

    assert_eq(d[8000:8200], x[8000:8200])


def test_h5py_newaxis():
    h5py = pytest.importorskip("h5py")

    with tmpfile("h5") as fn:
        with h5py.File(fn, mode="a") as f:
            x = f.create_dataset("/x", shape=(10, 10), dtype="f8")
            d = da.from_array(x, chunks=(5, 5))
            assert d[None, :, :].compute(scheduler="sync").shape == (1, 10, 10)
            assert d[:, None, :].compute(scheduler="sync").shape == (10, 1, 10)
            assert d[:, :, None].compute(scheduler="sync").shape == (10, 10, 1)
            assert same_keys(d[:, :, None], d[:, :, None])


def test_ellipsis_slicing():
    assert_eq(da.ones(4, chunks=2)[...], np.ones(4))


def test_point_slicing():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    result = d.vindex[[1, 2, 5, 5], [3, 1, 6, 1]]
    assert_eq(result, x[[1, 2, 5, 5], [3, 1, 6, 1]])

    result = d.vindex[[0, 1, 6, 0], [0, 1, 0, 7]]
    assert_eq(result, x[[0, 1, 6, 0], [0, 1, 0, 7]])
    assert same_keys(result, d.vindex[[0, 1, 6, 0], [0, 1, 0, 7]])


def test_point_slicing_with_full_slice():
    from dask.array.core import _get_axis

    x = np.arange(4 * 5 * 6 * 7).reshape((4, 5, 6, 7))
    d = da.from_array(x, chunks=(2, 3, 3, 4))

    inds = [
        [[1, 2, 3], None, [3, 2, 1], [5, 3, 4]],
        [[1, 2, 3], None, [4, 3, 2], None],
        [[1, 2, 3], [3, 2, 1]],
        [[1, 2, 3], [3, 2, 1], [3, 2, 1], [5, 3, 4]],
        [[], [], [], None],
        [np.array([1, 2, 3]), None, np.array([4, 3, 2]), None],
        [None, None, [1, 2, 3], [4, 3, 2]],
        [None, [0, 2, 3], None, [0, 3, 2]],
    ]

    for ind in inds:
        slc = [
            i if isinstance(i, (np.ndarray, list)) else slice(None, None) for i in ind
        ]
        result = d.vindex[tuple(slc)]

        # Rotate the expected result accordingly
        axis = _get_axis(ind)
        expected = x[tuple(slc)]
        expected = expected.transpose(
            [axis] + list(range(axis)) + list(range(axis + 1, expected.ndim))
        )

        assert_eq(result, expected)

        # Always have the first axis be the length of the points
        k = len(next(i for i in ind if isinstance(i, (np.ndarray, list))))
        assert result.shape[0] == k


def test_slice_with_floats():
    d = da.ones((5,), chunks=(3,))
    with pytest.raises(IndexError):
        d[1.5]
    with pytest.raises(IndexError):
        d[0:1.5]
    with pytest.raises(IndexError):
        d[[1, 1.5]]


@pytest.mark.parametrize("dtype", [np.int32, np.int64, np.uint32, np.uint64])
def test_slice_with_integer_types(dtype):
    x = np.arange(10)
    dx = da.from_array(x, chunks=5)
    inds = np.array([0, 3, 6], dtype=dtype)
    assert_eq(dx[inds], x[inds])


@pytest.mark.parametrize("cls", [int, np.int32, np.int64, np.uint32, np.uint64])
def test_index_with_integer_types(cls):
    x = np.arange(10)
    dx = da.from_array(x, chunks=5)
    inds = cls(3)
    assert_eq(dx[inds], x[inds])


def test_vindex_basic():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    # cases where basic and advanced indexing coincide
    result = d.vindex[0]
    assert_eq(result, x[0])

    result = d.vindex[0, 1]
    assert_eq(result, x[0, 1])

    result = d.vindex[[0, 1], ::-1]  # slices last
    assert_eq(result, x[:2, ::-1])


def test_vindex_nd():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    result = d.vindex[[[0, 1], [6, 0]], [[0, 1], [0, 7]]]
    assert_eq(result, x[[[0, 1], [6, 0]], [[0, 1], [0, 7]]])

    result = d.vindex[np.arange(7)[:, None], np.arange(8)[None, :]]
    assert_eq(result, x)

    result = d.vindex[np.arange(7)[None, :], np.arange(8)[:, None]]
    assert_eq(result, x.T)


@pytest.mark.parametrize("size", [0, 1])
def test_vindex_preserve_chunksize(size):
    np_arr = np.random.rand(10_000 * 40).reshape(100, 100, 40)
    arr = da.from_array(np_arr, chunks=(50, 50, 20))
    indices_2d = np.random.choice(np.arange(100), size=(10000 + size, 2))
    idx1 = indices_2d[:, 0]
    idx2 = indices_2d[:, 0]
    result = arr.vindex[idx1, idx2, slice(None)]
    assert result.chunks == (
        (2500, 2500, 2500, 2500) + ((1,) if size else ()),
        (20, 20),
    )
    assert_eq(result, np_arr[idx1, idx2, :])


def test_vindex_negative():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    result = d.vindex[np.array([0, -1])]
    assert_eq(result, x[np.array([0, -1])])


def test_vindex_errors():
    d = da.ones((5, 5, 5), chunks=(3, 3, 3))
    pytest.raises(IndexError, lambda: d.vindex[np.newaxis])
    pytest.raises(IndexError, lambda: d.vindex[[1, 2], [1, 2, 3]])
    pytest.raises(IndexError, lambda: d.vindex[[True] * 5])
    pytest.raises(IndexError, lambda: d.vindex[[0], [5]])
    pytest.raises(IndexError, lambda: d.vindex[[0], [-6]])
    with pytest.raises(IndexError, match="does not support indexing with dask objects"):
        d.vindex[[0], [0], da.array([0])]


def test_vindex_merge():
    from dask.array.core import _vindex_merge

    locations = [1], [2, 0]
    values = [np.array([[1, 2, 3]]), np.array([[10, 20, 30], [40, 50, 60]])]

    assert (
        _vindex_merge(locations, values)
        == np.array([[40, 50, 60], [1, 2, 3], [10, 20, 30]])
    ).all()


def test_vindex_identity():
    rng = da.random.default_rng(42)
    a, b = 10, 20

    x = rng.random(a, chunks=a // 2)
    assert x is x.vindex[:]
    assert x is x.vindex[:a]
    pytest.raises(IndexError, lambda: x.vindex[: a - 1])
    pytest.raises(IndexError, lambda: x.vindex[1:])
    pytest.raises(IndexError, lambda: x.vindex[0:a:2])

    x = rng.random((a, b), chunks=(a // 2, b // 2))
    assert x is x.vindex[:, :]
    assert x is x.vindex[:a, :b]
    pytest.raises(IndexError, lambda: x.vindex[:, : b - 1])
    pytest.raises(IndexError, lambda: x.vindex[:, 1:])
    pytest.raises(IndexError, lambda: x.vindex[:, 0:b:2])


def test_empty_array():
    assert_eq(np.arange(0), da.arange(0, chunks=5))


def test_memmap():
    with tmpfile("npy") as fn_1:
        with tmpfile("npy") as fn_2:
            try:
                x = da.arange(100, chunks=15)
                target = np.memmap(fn_1, shape=x.shape, mode="w+", dtype=x.dtype)

                x.store(target)

                assert_eq(target, x, check_type=False)

                np.save(fn_2, target)

                assert_eq(np.load(fn_2, mmap_mode="r"), x, check_type=False)
            finally:
                target._mmap.close()


def test_to_npy_stack():
    x = np.arange(5 * 10 * 10).reshape((5, 10, 10))
    d = da.from_array(x, chunks=(2, 4, 4))

    with tmpdir() as dirname:
        stackdir = os.path.join(dirname, "test")
        da.to_npy_stack(stackdir, d, axis=0)
        assert os.path.exists(os.path.join(stackdir, "0.npy"))
        assert (np.load(os.path.join(stackdir, "1.npy")) == x[2:4]).all()

        e = da.from_npy_stack(stackdir)
        assert_eq(d, e)


def test_view():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(2, 3))

    assert_eq(x.view(), d.view())
    assert_eq(x.view("i4"), d.view("i4"))
    assert_eq(x.view("i2"), d.view("i2"))
    assert all(isinstance(s, int) for s in d.shape)

    x = np.arange(8, dtype="i1")
    d = da.from_array(x, chunks=(4,))
    assert_eq(x.view("i4"), d.view("i4"))

    with pytest.raises(ValueError):
        x = np.arange(8, dtype="i1")
        d = da.from_array(x, chunks=(3,))
        d.view("i4")

    with pytest.raises(ValueError):
        d.view("i4", order="asdf")


def test_view_fortran():
    x = np.asfortranarray(np.arange(64).reshape((8, 8)))
    d = da.from_array(x, chunks=(2, 3))
    assert_eq(x.T.view("i4").T, d.view("i4", order="F"))
    assert_eq(x.T.view("i2").T, d.view("i2", order="F"))


def test_h5py_tokenize():
    h5py = pytest.importorskip("h5py")
    with tmpfile("hdf5") as fn1:
        with tmpfile("hdf5") as fn2:
            f = h5py.File(fn1, mode="a")
            g = h5py.File(fn2, mode="a")

            f["x"] = np.arange(10).astype(float)
            g["x"] = np.ones(10).astype(float)

            x1 = f["x"]
            x2 = g["x"]

            assert tokenize(x1) != tokenize(x2)


def test_map_blocks_with_changed_dimension():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(7, 4))

    e = d.map_blocks(lambda b: b.sum(axis=0), chunks=(4,), drop_axis=0, dtype=d.dtype)
    assert e.chunks == ((4, 4),)
    assert_eq(e, x.sum(axis=0))

    # Provided chunks have wrong shape
    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=0), chunks=(), drop_axis=0)

    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=0), chunks=((4, 4, 4),), drop_axis=0)

    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=1), chunks=((3, 4),), drop_axis=1)

    d = da.from_array(x, chunks=(4, 8))
    e = d.map_blocks(lambda b: b.sum(axis=1), drop_axis=1, dtype=d.dtype)
    assert e.chunks == ((4, 3),)
    assert_eq(e, x.sum(axis=1))

    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))
    e = d.map_blocks(
        lambda b: b[None, :, :, None],
        chunks=(1, 4, 4, 1),
        new_axis=[0, 3],
        dtype=d.dtype,
    )
    assert e.chunks == ((1,), (4, 4), (4, 4), (1,))
    assert_eq(e, x[None, :, :, None])

    e = d.map_blocks(lambda b: b[None, :, :, None], new_axis=[0, 3], dtype=d.dtype)
    assert e.chunks == ((1,), (4, 4), (4, 4), (1,))
    assert_eq(e, x[None, :, :, None])

    # Adding axis with a gap
    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b, new_axis=(3, 4))

    # Both new_axis and drop_axis
    d = da.from_array(x, chunks=(8, 4))
    e = d.map_blocks(
        lambda b: b.sum(axis=0)[:, None, None],
        drop_axis=0,
        new_axis=(1, 2),
        dtype=d.dtype,
    )
    assert e.chunks == ((4, 4), (1,), (1,))
    assert_eq(e, x.sum(axis=0)[:, None, None])

    d = da.from_array(x, chunks=(4, 8))
    e = d.map_blocks(
        lambda b: b.sum(axis=1)[:, None, None],
        drop_axis=1,
        new_axis=(1, 2),
        dtype=d.dtype,
    )
    assert e.chunks == ((4, 4), (1,), (1,))
    assert_eq(e, x.sum(axis=1)[:, None, None])


def test_map_blocks_with_negative_drop_axis():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(7, 4))

    for drop_axis in [0, -2]:
        # test with equivalent positive and negative drop_axis
        e = d.map_blocks(
            lambda b: b.sum(axis=0), chunks=(4,), drop_axis=drop_axis, dtype=d.dtype
        )
        assert e.chunks == ((4, 4),)
        assert_eq(e, x.sum(axis=0))


def test_map_blocks_with_invalid_drop_axis():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(7, 4))

    for drop_axis in [x.ndim, -x.ndim - 1]:
        with pytest.raises(ValueError):
            d.map_blocks(
                lambda b: b.sum(axis=0), chunks=(4,), drop_axis=drop_axis, dtype=d.dtype
            )


def test_map_blocks_with_changed_dimension_and_broadcast_chunks():
    # https://github.com/dask/dask/issues/4299
    a = da.from_array([1, 2, 3], 3)
    b = da.from_array(np.array([0, 1, 2, 0, 1, 2]), chunks=3)
    result = da.map_blocks(operator.add, a, b, chunks=b.chunks)
    expected = da.from_array(np.array([1, 3, 5, 1, 3, 5]), chunks=3)
    assert_eq(result, expected)


def test_broadcast_chunks():
    assert broadcast_chunks() == ()

    assert broadcast_chunks(((2, 3),)) == ((2, 3),)

    assert broadcast_chunks(((5, 5),), ((5, 5),)) == ((5, 5),)

    a = ((10, 10, 10), (5, 5))
    b = ((5, 5),)
    assert broadcast_chunks(a, b) == ((10, 10, 10), (5, 5))
    assert broadcast_chunks(b, a) == ((10, 10, 10), (5, 5))

    a = ((10, 10, 10), (5, 5))
    b = ((1,), (5, 5))
    assert broadcast_chunks(a, b) == ((10, 10, 10), (5, 5))

    a = ((10, 10, 10), (5, 5))
    b = ((3, 3), (5, 5))
    with pytest.raises(ValueError):
        broadcast_chunks(a, b)

    a = ((1,), (5, 5))
    b = ((1,), (5, 5))
    assert broadcast_chunks(a, b) == a

    a = ((1,), (np.nan, np.nan, np.nan))
    b = ((3, 3), (1,))
    r = broadcast_chunks(a, b)
    assert r[0] == b[0] and np.allclose(r[1], a[1], equal_nan=True)

    a = ((3, 3), (1,))
    b = ((1,), (np.nan, np.nan, np.nan))
    r = broadcast_chunks(a, b)
    assert r[0] == a[0] and np.allclose(r[1], b[1], equal_nan=True)

    a = ((3, 3), (5, 5))
    b = ((1,), (np.nan, np.nan, np.nan))
    with pytest.raises(ValueError):
        broadcast_chunks(a, b)


def test_chunks_error():
    x = np.ones((10, 10))
    with pytest.raises(ValueError):
        da.from_array(x, chunks=(5,))


def test_array_compute_forward_kwargs():
    x = da.arange(10, chunks=2).sum()
    x.compute(bogus_keyword=10)


def test_dont_fuse_outputs():
    dsk = {("x", 0): np.array([1, 2]), ("x", 1): (inc, ("x", 0))}

    a = da.Array(dsk, "x", chunks=(2,), shape=(4,), dtype=np.array([1]).dtype)
    assert_eq(a, np.array([1, 2, 2, 3], dtype=a.dtype))


def test_dont_dealias_outputs():
    dsk = {
        ("x", 0, 0): np.ones((2, 2)),
        ("x", 0, 1): np.ones((2, 2)),
        ("x", 1, 0): np.ones((2, 2)),
        ("x", 1, 1): ("x", 0, 0),
    }

    a = da.Array(dsk, "x", chunks=(2, 2), shape=(4, 4), dtype=np.ones(1).dtype)
    assert_eq(a, np.ones((4, 4)))


def test_timedelta_op():
    x = np.array([np.timedelta64(10, "h")])
    y = np.timedelta64(1, "h")
    a = da.from_array(x, chunks=(1,)) / y
    assert a.compute() == x / y


def test_to_delayed():
    x = da.random.default_rng().random((4, 4), chunks=(2, 2))
    y = x + 10

    [[a, b], [c, d]] = y.to_delayed()
    assert_eq(a.compute(), y[:2, :2])

    s = 2
    x = da.from_array(np.array(s), chunks=0)
    a = x.to_delayed()[tuple()]
    assert a.compute() == s


def test_to_delayed_optimize_graph():
    x = da.ones((4, 4), chunks=(2, 2))
    y = x[1:][1:][1:][:, 1:][:, 1:][:, 1:]

    # optimizations
    d = y.to_delayed().flatten().tolist()[0]
    assert len([k for k in d.dask if k[0].startswith("getitem")]) == 1
    assert d.key == (y.name, 0, 0)
    assert d.dask.layers.keys() == {"delayed-" + y.name}
    assert d.dask.dependencies == {"delayed-" + y.name: set()}
    assert d.__dask_layers__() == ("delayed-" + y.name,)

    # no optimizations
    d2 = y.to_delayed(optimize_graph=False).flatten().tolist()[0]
    assert d2.dask is y.dask
    assert d2.key == (y.name, 0, 0)
    assert d2.__dask_layers__() == y.__dask_layers__()

    assert (d.compute() == d2.compute()).all()


def test_cumulative():
    rng = np.random.default_rng(0)
    x = da.arange(20, chunks=5)
    assert_eq(x.cumsum(axis=0), np.arange(20).cumsum())
    assert_eq(x.cumprod(axis=0), np.arange(20).cumprod())

    assert_eq(da.nancumsum(x, axis=0), np.nancumsum(np.arange(20)))
    assert_eq(da.nancumprod(x, axis=0), np.nancumprod(np.arange(20)))

    a = rng.random(20)
    a[rng.random(a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=5)
    assert_eq(da.nancumsum(x, axis=0), np.nancumsum(a))
    assert_eq(da.nancumprod(x, axis=0), np.nancumprod(a))

    a = rng.random((20, 24))
    x = da.from_array(a, chunks=(6, 5))
    assert_eq(x.cumsum(axis=0), a.cumsum(axis=0))
    assert_eq(x.cumsum(axis=1), a.cumsum(axis=1))
    assert_eq(x.cumprod(axis=0), a.cumprod(axis=0))
    assert_eq(x.cumprod(axis=1), a.cumprod(axis=1))

    assert_eq(da.nancumsum(x, axis=0), np.nancumsum(a, axis=0))
    assert_eq(da.nancumsum(x, axis=1), np.nancumsum(a, axis=1))
    assert_eq(da.nancumprod(x, axis=0), np.nancumprod(a, axis=0))
    assert_eq(da.nancumprod(x, axis=1), np.nancumprod(a, axis=1))

    a = rng.random((20, 24))
    a[rng.random(a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=(6, 5))
    assert_eq(da.nancumsum(x, axis=0), np.nancumsum(a, axis=0))
    assert_eq(da.nancumsum(x, axis=1), np.nancumsum(a, axis=1))
    assert_eq(da.nancumprod(x, axis=0), np.nancumprod(a, axis=0))
    assert_eq(da.nancumprod(x, axis=1), np.nancumprod(a, axis=1))

    a = rng.random((20, 24, 13))
    x = da.from_array(a, chunks=(6, 5, 4))
    for axis in [0, 1, 2, -1, -2, -3]:
        assert_eq(x.cumsum(axis=axis), a.cumsum(axis=axis))
        assert_eq(x.cumprod(axis=axis), a.cumprod(axis=axis))

        assert_eq(da.nancumsum(x, axis=axis), np.nancumsum(a, axis=axis))
        assert_eq(da.nancumprod(x, axis=axis), np.nancumprod(a, axis=axis))

    a = rng.random((20, 24, 13))
    a[rng.random(a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=(6, 5, 4))
    for axis in [0, 1, 2, -1, -2, -3]:
        assert_eq(da.nancumsum(x, axis=axis), np.nancumsum(a, axis=axis))
        assert_eq(da.nancumprod(x, axis=axis), np.nancumprod(a, axis=axis))

    with pytest.raises(ValueError):
        x.cumsum(axis=3)

    with pytest.raises(ValueError):
        x.cumsum(axis=-4)


def test_from_delayed():
    v = delayed(np.ones)((5, 3))
    x = from_delayed(v, shape=(5, 3), dtype=np.ones(0).dtype)
    assert isinstance(x, Array)
    assert_eq(x, np.ones((5, 3)))


def test_from_delayed_meta():
    v = delayed(np.ones)((5, 3))
    x = from_delayed(v, shape=(5, 3), meta=np.ones(0))
    assert isinstance(x, Array)
    assert isinstance(x._meta, np.ndarray)


def test_A_property():
    x = da.ones(5, chunks=(2,))
    assert x.A is x


def test_copy_mutate():
    x = da.arange(5, chunks=(2,))
    y = x.copy()
    memo = {}
    y2 = copy.deepcopy(x, memo=memo)
    x[x % 2 == 0] = -1

    xx = np.arange(5)
    xx[xx % 2 == 0] = -1
    assert_eq(x, xx)

    assert_eq(y, np.arange(5))
    assert_eq(y2, np.arange(5))
    assert memo[id(x)] is y2


def test_npartitions():
    assert da.ones(5, chunks=(2,)).npartitions == 3
    assert da.ones((5, 5), chunks=(2, 3)).npartitions == 6


def test_elemwise_name():
    assert (da.ones(5, chunks=2) + 1).name.startswith("add-")


def test_map_blocks_name():
    assert da.ones(5, chunks=2).map_blocks(inc).name.startswith("inc-")


def test_map_blocks_token_deprecated():
    with pytest.warns(FutureWarning, match="use `name=` instead"):
        x = da.ones(5, chunks=2).map_blocks(inc, token="foo")
    assert x.name.startswith("foo-")


def test_from_array_names():
    x = np.ones(10)
    a = da.from_array(x, chunks=2)
    assert a.dask.keys() == {(a.name, i) for i in range(5)}


@pytest.mark.parametrize(
    "array", [da.arange(100, chunks=25), da.ones((10, 10), chunks=25)]
)
def test_array_picklable(array):
    from pickle import dumps, loads

    a2 = loads(dumps(array))
    assert_eq(array, a2)

    a3 = da.ma.masked_equal(array, 0)
    assert isinstance(a3._meta, np.ma.MaskedArray)
    a4 = loads(dumps(a3))
    assert_eq(a3, a4)
    assert isinstance(a4._meta, np.ma.MaskedArray)


def test_from_array_raises_on_bad_chunks():
    x = np.ones(10)

    with pytest.raises(ValueError):
        da.from_array(x, chunks=(5, 5, 5))

    # with pytest.raises(ValueError):
    #      da.from_array(x, chunks=100)

    with pytest.raises(ValueError):
        da.from_array(x, chunks=((5, 5, 5),))


def test_concatenate_axes():
    x = np.ones((2, 2, 2))

    assert_eq(concatenate_axes([x, x], axes=[0]), np.ones((4, 2, 2)))
    assert_eq(concatenate_axes([x, x, x], axes=[0]), np.ones((6, 2, 2)))
    assert_eq(concatenate_axes([x, x], axes=[1]), np.ones((2, 4, 2)))
    assert_eq(concatenate_axes([[x, x], [x, x]], axes=[0, 1]), np.ones((4, 4, 2)))
    assert_eq(concatenate_axes([[x, x], [x, x]], axes=[0, 2]), np.ones((4, 2, 4)))
    assert_eq(concatenate_axes([[x, x, x], [x, x, x]], axes=[1, 2]), np.ones((2, 4, 6)))

    with pytest.raises(ValueError):
        concatenate_axes(
            [[x, x], [x, x]], axes=[0]
        )  # not all nested lists accounted for
    with pytest.raises(ValueError):
        concatenate_axes([x, x], axes=[0, 1, 2, 3])  # too many axes


def test_blockwise_concatenate():
    x = da.ones((4, 4, 4), chunks=(2, 2, 2))
    y = da.ones((4, 4), chunks=(2, 2))

    def f(a, b):
        assert isinstance(a, np.ndarray)
        assert isinstance(b, np.ndarray)

        assert a.shape == (2, 4, 4)
        assert b.shape == (4, 4)

        return (a + b).sum(axis=(1, 2))

    z = da.blockwise(f, "i", x, "ijk", y, "jk", concatenate=True, dtype=x.dtype)
    assert_eq(z, np.ones(4) * 32)

    z = da.blockwise(add, "ij", y, "ij", y, "ij", concatenate=True, dtype=x.dtype)
    assert_eq(z, np.ones((4, 4)) * 2)

    def f(a, b, c):
        assert isinstance(a, np.ndarray)
        assert isinstance(b, np.ndarray)
        assert isinstance(c, np.ndarray)

        assert a.shape == (4, 2, 4)
        assert b.shape == (4, 4)
        assert c.shape == (4, 2)

        return np.ones(2)

    z = da.blockwise(
        f, "j", x, "ijk", y, "ki", y, "ij", concatenate=True, dtype=x.dtype
    )
    assert_eq(z, np.ones(4), check_shape=False)


def test_common_blockdim():
    assert common_blockdim([(5,), (5,)]) == (5,)
    assert common_blockdim([(5,), (2, 3)]) == (2, 3)
    assert common_blockdim([(5, 5), (2, 3, 5)]) == (2, 3, 5)
    assert common_blockdim([(5, 5), (2, 3, 5)]) == (2, 3, 5)
    assert common_blockdim([(5, 2, 3), (2, 3, 5)]) == (2, 3, 2, 3)

    assert common_blockdim([(1, 2), (2, 1)]) == (1, 1, 1)
    assert common_blockdim([(1, 2, 2), (2, 1, 2), (2, 2, 1)]) == (1, 1, 1, 1, 1)


def test_uneven_chunks_that_fit_neatly():
    x = da.arange(10, chunks=((5, 5),))
    y = da.ones(10, chunks=((5, 2, 3),))

    assert_eq(x + y, np.arange(10) + np.ones(10))

    z = x + y
    assert z.chunks == ((5, 2, 3),)


def test_elemwise_uneven_chunks():
    rng = da.random.default_rng()
    x = da.arange(10, chunks=((4, 6),))
    y = da.ones(10, chunks=((6, 4),))

    assert_eq(x + y, np.arange(10) + np.ones(10))

    z = x + y
    assert z.chunks == ((4, 2, 4),)

    x = rng.random((10, 10), chunks=((4, 6), (5, 2, 3)))
    y = rng.random((4, 10, 10), chunks=((2, 2), (6, 4), (2, 3, 5)))

    z = x + y
    assert_eq(x + y, x.compute() + y.compute())
    assert z.chunks == ((2, 2), (4, 2, 4), (2, 3, 2, 3))


def test_uneven_chunks_blockwise():
    rng = da.random.default_rng()
    x = rng.random((10, 10), chunks=((2, 3, 2, 3), (5, 5)))
    y = rng.random((10, 10), chunks=((4, 4, 2), (4, 2, 4)))
    z = da.blockwise(np.dot, "ik", x, "ij", y, "jk", dtype=x.dtype, concatenate=True)
    assert z.chunks == (x.chunks[0], y.chunks[1])

    assert_eq(z, x.compute().dot(y))


def test_warn_bad_rechunking():
    x = da.ones((20, 20), chunks=(20, 1))
    y = da.ones((20, 20), chunks=(1, 20))

    with pytest.warns(da.core.PerformanceWarning, match="factor of 20"):
        x + y


def test_concatenate_stack_dont_warn():
    with warnings.catch_warnings(record=True) as record:
        da.concatenate([da.ones(2, chunks=1)] * 62)
    assert not record

    with warnings.catch_warnings(record=True) as record:
        da.stack([da.ones(2, chunks=1)] * 62)
    assert not record


def test_map_blocks_delayed():
    x = da.ones((10, 10), chunks=(5, 5))
    y = np.ones((5, 5))

    z = x.map_blocks(add, y, dtype=x.dtype)
    z.dask.validate()
    dask.optimize(z)[0].dask.validate()

    yy = delayed(y)
    zz = x.map_blocks(add, yy, dtype=x.dtype)
    zz.dask.validate()
    dask.optimize(zz)[0].dask.validate()

    assert_eq(z, zz)

    assert yy.key in zz.dask


def test_no_chunks():
    X = np.arange(11)
    dsk = {("x", 0): np.arange(5), ("x", 1): np.arange(5, 11)}
    x = Array(dsk, "x", ((np.nan, np.nan),), np.arange(1).dtype)
    assert_eq(x + 1, X + 1)
    assert_eq(x.sum(), X.sum())
    assert_eq((x + 1).std(), (X + 1).std())
    assert_eq((x + x).std(), (X + X).std())
    assert_eq((x + x).std(keepdims=True), (X + X).std(keepdims=True))


def test_no_chunks_2d():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((np.nan, np.nan), (np.nan, np.nan, np.nan))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # divide by zero
        assert_eq(da.log(x), np.log(X))
    assert_eq(x.T, X.T)
    assert_eq(x.sum(axis=0, keepdims=True), X.sum(axis=0, keepdims=True))
    assert_eq(x.sum(axis=1, keepdims=True), X.sum(axis=1, keepdims=True))
    assert_eq(x.dot(x.T + 1), X.dot(X.T + 1))


def test_no_chunks_yes_chunks():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    assert (x + 1).chunks == ((2, 2), (np.nan, np.nan, np.nan))
    assert (x.T).chunks == ((np.nan, np.nan, np.nan), (2, 2))
    assert (x.dot(x.T)).chunks == ((2, 2), (2, 2))


def test_raise_informative_errors_no_chunks():
    X = np.arange(10)
    a = da.from_array(X, chunks=(5, 5))
    a._chunks = ((np.nan, np.nan),)

    b = da.from_array(X, chunks=(4, 4, 2))
    b._chunks = ((np.nan, np.nan, np.nan),)

    for op in [
        lambda: a + b,
        lambda: a[1],
        lambda: a[::2],
        lambda: a[-5],
        lambda: a.rechunk(3),
        lambda: a.reshape(2, 5),
    ]:
        with pytest.raises(ValueError) as e:
            op()
        if "chunk" not in str(e.value) or "unknown" not in str(e.value):
            op()


def test_no_chunks_slicing_2d():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    assert_eq(x[0], X[0])

    for op in [lambda: x[:, 4], lambda: x[:, ::2], lambda: x[0, 2:4]]:
        with pytest.raises(ValueError, match="chunk sizes are unknown"):
            op()


def test_index_array_with_array_1d():
    x = np.arange(10)
    dx = da.from_array(x, chunks=(5,))
    dx._chunks = ((np.nan, np.nan),)

    assert_eq(x[x > 6], dx[dx > 6])
    assert_eq(x[x % 2 == 0], dx[dx % 2 == 0])

    dy = da.ones(11, chunks=(3,))

    with pytest.raises(ValueError):
        dx[dy > 5]


def test_index_array_with_array_2d():
    x = np.arange(24).reshape((4, 6))
    dx = da.from_array(x, chunks=(2, 2))

    assert_eq(x[x > 6], dx[dx > 6])
    assert_eq(x[x % 2 == 0], dx[dx % 2 == 0])

    # Test with unknown chunks
    dx._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    with pytest.warns(UserWarning, match="different ordering") as record:
        assert sorted(x[x % 2 == 0].tolist()) == sorted(
            dx[dx % 2 == 0].compute().tolist()
        )
        assert sorted(x[x > 6].tolist()) == sorted(dx[dx > 6].compute().tolist())

    assert len(record) == 2


@pytest.mark.xfail(reason="Chunking does not align well")
def test_index_array_with_array_3d_2d():
    x = np.arange(4**3).reshape((4, 4, 4))
    dx = da.from_array(x, chunks=(2, 2, 2))

    ind = np.random.default_rng().random((4, 4)) > 0.5
    ind = np.arange(4**2).reshape((4, 4)) % 2 == 0
    dind = da.from_array(ind, (2, 2))

    assert_eq(x[ind], dx[dind])
    assert_eq(x[:, ind], dx[:, dind])


def test_setitem_1d():
    x = np.arange(10)
    dx = da.from_array(x.copy(), chunks=(5,))

    x[x > 6] = -1
    x[x % 2 == 0] = -2
    x[[2, 3]] = -3

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2
    dx[da.asarray([2, 3])] = -3

    assert_eq(x, dx)


def test_setitem_masked():
    # Test np.ma.masked assignment to object-type arrays
    x = np.ma.array(["a", 1, 3.14], dtype=object)
    dx = da.from_array(x.copy(), chunks=2)

    x[...] = np.ma.masked
    dx[...] = np.ma.masked

    assert_eq(x.mask, da.ma.getmaskarray(dx))


def test_setitem_hardmask():
    x = np.ma.array([1, 2, 3, 4], dtype=int)
    x.harden_mask()

    y = x.copy()
    assert y.hardmask

    x[0] = np.ma.masked
    x[0:2] = np.ma.masked

    dx = da.from_array(y)
    dx[0] = np.ma.masked
    dx[0:2] = np.ma.masked
    assert_eq(x, dx)


def test_setitem_slice_twice():
    x = np.array([1, 2, 3, 4, 5, 6], dtype=int)
    val = np.array([0, 0], dtype=int)
    y = x.copy()

    x[0:2] = val
    x[4:6] = val

    dx = da.from_array(y)
    dx[0:2] = val
    dx[4:6] = val
    assert_eq(x, dx)


def test_setitem_2d():
    x = np.arange(24).reshape((4, 6))
    dx = da.from_array(x.copy(), chunks=(2, 2))

    x[x > 6] = -1
    x[x % 2 == 0] = -2

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2

    assert_eq(x, dx)


def test_setitem_extended_API_0d():
    # 0-d array
    x = np.array(9)
    dx = da.from_array(9)

    x[()] = -1
    dx[()] = -1
    assert_eq(x, dx.compute())

    x[...] = -11
    dx[...] = -11
    assert_eq(x, dx.compute())


@pytest.mark.parametrize(
    "index, value",
    [
        [Ellipsis, -1],
        [slice(2, 8, 2), -2],
        [slice(8, None, 2), -3],
        [slice(8, None, 2), [-30]],
        [slice(1, None, -2), -4],
        [slice(1, None, -2), [-40]],
        [slice(3, None, 2), -5],
        [slice(-3, None, -2), -6],
        [slice(1, None, -2), -4],
        [slice(3, None, 2), -5],
        [slice(3, None, 2), [10, 11, 12, 13]],
        [slice(-4, None, -2), [14, 15, 16, 17]],
    ],
)
def test_setitem_extended_API_1d(index, value):
    # 1-d array
    x = np.arange(10)
    dx = da.from_array(x, chunks=(4, 6))
    dx[index] = value
    x[index] = value
    assert_eq(x, dx.compute())


@pytest.mark.parametrize(
    "index, value",
    [
        [Ellipsis, -1],
        [(slice(None, None, 2), slice(None, None, -1)), -1],
        [slice(1, None, 2), -1],
        [[4, 3, 1], -1],
        [(Ellipsis, 4), -1],
        [5, -1],
        [(slice(None), 2), range(6)],
        [3, range(10)],
        [(slice(None), [3, 5, 6]), [-30, -31, -32]],
        [([-1, 0, 1], 2), [-30, -31, -32]],
        [(slice(None, 2), slice(None, 3)), [-50, -51, -52]],
        [(slice(None), [6, 1, 3]), [-60, -61, -62]],
        [(slice(1, 3), slice(1, 4)), [[-70, -71, -72]]],
        [(slice(None), [9, 8, 8]), [-80, -81, 91]],
        [([True, False, False, False, True, False], 2), -1],
        [(3, [True, True, False, True, True, False, True, False, True, True]), -1],
        [(np.array([False, False, True, True, False, False]), slice(5, 7)), -1],
        [
            (
                4,
                da.from_array(
                    [False, False, True, True, False, False, True, False, False, True]
                ),
            ),
            -1,
        ],
        [
            (
                slice(2, 4),
                da.from_array(
                    [False, False, True, True, False, False, True, False, False, True]
                ),
            ),
            [[-100, -101, -102, -103], [-200, -201, -202, -203]],
        ],
        [slice(5, None, 2), -99],
        [slice(5, None, 2), range(1, 11)],
        [slice(1, None, -2), -98],
        [slice(1, None, -2), range(11, 21)],
    ],
)
def test_setitem_extended_API_2d(index, value):
    # 2-d array
    x = np.ma.arange(60).reshape((6, 10))
    dx = da.from_array(x, chunks=(2, 3))
    dx[index] = value
    x[index] = value
    assert_eq(x, dx.compute())


def test_setitem_extended_API_2d_rhs_func_of_lhs():
    # Cases:
    # * RHS and/or indices are a function of the LHS
    # * Indices have unknown chunk sizes
    # * RHS has extra leading size 1 dimensions compared to LHS
    x = np.arange(60).reshape((6, 10))
    chunks = (2, 3)

    dx = da.from_array(x, chunks=chunks)
    dx[2:4, dx[0] > 3] = -5
    x[2:4, x[0] > 3] = -5
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[2, dx[0] < -2] = -7
    x[2, x[0] < -2] = -7
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[dx % 2 == 0] = -8
    x[x % 2 == 0] = -8
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[dx % 2 == 0] = -8
    x[x % 2 == 0] = -8
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[3:5, 5:1:-2] = -dx[:2, 4:1:-2]
    x[3:5, 5:1:-2] = -x[:2, 4:1:-2]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[0, 1:3] = -dx[0, 4:2:-1]
    x[0, 1:3] = -x[0, 4:2:-1]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[...] = dx
    x[...] = x
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[...] = dx[...]
    x[...] = x[...]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[0] = dx[-1]
    x[0] = x[-1]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[0, :] = dx[-2, :]
    x[0, :] = x[-2, :]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[:, 1] = dx[:, -3]
    x[:, 1] = x[:, -3]
    assert_eq(x, dx.compute())

    index = da.from_array([0, 2], chunks=(2,))
    dx = da.from_array(x, chunks=chunks)
    dx[index, 8] = [99, 88]
    x[[0, 2], 8] = [99, 88]
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=chunks)
    dx[:, index] = dx[:, :2]
    x[:, [0, 2]] = x[:, :2]
    assert_eq(x, dx.compute())

    index = da.where(da.arange(3, chunks=(1,)) < 2)[0]
    dx = da.from_array(x, chunks=chunks)
    dx[index, 7] = [-23, -33]
    x[index.compute(), 7] = [-23, -33]
    assert_eq(x, dx.compute())

    index = da.where(da.arange(3, chunks=(1,)) < 2)[0]
    dx = da.from_array(x, chunks=chunks)
    dx[(index,)] = -34
    x[(index.compute(),)] = -34
    assert_eq(x, dx.compute())

    index = index - 4
    dx = da.from_array(x, chunks=chunks)
    dx[index, 7] = [-43, -53]
    x[index.compute(), 7] = [-43, -53]
    assert_eq(x, dx.compute())

    index = da.from_array([0, -1], chunks=(1,))
    x[[0, -1]] = 9999
    dx[(index,)] = 9999
    assert_eq(x, dx.compute())

    dx = da.from_array(x, chunks=(-1, -1))
    dx[...] = da.from_array(x, chunks=chunks)
    assert_eq(x, dx.compute())

    # RHS has extra leading size 1 dimensions compared to LHS
    dx = da.from_array(x.copy(), chunks=(2, 3))
    v = x.reshape((1, 1) + x.shape)
    x[...] = v
    dx[...] = v
    assert_eq(x, dx.compute())

    index = da.where(da.arange(3, chunks=(1,)) < 2)[0]
    v = -np.arange(12).reshape(1, 1, 6, 2)
    x[:, [0, 1]] = v
    dx[:, index] = v
    assert_eq(x, dx.compute())


@pytest.mark.parametrize(
    "index, value",
    [
        [(1, slice(1, 7, 2)), np.ma.masked],
        [(slice(1, 5, 2), [7, 5]), np.ma.masked_all((2, 2))],
    ],
)
def test_setitem_extended_API_2d_mask(index, value):
    x = np.ma.arange(60).reshape((6, 10))
    dx = da.from_array(x.data, chunks=(2, 3))
    # See https://github.com/numpy/numpy/issues/23000 for the `RuntimeWarning`
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=RuntimeWarning,
            message="invalid value encountered in cast",
        )
        x[index] = value
        dx[index] = value
    dx = dx.persist()
    assert_eq(x, dx.compute())
    assert_eq(x.mask, da.ma.getmaskarray(dx).compute())


def test_setitem_on_read_only_blocks():
    # Outputs of broadcast_trick-style functions contain read-only
    # arrays
    dx = da.empty((4, 6), dtype=float, chunks=(2, 2))
    dx[0] = 99

    assert_eq(dx[0, 0], 99.0)

    dx[0:2] = 88

    assert_eq(dx[0, 0], 88.0)


def test_setitem_errs():
    x = da.ones((4, 4), chunks=(2, 2))

    with pytest.raises(ValueError):
        x[x > 1] = x

    # Shape mismatch
    with pytest.raises(ValueError):
        x[[True, True, False, False], 0] = [2, 3, 4]

    with pytest.raises(ValueError):
        x[[True, True, True, False], 0] = [2, 3]

    with pytest.raises(ValueError):
        x[0, [True, True, True, False]] = [2, 3]

    with pytest.raises(ValueError):
        x[0, [True, True, True, False]] = [1, 2, 3, 4, 5]

    with pytest.raises(ValueError):
        x[da.from_array([True, True, True, False]), 0] = [1, 2, 3, 4, 5]

    with pytest.raises(ValueError):
        x[0, da.from_array([True, False, False, True])] = [1, 2, 3, 4, 5]

    with pytest.raises(ValueError):
        x[:, 0] = [2, 3, 4]

    with pytest.raises(ValueError):
        x[0, :] = [1, 2, 3, 4, 5]

    x = da.ones((4, 4), chunks=(2, 2))

    # Too many indices
    with pytest.raises(IndexError):
        x[:, :, :] = 2

    # 2-d boolean indexing a single dimension
    with pytest.raises(IndexError):
        x[[[True, True, False, False]], 0] = 5

    # 2-d indexing a single dimension
    with pytest.raises(IndexError):
        x[[[1, 2, 3]], 0] = 5

    # Multiple 1-d boolean/integer arrays
    with pytest.raises(NotImplementedError):
        x[[1, 2], [2, 3]] = 6

    with pytest.raises(NotImplementedError):
        x[[True, True, False, False], [2, 3]] = 5

    with pytest.raises(NotImplementedError):
        x[[True, True, False, False], [False, True, False, False]] = 7

    # scalar boolean indexing
    with pytest.raises(NotImplementedError):
        x[True] = 5

    with pytest.raises(NotImplementedError):
        x[np.array(True)] = 5

    with pytest.raises(NotImplementedError):
        x[0, da.from_array(True)] = 5

    # Scalar arrays
    y = da.from_array(np.array(1))
    with pytest.raises(IndexError):
        y[:] = 2

    # RHS has non-brodacastable extra leading dimensions
    x = np.arange(12).reshape((3, 4))
    dx = da.from_array(x, chunks=(2, 2))
    with pytest.raises(ValueError):
        dx[...] = np.arange(24).reshape((2, 1, 3, 4))

    # RHS doesn't have chunks set
    dx = da.unique(da.random.default_rng().random([10]))
    with pytest.raises(ValueError, match="Arrays chunk sizes are unknown"):
        dx[0] = 0

    # np.nan assigned to integer array
    x = da.ones((3, 3), dtype=int)
    with pytest.raises(ValueError, match="cannot convert float NaN to integer"):
        x[:, 1] = np.nan
    with pytest.raises(ValueError, match="cannot convert float infinity to integer"):
        x[:, 1] = np.inf
    with pytest.raises(ValueError, match="cannot convert float infinity to integer"):
        x[:, 1] = -np.inf


@pytest.mark.parametrize("idx_namespace", [np, da])
def test_setitem_bool_index_errs(idx_namespace):
    x = da.ones((3, 4), chunks=(2, 2))
    y = da.ones(4, chunks=2)
    array = idx_namespace.array

    # Shape mismatch
    with pytest.raises(ValueError):
        y[array([True, True, True, False])] = [2, 3]

    with pytest.raises(ValueError):
        # A naive where(idx, val, x) would produce a result
        y[array([True, True, True, False])] = [1, 2, 3, 4]

    with pytest.raises(ValueError):
        y[array([True, True, True, False])] = [1, 2, 3, 4, 5]

    with pytest.raises(ValueError):
        y[array([True, False, False, True])] = [1, 2, 3, 4, 5]

    # Too many/not enough booleans
    with pytest.raises(IndexError):
        y[array([True, False, True])] = 5

    with pytest.raises(IndexError):
        y[array([False, True, True, True, False])] = 5

    # Situations where a naive da.where(idx, val, x) would produce a result
    with pytest.raises(IndexError):
        x[array([True, False, False, True])] = 1

    with pytest.raises(IndexError):
        y[array([[True], [False]])] = 1  # da.where would broadcast to ndim=2


def test_zero_slice_dtypes():
    x = da.arange(5, chunks=1)
    y = x[[]]
    assert y.dtype == x.dtype
    assert y.shape == (0,)
    assert_eq(x[[]], np.arange(5)[[]])


def test_zero_sized_array_rechunk():
    x = da.arange(5, chunks=1)[:0]
    y = da.blockwise(identity, "i", x, "i", dtype=x.dtype)
    assert_eq(x, y)


def test_blockwise_zero_shape():
    da.blockwise(
        lambda x: x,
        "i",
        da.arange(10, chunks=10),
        "i",
        da.from_array(np.ones((0, 2)), ((0,), 2)),
        "ab",
        da.from_array(np.ones((0,)), ((0,),)),
        "a",
        dtype="float64",
    )


def test_blockwise_zero_shape_new_axes():
    da.blockwise(
        lambda x: np.ones(42),
        "i",
        da.from_array(np.ones((0, 2)), ((0,), 2)),
        "ab",
        da.from_array(np.ones((0,)), ((0,),)),
        "a",
        dtype="float64",
        new_axes={"i": 42},
    )


def test_broadcast_against_zero_shape():
    assert_eq(da.arange(1, chunks=1)[:0] + 0, np.arange(1)[:0] + 0)
    assert_eq(da.arange(1, chunks=1)[:0] + 0.1, np.arange(1)[:0] + 0.1)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:0] + 0, np.ones((5, 5))[:0] + 0)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:0] + 0.1, np.ones((5, 5))[:0] + 0.1)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:, :0] + 0, np.ones((5, 5))[:, :0] + 0)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:, :0] + 0.1, np.ones((5, 5))[:, :0] + 0.1)


def test_from_array_name():
    x = np.array([1, 2, 3, 4, 5])
    chunks = x.shape
    # Default is tokenize the array
    dx = da.from_array(x, chunks=chunks)
    hashed_name = dx.name
    assert da.from_array(x, chunks=chunks).name == hashed_name
    # Specify name directly
    assert da.from_array(x, chunks=chunks, name="x").name == "x"
    # False gives a random name
    dx2 = da.from_array(x, chunks=chunks, name=False)
    dx3 = da.from_array(x, chunks=chunks, name=False)
    assert dx2.name != hashed_name
    assert dx3.name != hashed_name
    assert dx2.name != dx3.name


def test_concatenate_errs():
    with pytest.raises(ValueError, match=r"Shapes.*\(2, 1\)"):
        da.concatenate(
            [da.zeros((2, 1), chunks=(2, 1)), da.zeros((2, 3), chunks=(2, 3))]
        )

    with pytest.raises(ValueError):
        da.concatenate(
            [da.zeros((1, 2), chunks=(1, 2)), da.zeros((3, 2), chunks=(3, 2))], axis=1
        )


def test_stack_errs():
    with pytest.raises(ValueError) as e:
        da.stack([da.zeros((2,), chunks=2)] * 10 + [da.zeros((3,), chunks=3)] * 10)

    assert (
        str(e.value)
        == "Stacked arrays must have the same shape. The first array had shape (2,), while array 11 has shape (3,)."
    )
    assert len(str(e.value)) < 105


def test_blockwise_with_numpy_arrays():
    x = np.ones(10)
    y = da.ones(10, chunks=(5,))

    assert_eq(x + y, x + x)

    s = da.sum(x)
    assert any(isinstance(v, np.ndarray) for v in s.dask.values())


@pytest.mark.parametrize("chunks", (100, 6))
@pytest.mark.parametrize("other", [[0, 0, 1], [2, 1, 3], (0, 0, 1)])
def test_elemwise_with_lists(chunks, other):
    x = np.arange(12).reshape((4, 3))
    d = da.arange(12, chunks=chunks).reshape((4, 3))

    x2 = np.vstack([x[:, 0], x[:, 1], x[:, 2]]).T
    d2 = da.vstack([d[:, 0], d[:, 1], d[:, 2]]).T

    assert_eq(x2, d2)

    x3 = x2 * other
    d3 = d2 * other

    assert_eq(x3, d3)


def test_constructor_plugin():
    L = []
    L2 = []
    with dask.config.set(array_plugins=[L.append, L2.append]):
        x = da.ones(10, chunks=5)
        y = x + 1

    assert L == L2 == [x, y]

    with dask.config.set(array_plugins=[lambda x: x.compute()]):
        x = da.ones(10, chunks=5)
        y = x + 1

    assert isinstance(y, np.ndarray)
    assert len(L) == 2


def test_no_warnings_on_metadata():
    x = da.ones(5, chunks=3)
    with warnings.catch_warnings(record=True) as record:
        da.arccos(x)

    assert not record


def test_delayed_array_key_hygeine():
    a = da.zeros((1,), chunks=(1,))
    d = delayed(identity)(a)
    b = da.from_delayed(d, shape=a.shape, dtype=a.dtype)
    assert_eq(a, b)


def test_empty_chunks_in_array_len():
    x = da.ones((), chunks=())
    with pytest.raises(TypeError) as exc_info:
        len(x)

    err_msg = "len() of unsized object"
    assert err_msg in str(exc_info.value)


@pytest.mark.parametrize("dtype", [None, [("a", "f4"), ("b", object)]])
def test_meta(dtype):
    a = da.zeros((1,), chunks=(1,))
    assert a._meta.dtype == a.dtype
    assert isinstance(a._meta, np.ndarray)
    assert a.nbytes < 1000


@pytest.mark.parametrize(
    "shape,limit,expected",
    [
        (100, 10, (10,) * 10),
        (20, 10, (10, 10)),
        (20, 5, (5, 5, 5, 5)),
        (24, 5, (5, 5, 5, 5, 4)),
        (23, 5, (5, 5, 5, 5, 3)),  # relatively prime, don't use 1s
        (1000, 167, (167, 167, 167, 167, 167, 165)),
    ],
)
def test_normalize_chunks_auto_1d(shape, limit, expected):
    result = normalize_chunks("auto", (shape,), limit=limit, dtype=np.uint8)
    assert result == (expected,)


@pytest.mark.parametrize(
    "shape,chunks,limit,expected",
    [
        ((20, 20), ("auto", 2), 20, ((10, 10), (2,) * 10)),
        (
            (20, 20),
            ("auto", (2, 2, 2, 2, 2, 5, 5)),
            20,
            ((4, 4, 4, 4, 4), (2, 2, 2, 2, 2, 5, 5)),
        ),
        ((1, 20), "auto", 10, ((1,), (10, 10))),
    ],
)
def test_normalize_chunks_auto_2d(shape, chunks, limit, expected):
    result = normalize_chunks(chunks, shape, limit=limit, dtype="uint8")
    assert result == expected


def test_normalize_chunks_auto_3d():
    result = normalize_chunks(
        ("auto", "auto", 2), (20, 20, 20), limit=200, dtype="uint8"
    )
    expected = ((10, 10), (10, 10), (2,) * 10)
    assert result == expected

    result = normalize_chunks("auto", (20, 20, 20), limit=8, dtype="uint8")
    expected = ((2,) * 10,) * 3
    assert result == expected


def test_constructors_chunks_dict():
    x = da.ones((20, 20), chunks={0: 10, 1: 5})
    assert x.chunks == ((10, 10), (5, 5, 5, 5))

    x = da.ones((20, 20), chunks={0: 10, 1: "auto"})
    assert x.chunks == ((10, 10), (20,))


def test_from_array_chunks_dict():
    with dask.config.set({"array.chunk-size": "128kiB"}):
        x = np.empty((100, 100, 100))
        y = da.from_array(x, chunks={0: 10, 1: -1, 2: "auto"})
        z = da.from_array(x, chunks=(10, 100, (16,) * 6 + (4,)))
        assert y.chunks == z.chunks


@pytest.mark.parametrize("dtype", [object, [("a", object), ("b", int)]])
def test_normalize_chunks_object_dtype(dtype):
    x = np.array(["a", "abc"], dtype=object)
    with pytest.raises(NotImplementedError):
        da.from_array(x, chunks="auto")


def test_normalize_chunks_tuples_of_tuples():
    result = normalize_chunks(((2, 3, 5), "auto"), (10, 10), limit=10, dtype=np.uint8)
    expected = ((2, 3, 5), (2, 2, 2, 2, 2))
    assert result == expected


def test_normalize_chunks_nan():
    with pytest.raises(ValueError) as info:
        normalize_chunks("auto", (np.nan,), limit=10, dtype=np.uint8)
    assert "auto" in str(info.value)
    with pytest.raises(ValueError) as info:
        normalize_chunks(((np.nan, np.nan), "auto"), (10, 10), limit=10, dtype=np.uint8)
    assert "auto" in str(info.value)


def test_pandas_from_dask_array():
    pd = pytest.importorskip("pandas")
    a = da.ones((12,), chunks=4)
    s = pd.Series(a, index=range(12))
    assert s.dtype == a.dtype
    assert_eq(s.values, a)


def test_from_zarr_unique_name():
    zarr = pytest.importorskip("zarr")
    a = zarr.array([1, 2, 3])
    b = zarr.array([4, 5, 6])

    assert da.from_zarr(a).name != da.from_zarr(b).name


def test_from_zarr_name():
    zarr = pytest.importorskip("zarr")
    a = zarr.array([1, 2, 3])
    assert da.from_zarr(a, name="foo").name == "foo"


def test_zarr_roundtrip():
    pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.zeros((3, 3), chunks=(1, 1))
        a.to_zarr(d)
        a2 = da.from_zarr(d)
        assert_eq(a, a2)
        assert a2.chunks == a.chunks


def test_zarr_roundtrip_with_path_like():
    pytest.importorskip("zarr")
    with tmpdir() as d:
        path = pathlib.Path(d)
        a = da.zeros((3, 3), chunks=(1, 1))
        a.to_zarr(path)
        a2 = da.from_zarr(path)
        assert_eq(a, a2)
        assert a2.chunks == a.chunks


def test_to_zarr_accepts_empty_array_without_exception_raised():
    pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.from_array(np.arange(0))
        a.to_zarr(d)


@pytest.mark.parametrize("compute", [False, True])
def test_zarr_return_stored(compute):
    pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.zeros((3, 3), chunks=(1, 1))
        a2 = a.to_zarr(d, compute=compute, return_stored=True)
        assert isinstance(a2, Array)
        assert_eq(a, a2, check_graph=False)
        assert a2.chunks == a.chunks


@pytest.mark.parametrize("inline_array", [True, False])
def test_zarr_inline_array(inline_array):
    zarr = pytest.importorskip("zarr")
    a = zarr.array([1, 2, 3])
    dsk = dict(da.from_zarr(a, inline_array=inline_array).dask)
    assert len(dsk) == (0 if inline_array else 1) + 1
    assert (a in dsk.values()) is not inline_array


def test_zarr_existing_array():
    zarr = pytest.importorskip("zarr")
    c = (1, 1)
    a = da.ones((3, 3), chunks=c)
    z = zarr.zeros_like(a, chunks=c)
    a.to_zarr(z)
    a2 = da.from_zarr(z)
    assert_eq(a, a2)
    assert a2.chunks == a.chunks


def test_to_zarr_unknown_chunks_raises():
    pytest.importorskip("zarr")
    a = da.random.default_rng().random((10,), chunks=(3,))
    a = a[a > 0.5]
    with pytest.raises(ValueError, match="unknown chunk sizes"):
        a.to_zarr({})


def test_read_zarr_chunks():
    pytest.importorskip("zarr")
    a = da.zeros((9,), chunks=(3,))
    with tmpdir() as d:
        a.to_zarr(d)
        arr = da.from_zarr(d, chunks=(5,))
        assert arr.chunks == ((5, 4),)


def test_zarr_pass_store():
    zarr = pytest.importorskip("zarr")

    with tmpdir() as d:
        if Version(zarr.__version__) < Version("3.0.0.a0"):
            store = zarr.storage.DirectoryStore(d)
        else:
            store = zarr.storage.LocalStore(d, read_only=False)
        a = da.zeros((3, 3), chunks=(1, 1))
        a.to_zarr(store)
        a2 = da.from_zarr(store)
        assert_eq(a, a2)
        assert a2.chunks == a.chunks


def test_zarr_group():
    zarr = pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.zeros((3, 3), chunks=(1, 1))
        a.to_zarr(d, component="test")
        with pytest.raises((OSError, ValueError)):
            a.to_zarr(d, component="test", overwrite=False)
        a.to_zarr(d, component="test", overwrite=True)

        # second time is fine, group exists
        a.to_zarr(d, component="test2", overwrite=False)
        a.to_zarr(d, component="nested/test", overwrite=False)

        group = zarr.open_group(store=d, mode="r")
        assert set(group) == {"nested", "test", "test2"}
        assert "test" in group["nested"]

        a2 = da.from_zarr(d, component="test")
        assert_eq(a, a2)
        assert a2.chunks == a.chunks


@pytest.mark.parametrize(
    "shape, chunks, expect_rechunk",
    [
        ((6, 2), ((2, 1, 1, 2), 1), True),
        ((6, 2), ((2, 1, 2, 1), 1), True),
        ((7, 2), ((2, 2, 2, 1), 1), False),
        ((2, 7), (1, (2, 2, 2, 1)), False),
        ((2, 6), (1, (2, 1, 2, 1)), True),
    ],
)
def test_zarr_irregular_chunks(shape, chunks, expect_rechunk):
    pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.zeros(shape, chunks=chunks)  # ((2, 1, 1, 2), 1))
        store_delayed = a.to_zarr(d, component="test", compute=False)
        assert (
            any("rechunk" in key_split(k) for k in dict(store_delayed.dask))
            is expect_rechunk
        )
        store_delayed.compute()


@pytest.mark.parametrize(
    "data",
    [
        [(), True],
        [((1,),), True],
        [((1, 1, 1),), True],
        [((1,), (1,)), True],
        [((2, 2, 1),), True],
        [((2, 2, 3),), False],
        [((1, 1, 1), (2, 2, 3)), False],
        [((1, 2, 1),), False],
    ],
)
def test_regular_chunks(data):
    from dask.array.core import _check_regular_chunks

    chunkset, expected = data
    assert _check_regular_chunks(chunkset) == expected


def test_zarr_nocompute():
    pytest.importorskip("zarr")
    with tmpdir() as d:
        a = da.zeros((3, 3), chunks=(1, 1))
        out = a.to_zarr(d, compute=False)
        assert isinstance(out, Delayed)
        dask.compute(out)
        a2 = da.from_zarr(d)
        assert_eq(a, a2)
        assert a2.chunks == a.chunks


def test_zarr_regions():
    zarr = pytest.importorskip("zarr")

    a = da.arange(16).reshape((4, 4)).rechunk(2)
    z = zarr.zeros_like(a, chunks=2)

    a[:2, :2].to_zarr(z, region=(slice(2), slice(2)))
    a2 = da.from_zarr(z)
    expected = [[0, 1, 0, 0], [4, 5, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    assert_eq(a2, expected)
    assert a2.chunks == a.chunks

    a[:3, 3:4].to_zarr(z, region=(slice(1, 4), slice(2, 3)))
    a2 = da.from_zarr(z)
    expected = [[0, 1, 0, 0], [4, 5, 3, 0], [0, 0, 7, 0], [0, 0, 11, 0]]
    assert_eq(a2, expected)
    assert a2.chunks == a.chunks

    a[3:, 3:].to_zarr(z, region=(slice(2, 3), slice(1, 2)))
    a2 = da.from_zarr(z)
    expected = [[0, 1, 0, 0], [4, 5, 3, 0], [0, 15, 7, 0], [0, 0, 11, 0]]
    assert_eq(a2, expected)
    assert a2.chunks == a.chunks

    with pytest.raises(ValueError):
        with tmpdir() as d:
            a.to_zarr(d, region=(slice(2), slice(2)))


def test_tiledb_roundtrip():
    tiledb = pytest.importorskip("tiledb")
    # 1) load with default chunking
    # 2) load from existing tiledb.DenseArray
    # 3) write to existing tiledb.DenseArray
    rng = da.random.default_rng()
    a = rng.random((3, 3))
    with tmpdir() as uri:
        da.to_tiledb(a, uri)
        tdb = da.from_tiledb(uri)

        assert_eq(a, tdb)
        assert a.chunks == tdb.chunks

        # from tiledb.array
        with tiledb.open(uri) as t:
            tdb2 = da.from_tiledb(t)
            assert_eq(a, tdb2)

    with tmpdir() as uri2:
        with tiledb.empty_like(uri2, a) as t:
            a.to_tiledb(t)
            assert_eq(da.from_tiledb(uri2), a)

    # specific chunking
    with tmpdir() as uri:
        a = rng.random((3, 3), chunks=(1, 1))
        a.to_tiledb(uri)
        tdb = da.from_tiledb(uri)

        assert_eq(a, tdb)
        assert a.chunks == tdb.chunks


def test_tiledb_multiattr():
    tiledb = pytest.importorskip("tiledb")
    dom = tiledb.Domain(
        tiledb.Dim("x", (0, 1000), tile=100), tiledb.Dim("y", (0, 1000), tile=100)
    )
    schema = tiledb.ArraySchema(
        attrs=(tiledb.Attr("attr1"), tiledb.Attr("attr2")), domain=dom
    )

    with tmpdir() as uri:
        tiledb.DenseArray.create(uri, schema)
        tdb = tiledb.DenseArray(uri, "w")

        rng = np.random.default_rng()
        ar1 = rng.standard_normal(tdb.schema.shape)
        ar2 = rng.standard_normal(tdb.schema.shape)

        tdb[:] = {"attr1": ar1, "attr2": ar2}
        tdb = tiledb.DenseArray(uri, "r")

        # basic round-trip from dask.array
        d = da.from_tiledb(uri, attribute="attr2")
        assert_eq(d, ar2)

        # smoke-test computation directly on the TileDB view
        d = da.from_tiledb(uri, attribute="attr2")
        assert_eq(np.mean(ar2), d.mean().compute(scheduler="threads"))


def test_blockview():
    x = da.arange(10, chunks=2)
    blockview = BlockView(x)
    assert x.blocks == blockview
    assert isinstance(blockview[0], da.Array)

    assert_eq(blockview[0], x[:2])
    assert_eq(blockview[-1], x[-2:])
    assert_eq(blockview[:3], x[:6])
    assert_eq(blockview[[0, 1, 2]], x[:6])
    assert_eq(blockview[[3, 0, 2]], np.array([6, 7, 0, 1, 4, 5]))
    assert_eq(blockview.shape, tuple(map(len, x.chunks)))
    assert_eq(blockview.size, math.prod(blockview.shape))
    assert_eq(
        blockview.ravel(), [blockview[idx] for idx in np.ndindex(blockview.shape)]
    )

    x = da.random.default_rng().random((20, 20), chunks=(4, 5))
    blockview = BlockView(x)
    assert_eq(blockview[0], x[:4])
    assert_eq(blockview[0, :3], x[:4, :15])
    assert_eq(blockview[:, :3], x[:, :15])
    assert_eq(blockview.shape, tuple(map(len, x.chunks)))
    assert_eq(blockview.size, math.prod(blockview.shape))
    assert_eq(
        blockview.ravel(), [blockview[idx] for idx in np.ndindex(blockview.shape)]
    )

    x = da.ones((40, 40, 40), chunks=(10, 10, 10))
    blockview = BlockView(x)
    assert_eq(blockview[0, :, 0], np.ones((10, 40, 10)))
    assert_eq(blockview.shape, tuple(map(len, x.chunks)))
    assert_eq(blockview.size, math.prod(blockview.shape))
    assert_eq(
        blockview.ravel(), [blockview[idx] for idx in np.ndindex(blockview.shape)]
    )

    x = da.ones((2, 2), chunks=1)
    with pytest.raises(ValueError):
        blockview[[0, 1], [0, 1]]
    with pytest.raises(ValueError):
        blockview[np.array([0, 1]), [0, 1]]
    with pytest.raises(ValueError) as info:
        blockview[np.array([0, 1]), np.array([0, 1])]
    assert "list" in str(info.value)
    with pytest.raises(ValueError) as info:
        blockview[None, :, :]
    assert "newaxis" in str(info.value) and "not supported" in str(info.value)
    with pytest.raises(IndexError) as info:
        blockview[100, 100]


def test_blocks_indexer():
    x = da.arange(10, chunks=2)

    assert isinstance(x.blocks[0], da.Array)

    assert_eq(x.blocks[0], x[:2])
    assert_eq(x.blocks[-1], x[-2:])
    assert_eq(x.blocks[:3], x[:6])
    assert_eq(x.blocks[[0, 1, 2]], x[:6])
    assert_eq(x.blocks[[3, 0, 2]], np.array([6, 7, 0, 1, 4, 5]))

    x = da.random.default_rng().random((20, 20), chunks=(4, 5))
    assert_eq(x.blocks[0], x[:4])
    assert_eq(x.blocks[0, :3], x[:4, :15])
    assert_eq(x.blocks[:, :3], x[:, :15])

    x = da.ones((40, 40, 40), chunks=(10, 10, 10))
    assert_eq(x.blocks[0, :, 0], np.ones((10, 40, 10)))

    x = da.ones((2, 2), chunks=1)
    with pytest.raises(ValueError):
        x.blocks[[0, 1], [0, 1]]
    with pytest.raises(ValueError):
        x.blocks[np.array([0, 1]), [0, 1]]
    with pytest.raises(ValueError) as info:
        x.blocks[np.array([0, 1]), np.array([0, 1])]
    assert "list" in str(info.value)
    with pytest.raises(ValueError) as info:
        x.blocks[None, :, :]
    assert "newaxis" in str(info.value) and "not supported" in str(info.value)
    with pytest.raises(IndexError) as info:
        x.blocks[100, 100]


def test_partitions_indexer():
    # .partitions is an alias of .blocks for dask arrays
    x = da.arange(10, chunks=2)

    assert isinstance(x.partitions[0], da.Array)

    assert_eq(x.partitions[0], x[:2])
    assert_eq(x.partitions[-1], x[-2:])
    assert_eq(x.partitions[:3], x[:6])
    assert_eq(x.partitions[[0, 1, 2]], x[:6])
    assert_eq(x.partitions[[3, 0, 2]], np.array([6, 7, 0, 1, 4, 5]))

    x = da.random.default_rng().random((20, 20), chunks=(4, 5))
    assert_eq(x.partitions[0], x[:4])
    assert_eq(x.partitions[0, :3], x[:4, :15])
    assert_eq(x.partitions[:, :3], x[:, :15])

    x = da.ones((40, 40, 40), chunks=(10, 10, 10))
    assert_eq(x.partitions[0, :, 0], np.ones((10, 40, 10)))

    x = da.ones((2, 2), chunks=1)
    with pytest.raises(ValueError):
        x.partitions[[0, 1], [0, 1]]
    with pytest.raises(ValueError):
        x.partitions[np.array([0, 1]), [0, 1]]
    with pytest.raises(ValueError) as info:
        x.partitions[np.array([0, 1]), np.array([0, 1])]
    assert "list" in str(info.value)
    with pytest.raises(ValueError) as info:
        x.partitions[None, :, :]
    assert "newaxis" in str(info.value) and "not supported" in str(info.value)
    with pytest.raises(IndexError) as info:
        x.partitions[100, 100]


@pytest.mark.filterwarnings("ignore:the matrix subclass:PendingDeprecationWarning")
@pytest.mark.parametrize(
    "container", [pytest.param("array", marks=skip_if_no_sparray()), "matrix"]
)
def test_dask_array_holds_scipy_sparse_containers(container):
    pytest.importorskip("scipy.sparse")
    import scipy.sparse

    cls = scipy.sparse.csr_matrix if container == "matrix" else scipy.sparse.csr_array
    kind = scipy.sparse.spmatrix if container == "matrix" else scipy.sparse.sparray

    x = da.random.default_rng().random((1000, 10), chunks=(100, 10))
    x[x < 0.9] = 0
    xx = x.compute()
    y = x.map_blocks(cls)

    vs = y.to_delayed().flatten().tolist()
    values = dask.compute(*vs, scheduler="single-threaded")
    assert all(isinstance(v, cls) for v in values)

    yy = y.compute(scheduler="single-threaded")
    assert isinstance(yy, kind)
    assert (yy == xx).all()

    z = x.T.map_blocks(cls)
    zz = z.compute(scheduler="single-threaded")
    assert isinstance(zz, kind)
    assert (zz == xx.T).all()


@pytest.mark.parametrize(
    "index",
    [
        [5, 8],
        0,
        slice(5, 8),
        np.array([5, 8]),
        np.array([True, False] * 500),
        [True, False] * 500,
    ],
)
@pytest.mark.parametrize(
    ("sparse_module_path", "container"),
    [
        ("scipy.sparse", "csr_matrix"),
        pytest.param("scipy.sparse", "csr_array", marks=skip_if_no_sparray()),
        ("cupyx.scipy.sparse", "csr_matrix"),
    ],
)
def test_scipy_sparse_indexing(index, sparse_module_path, container):
    sp = pytest.importorskip(sparse_module_path)

    if sparse_module_path == "cupyx.scipy.sparse":
        backend = "cupy"
    else:
        backend = "numpy"

    with dask.config.set({"array.backend": backend}):
        x = da.random.default_rng().random((1000, 10), chunks=(100, 10))
        x[x < 0.9] = 0
        y = x.map_blocks(getattr(sp, container))

    assert not (
        x[index, :].compute(scheduler="single-threaded")
        != y[index, :].compute(scheduler="single-threaded")
    ).sum()


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize(
    "container", [pytest.param("array", marks=skip_if_no_sparray()), "matrix"]
)
def test_scipy_sparse_concatenate(axis, container):
    pytest.importorskip("scipy.sparse")
    import scipy.sparse

    cls = scipy.sparse.csr_matrix if container == "matrix" else scipy.sparse.csr_array

    rng = da.random.default_rng()

    xs = []
    ys = []
    for _ in range(2):
        x = rng.random((1000, 10), chunks=(100, 10))
        x[x < 0.9] = 0
        xs.append(x)
        ys.append(x.map_blocks(cls))

    z = da.concatenate(ys, axis=axis)
    z = z.compute()

    if axis == 0:
        sp_concatenate = scipy.sparse.vstack
    elif axis == 1:
        sp_concatenate = scipy.sparse.hstack
    z_expected = sp_concatenate([cls(e.compute()) for e in xs])

    assert (z != z_expected).nnz == 0


@pytest.mark.parametrize("func", [da.asarray, da.asanyarray])
@pytest.mark.parametrize("src", [[[1, 2]], np.asarray([[1, 2]]), da.asarray([[1, 2]])])
def test_scipy_sparse_asarray_like(src, func):
    """scipy.sparse.csr_matrix objects are not a valid argument for
    np.asarray(..., like=...) and require special-casing.
    """
    pytest.importorskip("scipy.sparse")
    import scipy.sparse

    mtx = scipy.sparse.csr_matrix([[3, 4, 5], [6, 7, 8]])
    like = da.from_array(mtx)

    a = func(src, like=like)
    assert isinstance(a._meta, type(mtx))
    assert isinstance(a.compute(), type(mtx))

    # Respect dtype; quietly disregard order
    a = func(src, dtype=np.float32, order="C", like=like)
    assert a.dtype == np.float32
    assert a.compute().dtype == np.float32
    assert isinstance(a._meta, type(mtx))
    assert isinstance(a.compute(), type(mtx))


def test_3851():
    with warnings.catch_warnings(record=True) as record:
        Y = da.random.default_rng().random((10, 10), chunks="auto")
        da.argmax(Y, axis=0).compute()
    assert not record


def test_3925():
    x = da.from_array(np.array(["a", "b", "c"], dtype=object), chunks=-1)
    assert (x[0] == x[0]).compute(scheduler="sync")


def test_map_blocks_large_inputs_delayed():
    a = da.ones(10, chunks=(5,))
    b = np.ones(1000000)

    c = a.map_blocks(add, b)
    assert any(b is v for v in c.dask.values())
    assert repr(dict(c.dask)).count(repr(b)[:10]) == 1  # only one occurrence

    d = a.map_blocks(lambda x, y: x + y.sum(), y=b)
    assert_eq(d, d)
    assert any(b is v for v in d.dask.values())
    assert repr(dict(c.dask)).count(repr(b)[:10]) == 1  # only one occurrence


def test_blockwise_large_inputs_delayed():
    def func(a, b):
        return a

    a = da.ones(10, chunks=(5,))
    b = np.ones(1000000)

    c = da.blockwise(func, "i", a, "i", b, None, dtype=a.dtype)
    assert any(b is v for v in c.dask.values())
    assert repr(dict(c.dask)).count(repr(b)[:10]) == 1  # only one occurrence
    assert_eq(c, c)

    d = da.blockwise(lambda x, y: x, "i", a, "i", y=b, dtype=a.dtype)
    assert any(b is v for v in d.dask.values())
    assert repr(dict(c.dask)).count(repr(b)[:10]) == 1  # only one occurrence
    assert_eq(d, d)


def test_slice_reversed():
    x = da.ones(10, chunks=-1)
    y = x[6:3]

    assert_eq(y, np.ones(0))


def test_map_blocks_chunks():
    x = da.arange(400, chunks=(100,))
    y = da.arange(40, chunks=(10,))

    def func(a, b):
        return np.array([a.max(), b.max()])

    assert_eq(
        da.map_blocks(func, x, y, chunks=(2,), dtype=x.dtype),
        np.array([99, 9, 199, 19, 299, 29, 399, 39]),
    )


def test_nbytes_auto():
    chunks = normalize_chunks("800B", shape=(500,), dtype="float64")
    assert chunks == ((100, 100, 100, 100, 100),)
    chunks = normalize_chunks("200B", shape=(10, 10), dtype="float64")
    assert chunks == ((5, 5), (5, 5))
    chunks = normalize_chunks((5, "200B"), shape=(10, 10), dtype="float64")
    assert chunks == ((5, 5), (5, 5))
    chunks = normalize_chunks("33B", shape=(10, 10), dtype="float64")
    assert chunks == ((2, 2, 2, 2, 2), (2, 2, 2, 2, 2))
    chunks = normalize_chunks("1800B", shape=(10, 20, 30), dtype="float64")
    assert chunks == ((6, 4), (6, 6, 6, 2), (6, 6, 6, 6, 6))

    with pytest.raises(ValueError):
        normalize_chunks("10B", shape=(10,), limit=20, dtype="float64")
    with pytest.raises(ValueError):
        normalize_chunks("100B", shape=(10, 10), limit=20, dtype="float64")
    with pytest.raises(ValueError):
        normalize_chunks(("100B", "10B"), shape=(10, 10), dtype="float64")
    with pytest.raises(ValueError):
        normalize_chunks(("10B", "10B"), shape=(10, 10), limit=20, dtype="float64")


def test_auto_chunks():
    chunks = ((1264, 1264, 1264, 1264, 1264, 1264, 1045), (1264, 491))
    shape = sum(chunks[0]), sum(chunks[1])
    result = normalize_chunks(
        ("auto", "auto"), shape=shape, dtype="int32", previous_chunks=chunks
    )
    assert result == ((8629,), (1755,))


def test_auto_chunks_h5py():
    h5py = pytest.importorskip("h5py")

    with tmpfile(".hdf5") as fn:
        with h5py.File(fn, mode="a") as f:
            d = f.create_dataset(
                "/x", shape=(1000, 1000), chunks=(32, 64), dtype="float64"
            )
            d[:] = 1

        with h5py.File(fn, mode="a") as f:
            d = f["x"]
            with dask.config.set({"array.chunk-size": "1 MiB"}):
                x = da.from_array(d)
                assert isinstance(x._meta, np.ndarray)
                assert x.chunks == ((256, 256, 256, 232), (512, 488))


def test_no_warnings_from_blockwise():
    with warnings.catch_warnings(record=True) as record:
        x = da.ones((3, 10, 10), chunks=(3, 2, 2))
        da.map_blocks(lambda y: np.mean(y, axis=0), x, dtype=x.dtype, drop_axis=0)
    assert not record

    with warnings.catch_warnings(record=True) as record:
        x = da.ones((15, 15), chunks=(5, 5))
        (x.dot(x.T + 1) - x.mean(axis=0)).std()
    assert not record

    with warnings.catch_warnings(record=True) as record:
        x = da.ones((1,), chunks=(1,))
        1 / x[0]
    assert not record


def test_from_array_meta():
    sparse = pytest.importorskip("sparse")
    x = np.ones(10)
    meta = sparse.COO.from_numpy(x)
    y = da.from_array(x, meta=meta)
    assert isinstance(y._meta, sparse.COO)


def test_compute_chunk_sizes():
    x = da.from_array(np.linspace(-1, 1, num=50), chunks=10)
    y = x[x < 0]
    assert np.isnan(y.shape[0])
    assert y.chunks == ((np.nan,) * 5,)

    z = y.compute_chunk_sizes()
    assert y is z
    assert z.chunks == ((10, 10, 5, 0, 0),)
    assert len(z) == 25

    # check that dtype of chunk dimensions is `int`
    assert isinstance(z.chunks[0][0], int)


def test_compute_chunk_sizes_2d_array():
    X = np.linspace(-1, 1, num=9 * 4).reshape(9, 4)
    X = da.from_array(X, chunks=(3, 4))
    idx = X.sum(axis=1) > 0
    Y = X[idx]

    # This is very similar to the DataFrame->Array conversion
    assert np.isnan(Y.shape[0]) and Y.shape[1] == 4
    assert Y.chunks == ((np.nan, np.nan, np.nan), (4,))

    Z = Y.compute_chunk_sizes()
    assert Y is Z
    assert Z.chunks == ((0, 1, 3), (4,))
    assert Z.shape == (4, 4)


def test_compute_chunk_sizes_3d_array(N=8):
    X = np.linspace(-1, 2, num=8 * 8 * 8).reshape(8, 8, 8)
    X = da.from_array(X, chunks=(4, 4, 4))
    idx = X.sum(axis=0).sum(axis=0) > 0
    Y = X[idx]
    idx = X.sum(axis=1).sum(axis=1) < 0
    Y = Y[:, idx]
    idx = X.sum(axis=2).sum(axis=1) > 0.1
    Y = Y[:, :, idx]

    # Checking to make sure shapes are different on outputs
    assert Y.compute().shape == (8, 3, 5)
    assert X.compute().shape == (8, 8, 8)

    assert Y.chunks == ((np.nan, np.nan),) * 3
    assert all(np.isnan(s) for s in Y.shape)
    Z = Y.compute_chunk_sizes()
    assert Z is Y
    assert Z.shape == (8, 3, 5)
    assert Z.chunks == ((4, 4), (3, 0), (1, 4))


def _known(num=50):
    return da.from_array(np.linspace(-1, 1, num=num), chunks=10)


@pytest.fixture()
def unknown():
    x = _known()
    y = x[x < 0]
    assert y.chunks == ((np.nan,) * 5,)
    return y


def test_compute_chunk_sizes_warning_fixes_rechunk(unknown):
    y = unknown
    with pytest.raises(ValueError, match="compute_chunk_sizes"):
        y.rechunk("auto")
    y.compute_chunk_sizes()
    y.rechunk("auto")


def test_compute_chunk_sizes_warning_fixes_to_zarr(unknown):
    pytest.importorskip("zarr")
    y = unknown
    with tmpdir() as d:
        with pytest.raises(ValueError, match="compute_chunk_sizes"):
            y.to_zarr(d)
        y.compute_chunk_sizes()
        y.to_zarr(d)


def test_compute_chunk_sizes_warning_fixes_to_svg(unknown):
    y = unknown
    with pytest.raises(NotImplementedError, match="compute_chunk_sizes"):
        y.to_svg()
    y.compute_chunk_sizes()
    y.to_svg()


def test_compute_chunk_sizes_warning_fixes_concatenate():
    x = _known(num=100).reshape(10, 10)
    idx = x.sum(axis=0) > 0
    y1 = x[idx]
    y2 = x[idx]
    with pytest.raises(ValueError, match="compute_chunk_sizes"):
        da.concatenate((y1, y2), axis=1)
    y1.compute_chunk_sizes()
    y2.compute_chunk_sizes()
    da.concatenate((y1, y2), axis=1)


def test_compute_chunk_sizes_warning_fixes_reduction(unknown):
    y = unknown
    with pytest.raises(ValueError, match="compute_chunk_sizes"):
        da.argmin(y)
    y.compute_chunk_sizes()
    da.argmin(y)


def test_compute_chunk_sizes_warning_fixes_reshape(unknown):
    y = unknown
    with pytest.raises(ValueError, match="compute_chunk_sizes"):
        da.reshape(y, (5, 5))
    y.compute_chunk_sizes()
    da.reshape(y, (5, 5))


def test_compute_chunk_sizes_warning_fixes_slicing():
    x = _known(num=100).reshape(10, 10)
    y = x[x.sum(axis=0) < 0]
    with pytest.raises(ValueError, match="compute_chunk_sizes"):
        y[:3, :]
    y.compute_chunk_sizes()
    y[:3, :]


def test_rechunk_auto():
    x = da.ones(10, chunks=(1,))
    y = x.rechunk()

    assert y.npartitions == 1


def test_chunk_assignment_invalidates_cached_properties():
    x = da.ones((4,), chunks=(1,))
    y = x.copy()
    # change chunks directly, which should change all of the tested properties
    y._chunks = ((2, 2), (0, 0, 0, 0))
    assert not x.ndim == y.ndim
    assert not x.shape == y.shape
    assert not x.size == y.size
    assert not x.numblocks == y.numblocks
    assert not x.npartitions == y.npartitions
    assert not x.__dask_keys__() == y.__dask_keys__()
    assert not np.array_equal(x._key_array, y._key_array)


def test_map_blocks_series():
    pd = pytest.importorskip("pandas")
    import dask.dataframe as dd

    pytest.skip("array roundtrips don't work yet")
    from dask.dataframe.utils import assert_eq as dd_assert_eq

    x = da.ones(10, chunks=(5,))
    s = x.map_blocks(pd.Series)
    assert isinstance(s, dd.Series)
    assert s.npartitions == x.npartitions

    dd_assert_eq(s, s)


@pytest.mark.xfail(reason="need to remove singleton index dimension")
def test_map_blocks_dataframe():
    pd = pytest.importorskip("pandas")
    import dask.dataframe as dd
    from dask.dataframe.utils import assert_eq as dd_assert_eq

    x = da.ones((10, 2), chunks=(5, 2))
    s = x.map_blocks(pd.DataFrame)
    assert isinstance(s, dd.DataFrame)
    assert s.npartitions == x.npartitions
    dd_assert_eq(s, s)


def test_dask_layers():
    a = da.ones(1)
    assert a.dask.layers.keys() == {a.name}
    assert a.dask.dependencies == {a.name: set()}
    assert a.__dask_layers__() == (a.name,)
    b = a + 1
    assert b.dask.layers.keys() == {a.name, b.name}
    assert b.dask.dependencies == {a.name: set(), b.name: {a.name}}
    assert b.__dask_layers__() == (b.name,)


def test_len_object_with_unknown_size():
    a = da.random.default_rng().random(size=(20, 2))
    b = a[a < 0.5]
    with pytest.raises(ValueError, match="on object with unknown chunk size"):
        assert len(b)


@pytest.mark.parametrize("ndim", [0, 1, 3, 8])
def test_chunk_shape_broadcast(ndim):
    from functools import partial

    def f(x, ndim=0):
        # Ignore `x` and return arbitrary one-element array of dimensionality `ndim`
        # For example,
        # f(x, 0) = array(5)
        # f(x, 1) = array([5])
        # f(x, 2) = array([[5]])
        # f(x, 3) = array([[[5]]])
        return np.array(5)[(np.newaxis,) * ndim]

    array = da.from_array([1] + [2, 2] + [3, 3, 3], chunks=((1, 2, 3),))
    out_chunks = ((1, 1, 1),)

    # check ``enforce_ndim`` keyword parameter of ``map_blocks()``
    out = array.map_blocks(partial(f, ndim=ndim), chunks=out_chunks, enforce_ndim=True)

    if ndim != 1:
        with pytest.raises(ValueError, match="Dimension mismatch:"):
            out.compute()
    else:
        out.compute()  # should not raise an exception

    # check ``check_ndim`` keyword parameter of ``assert_eq()``
    out = array.map_blocks(partial(f, ndim=ndim), chunks=out_chunks)
    expected = np.array([5, 5, 5])
    try:
        assert_eq(out, expected)
    except AssertionError:
        assert_eq(out, expected, check_ndim=False)
    else:
        if ndim != 1:
            raise AssertionError("Expected a ValueError: Dimension mismatch")


def test_chunk_non_array_like():
    array = da.from_array([1] + [2, 2] + [3, 3, 3], chunks=((1, 2, 3),))
    out_chunks = ((1, 1, 1),)

    # check ``enforce_ndim`` keyword parameter of ``map_blocks()``
    out = array.map_blocks(lambda x: 5, chunks=out_chunks, enforce_ndim=True)

    with pytest.raises(ValueError, match="Dimension mismatch:"):
        out.compute()

    expected = np.array([5, 5, 5])
    # check ``check_ndim`` keyword parameter of ``assert_eq()``
    out = array.map_blocks(lambda x: 5, chunks=out_chunks)
    try:
        assert_eq(out, expected)
    except AssertionError:
        assert_eq(out, expected, check_chunks=False)
    else:
        raise AssertionError("Expected a ValueError: Dimension mismatch")


def test_to_backend():
    # Test that `Array.to_backend` works as expected
    with dask.config.set({"array.backend": "numpy"}):
        # Start with numpy-backed array
        x = da.ones(10)
        assert isinstance(x._meta, np.ndarray)

        # Default `to_backend` shouldn't change data
        assert_eq(x, x.to_backend())

        # Moving to a "missing" backend should raise an error
        with pytest.raises(ValueError, match="No backend dispatch registered"):
            x.to_backend("missing")


def test_from_array_copies():
    x = np.arange(60).reshape((6, 10))
    original_array = x.copy()
    chunks = (2, 3)
    dx = da.from_array(x, chunks=chunks)
    x[2:4, x[0] > 3] = -5
    assert_eq(original_array, dx)


def test_from_array_xarray_dataarray():
    xr = pytest.importorskip("xarray")
    arr = xr.DataArray(da.random.random((1000, 1000), chunks=(50, 50)))
    dask_array = da.from_array(arr)
    dsk = collections_to_dsk([dask_array])
    assert len(dsk) == 400
    assert all(k[0].startswith("random_sample") for k in dsk)
    assert_eq(dask_array, arr.data)

    arr = xr.DataArray(np.random.random((100, 100)))
    dask_array = da.from_array(arr)
    assert_eq(dask_array.compute().data, arr.data)


def test_load_store_chunk():
    actual = np.array([0, 0, 0, 0, 0, 0])
    load_store_chunk(
        x=np.array([1, 2, 3]),
        out=actual,
        index=slice(2, 5),
        lock=False,
        return_stored=False,
        load_stored=False,
    )
    expected = np.array([0, 0, 1, 2, 3, 0])
    assert all(actual == expected)
    # index should not be used on empty array
    actual = load_store_chunk(
        x=np.array([]),
        out=np.array([]),
        index=2,
        lock=False,
        return_stored=True,
        load_stored=False,
    )
    expected = np.array([])
    assert all(actual == expected)


def test_scalar_setitem():
    """After a da.Array.__getitem__ call that returns a scalar, the chunk contains a
    read-only np.generic instead of a writeable np.ndarray. This is a specific quirk of
    numpy; cupy and other backends always return a 0-dimensional array.
    Make sure that __setitem__ still works.
    """
    x = da.zeros(1)
    y = x[0]
    assert isinstance(y.compute(), np.generic)
    y[()] = 2
    assert_eq(y, 2.0)
    assert isinstance(y.compute(), np.ndarray)


@pytest.mark.parametrize(
    "idx", [[0], [True, False], da.array([0]), da.array([True, False])]
)
@pytest.mark.parametrize(
    "val",
    [3.3, np.float64(3.3), np.int64(3), da.array(3.3), da.array(3, dtype=np.int64)],
)
def test_setitem_no_dtype_broadcast(idx, val):
    x = da.array([1, 2], dtype=np.int32)
    x[idx] = val
    assert_eq(x, da.array([3, 2], dtype=np.int32))
