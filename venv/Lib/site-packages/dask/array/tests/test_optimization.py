from __future__ import annotations

import pytest

from dask._task_spec import Alias, DataNode

pytest.importorskip("numpy")

import numpy as np

import dask
import dask.array as da
from dask.array.chunk import getitem
from dask.array.core import getter
from dask.array.optimization import fuse_slice, optimize, optimize_blockwise
from dask.array.utils import assert_eq
from dask.highlevelgraph import HighLevelGraph


def _check_get_task_eq(a, b) -> bool:
    """
    Check that two tasks (possibly containing nested tasks) are equal, where
    equality is lax by allowing the callable in a SubgraphCallable to be the same
    as a non-wrapped task.
    """
    if len(a) < 1 or len(a) != len(b):
        return False

    a_callable = a[0]
    b_callable = b[0]
    if a_callable != b_callable:
        return False

    for ae, be in zip(a[1:], b[1:]):
        if dask.core.istask(ae):
            if not _check_get_task_eq(ae, be):
                return False
        elif ae != be:
            return False
    return True


def test_optimize_with_getitem_fusion():
    dsk = {
        "a": "some-array",
        "b": (getter, "a", (slice(10, 20), slice(100, 200))),
        "c": (getter, "b", (5, slice(50, 60))),
    }

    result = optimize(dsk, ["c"])
    assert isinstance(result["c"], Alias)
    assert len(result) < len(dsk)


def test_fuse_slice():
    assert fuse_slice(slice(10, 15), slice(0, 5, 2)) == slice(10, 15, 2)

    assert fuse_slice((slice(100, 200),), (None, slice(10, 20))) == (
        None,
        slice(110, 120),
    )
    assert fuse_slice((slice(100, 200),), (slice(10, 20), None)) == (
        slice(110, 120),
        None,
    )
    assert fuse_slice((1,), (None,)) == (1, None)
    assert fuse_slice((1, slice(10, 20)), (None, None, 3, None)) == (
        1,
        None,
        None,
        13,
        None,
    )

    with pytest.raises(NotImplementedError):
        fuse_slice(slice(10, 15, 2), -1)
    # Regression test for #3076
    with pytest.raises(NotImplementedError):
        fuse_slice(None, np.array([0, 0]))


def test_fuse_slice_with_lists():
    assert fuse_slice(slice(10, 20, 2), [1, 2, 3]) == [12, 14, 16]
    assert fuse_slice([10, 20, 30, 40, 50], [3, 1, 2]) == [40, 20, 30]
    assert fuse_slice([10, 20, 30, 40, 50], 3) == 40
    assert fuse_slice([10, 20, 30, 40, 50], -1) == 50
    assert fuse_slice([10, 20, 30, 40, 50], slice(1, None, 2)) == [20, 40]
    assert fuse_slice(
        (slice(None), slice(0, 10), [1, 2, 3]), (slice(None), slice(1, 5), slice(None))
    ) == (slice(0, None), slice(1, 5), [1, 2, 3])
    assert fuse_slice(
        (slice(None), slice(None), [1, 2, 3]), (slice(None), slice(1, 5), 1)
    ) == (slice(0, None), slice(1, 5), 2)


def test_nonfusible_fancy_indexing():
    nil = slice(None)
    cases = [  # x[:, list, :][int, :, :]
        ((nil, [1, 2, 3], nil), (0, nil, nil)),
        # x[int, :, :][:, list, :]
        ((0, nil, nil), (nil, [1, 2, 3], nil)),
        # x[:, list, :, :][:, :, :, int]
        ((nil, [1, 2], nil, nil), (nil, nil, nil, 0)),
    ]

    for a, b in cases:
        with pytest.raises(NotImplementedError):
            fuse_slice(a, b)


def test_dont_fuse_numpy_arrays():
    x = np.ones(10)
    for _ in [(5,), (10,)]:
        y = da.from_array(x, chunks=(10,))

        dsk = y.__dask_optimize__(y.dask, y.__dask_keys__())
        assert (
            sum(
                isinstance(v, DataNode) and isinstance(v.value, np.ndarray)
                for v in dsk.values()
            )
            == 1
        )


def test_fuse_slices_with_alias():
    dsk = {
        "x": np.arange(16).reshape((4, 4)),
        ("dx", 0, 0): (getter, "x", (slice(0, 4), slice(0, 4))),
        ("alias", 0, 0): ("dx", 0, 0),
        ("dx2", 0): (getitem, ("alias", 0, 0), (slice(None), 0)),
    }
    keys = [("dx2", 0)]
    dsk2 = optimize(dsk, keys)
    assert len(dsk2) == 2


@pytest.mark.parametrize("chunks", [10, 5, 3])
def test_fuse_getter_with_asarray(chunks):
    x = np.ones(10) * 1234567890
    y = da.ones(10, chunks=chunks)
    z = x + y
    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())
    if chunks == 10:
        assert len(dsk) == 2 and any(isinstance(t, Alias) for t in dsk.values())
    else:
        assert any(
            isinstance(v, DataNode) and isinstance(v.value, np.ndarray)
            for v in dsk.values()
        )
    assert_eq(z, x + 1)


@pytest.mark.xfail(reason="blockwise fusion does not respect this, which is ok")
def test_turn_off_fusion():
    x = da.ones(10, chunks=(5,))
    y = da.sum(x + 1 + 2 + 3)

    a = y.__dask_optimize__(y.dask, y.__dask_keys__())

    with dask.config.set({"optimization.fuse.ave-width": 0}):
        b = y.__dask_optimize__(y.dask, y.__dask_keys__())

    assert dask.get(a, y.__dask_keys__()) == dask.get(b, y.__dask_keys__())
    assert len(a) < len(b)


def test_disable_lowlevel_fusion():
    """Check that by disabling fusion, the HLG survives through optimizations"""

    with dask.config.set({"optimization.fuse.active": False}):
        y = da.ones(3, chunks=(3,), dtype="int")
        optimize = y.__dask_optimize__
        dsk1 = y.__dask_graph__()
        dsk2 = optimize(dsk1, y.__dask_keys__())
        assert isinstance(dsk1, HighLevelGraph)
        assert isinstance(dsk2, HighLevelGraph)
        assert dsk1 == dsk2
        y = y.persist()
        assert isinstance(y.__dask_graph__(), HighLevelGraph)
        assert_eq(y, [1] * 3)


def test_array_creation_blockwise_fusion():
    """
    Check that certain array creation routines work with blockwise and can be
    fused with other blockwise operations.
    """
    x = da.ones(3, chunks=(3,))
    y = da.zeros(3, chunks=(3,))
    z = da.full(3, fill_value=2, chunks=(3,))
    a = x + y + z
    dsk1 = a.__dask_graph__()
    assert len(dsk1) == 5
    dsk2 = optimize_blockwise(dsk1)
    assert len(dsk2) == 1
    assert_eq(a, np.full(3, 3.0))


def test_gh3937():
    # test for github issue #3937
    x = da.from_array([1, 2, 3.0], (2,))
    x = da.concatenate((x, [x[-1]]))
    y = x.rechunk((2,))
    # This will produce Integral type indices that are not ints (np.int64), failing
    # the optimizer
    y = da.coarsen(np.sum, y, {0: 2})
    # How to trigger the optimizer explicitly?
    y.compute()


def test_double_dependencies():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(4, 4))
    X = d + 1
    X = da.dot(X, X.T)

    assert_eq(X.compute(optimize_graph=False), X)


def test_fuse_roots():
    x = da.ones(10, chunks=(2,))
    y = da.zeros(10, chunks=(2,))
    z = (x + 1) + (2 * y**2)
    (zz,) = dask.optimize(z)
    # assert len(zz.dask) == 5
    assert sum(map(dask.istask, zz.dask.values())) == 5  # there are some aliases
    assert_eq(zz, z)


def test_fuse_roots_annotations():
    x = da.ones(10, chunks=(2,))
    y = da.zeros(10, chunks=(2,))

    with dask.annotate(foo="bar"):
        y = y**2

    z = (x + 1) + (2 * y)
    hlg = dask.blockwise.optimize_blockwise(z.dask)
    assert len(hlg.layers) == 3
    assert {"foo": "bar"} in [l.annotations for l in hlg.layers.values()]
    za = da.Array(hlg, z.name, z.chunks, z.dtype)
    assert_eq(za, z)


@pytest.mark.parametrize("optimize_graph", [True, False])
def test_optimize_blockwise_duplicate_dependency(optimize_graph):
    # Two blockwise operations in a row with duplicate name
    # (See: https://github.com/dask/dask/issues/8535)
    xx = da.from_array(np.array([[1, 1], [2, 2]]), chunks=1)
    xx = xx * 2
    z = da.matmul(xx, xx)

    # Compare to known answer
    result = z.compute(optimize_graph=optimize_graph)
    assert assert_eq(result, [[12, 12], [24, 24]])
