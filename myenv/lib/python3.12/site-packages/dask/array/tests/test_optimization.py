from __future__ import annotations

import pytest

pytest.importorskip("numpy")

import numpy as np

import dask
import dask.array as da
from dask.array.chunk import getitem as da_getitem
from dask.array.core import getter as da_getter
from dask.array.core import getter_nofancy as da_getter_nofancy
from dask.array.optimization import (
    _is_getter_task,
    fuse_slice,
    optimize,
    optimize_blockwise,
    optimize_slices,
)
from dask.array.utils import assert_eq
from dask.highlevelgraph import HighLevelGraph
from dask.optimization import SubgraphCallable, fuse
from dask.utils import SerializableLock


def _wrap_getter(func, wrap):
    """
    Getters generated from a Blockwise layer might be wrapped in a SubgraphCallable.
    Make sure that the optimization functions can still work if that is the case.
    """
    if wrap:
        return SubgraphCallable({"key": (func, "index")}, outkey="key", inkeys="index")
    else:
        return func


@pytest.fixture(params=[True, False])
def getter(request):
    """
    Parameterized fixture for dask.array.core.getter both alone (False)
    and wrapped in a SubgraphCallable (True).
    """
    yield _wrap_getter(da_getter, request.param)


@pytest.fixture(params=[True, False])
def getitem(request):
    """
    Parameterized fixture for dask.array.chunk.getitem both alone (False)
    and wrapped in a SubgraphCallable (True).
    """
    yield _wrap_getter(da_getitem, request.param)


@pytest.fixture(params=[True, False])
def getter_nofancy(request):
    """
    Parameterized fixture for dask.array.chunk.getter_nofancy both alone (False)
    and wrapped in a SubgraphCallable (True).
    """
    yield _wrap_getter(da_getter_nofancy, request.param)


def _check_get_task_eq(a, b) -> bool:
    """
    Check that two tasks (possibly containing nested tasks) are equal, where
    equality is lax by allowing the callable in a SubgraphCallable to be the same
    as a non-wrapped task.
    """
    if len(a) < 1 or len(a) != len(b):
        return False

    a_callable = (
        list(a[0].dsk.values())[0][0] if isinstance(a[0], SubgraphCallable) else a[0]
    )
    b_callable = (
        list(b[0].dsk.values())[0][0] if isinstance(b[0], SubgraphCallable) else b[0]
    )
    if a_callable != b_callable:
        return False

    for ae, be in zip(a[1:], b[1:]):
        if dask.core.istask(ae):
            if not _check_get_task_eq(ae, be):
                return False
        elif ae != be:
            return False
    return True


def _assert_getter_dsk_eq(a, b):
    """
    Compare two getter dsks.

    TODO: this is here to support the fact that low-level array slice fusion needs to be
    able to introspect slicing tasks. But some slicing tasks (e.g. `from_array`) could
    be hidden within SubgraphCallables. This and _check_get_task_eq should be removed
    when high-level slicing lands, and replaced with basic equality checks.
    """
    assert a.keys() == b.keys()
    for k, av in a.items():
        bv = b[k]
        if dask.core.istask(av):
            assert _check_get_task_eq(av, bv)
        else:
            assert av == bv


def test_fuse_getitem(getter, getter_nofancy, getitem):
    pairs = [
        (
            (getter, (getter, "x", slice(1000, 2000)), slice(15, 20)),
            (getter, "x", slice(1015, 1020)),
        ),
        (
            (
                getitem,
                (getter, "x", (slice(1000, 2000), slice(100, 200))),
                (slice(15, 20), slice(50, 60)),
            ),
            (getter, "x", (slice(1015, 1020), slice(150, 160))),
        ),
        (
            (
                getitem,
                (getter_nofancy, "x", (slice(1000, 2000), slice(100, 200))),
                (slice(15, 20), slice(50, 60)),
            ),
            (getter_nofancy, "x", (slice(1015, 1020), slice(150, 160))),
        ),
        ((getter, (getter, "x", slice(1000, 2000)), 10), (getter, "x", 1010)),
        (
            (getitem, (getter, "x", (slice(1000, 2000), 10)), (slice(15, 20),)),
            (getter, "x", (slice(1015, 1020), 10)),
        ),
        (
            (getitem, (getter_nofancy, "x", (slice(1000, 2000), 10)), (slice(15, 20),)),
            (getter_nofancy, "x", (slice(1015, 1020), 10)),
        ),
        (
            (getter, (getter, "x", (10, slice(1000, 2000))), (slice(15, 20),)),
            (getter, "x", (10, slice(1015, 1020))),
        ),
        (
            (
                getter,
                (getter, "x", (slice(1000, 2000), slice(100, 200))),
                (slice(None, None), slice(50, 60)),
            ),
            (getter, "x", (slice(1000, 2000), slice(150, 160))),
        ),
        (
            (getter, (getter, "x", (None, slice(None, None))), (slice(None, None), 5)),
            (getter, "x", (None, 5)),
        ),
        (
            (
                getter,
                (getter, "x", (slice(1000, 2000), slice(10, 20))),
                (slice(5, 10),),
            ),
            (getter, "x", (slice(1005, 1010), slice(10, 20))),
        ),
        (
            (
                getitem,
                (getitem, "x", (slice(1000, 2000),)),
                (slice(5, 10), slice(10, 20)),
            ),
            (getitem, "x", (slice(1005, 1010), slice(10, 20))),
        ),
        (
            (getter, (getter, "x", slice(1000, 2000), False, False), slice(15, 20)),
            (getter, "x", slice(1015, 1020)),
        ),
        (
            (getter, (getter, "x", slice(1000, 2000)), slice(15, 20), False, False),
            (getter, "x", slice(1015, 1020)),
        ),
        (
            (
                getter,
                (getter_nofancy, "x", slice(1000, 2000), False, False),
                slice(15, 20),
                False,
                False,
            ),
            (getter_nofancy, "x", slice(1015, 1020), False, False),
        ),
    ]

    for inp, expected in pairs:
        result = optimize_slices({"y": inp})
        _assert_getter_dsk_eq(result, {"y": expected})


def test_fuse_getitem_lock(getter, getter_nofancy, getitem):
    lock1 = SerializableLock()
    lock2 = SerializableLock()

    pairs = [
        (
            (getter, (getter, "x", slice(1000, 2000), True, lock1), slice(15, 20)),
            (getter, "x", slice(1015, 1020), True, lock1),
        ),
        (
            (
                getitem,
                (getter, "x", (slice(1000, 2000), slice(100, 200)), True, lock1),
                (slice(15, 20), slice(50, 60)),
            ),
            (getter, "x", (slice(1015, 1020), slice(150, 160)), True, lock1),
        ),
        (
            (
                getitem,
                (
                    getter_nofancy,
                    "x",
                    (slice(1000, 2000), slice(100, 200)),
                    True,
                    lock1,
                ),
                (slice(15, 20), slice(50, 60)),
            ),
            (getter_nofancy, "x", (slice(1015, 1020), slice(150, 160)), True, lock1),
        ),
        (
            (
                getter,
                (getter, "x", slice(1000, 2000), True, lock1),
                slice(15, 20),
                True,
                lock2,
            ),
            (
                getter,
                (getter, "x", slice(1000, 2000), True, lock1),
                slice(15, 20),
                True,
                lock2,
            ),
        ),
    ]

    for inp, expected in pairs:
        result = optimize_slices({"y": inp})
        _assert_getter_dsk_eq(result, {"y": expected})


def test_optimize_with_getitem_fusion(getter):
    dsk = {
        "a": "some-array",
        "b": (getter, "a", (slice(10, 20), slice(100, 200))),
        "c": (getter, "b", (5, slice(50, 60))),
    }

    result = optimize(dsk, ["c"])
    expected_task = (getter, "some-array", (15, slice(150, 160)))
    assert any(_check_get_task_eq(v, expected_task) for v in result.values())
    assert len(result) < len(dsk)


def test_optimize_slicing(getter):
    dsk = {
        "a": (range, 10),
        "b": (getter, "a", (slice(None, None, None),)),
        "c": (getter, "b", (slice(None, None, None),)),
        "d": (getter, "c", (slice(0, 5, None),)),
        "e": (getter, "d", (slice(None, None, None),)),
    }

    expected = {"e": (getter, (range, 10), (slice(0, 5, None),))}
    result = optimize_slices(fuse(dsk, [], rename_keys=False)[0])
    _assert_getter_dsk_eq(result, expected)

    # protect output keys
    expected = {
        "c": (getter, (range, 10), (slice(0, None, None),)),
        "d": (getter, "c", (slice(0, 5, None),)),
        "e": (getter, "d", (slice(None, None, None),)),
    }
    result = optimize_slices(fuse(dsk, ["c", "d", "e"], rename_keys=False)[0])

    _assert_getter_dsk_eq(result, expected)


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


def test_hard_fuse_slice_cases(getter):
    dsk = {
        "x": (getter, (getter, "x", (None, slice(None, None))), (slice(None, None), 5))
    }
    _assert_getter_dsk_eq(optimize_slices(dsk), {"x": (getter, "x", (None, 5))})


def test_dont_fuse_numpy_arrays():
    x = np.ones(10)
    for _ in [(5,), (10,)]:
        y = da.from_array(x, chunks=(10,))

        dsk = y.__dask_optimize__(y.dask, y.__dask_keys__())
        assert sum(isinstance(v, np.ndarray) for v in dsk.values()) == 1


def test_fuse_slices_with_alias(getter, getitem):
    dsk = {
        "x": np.arange(16).reshape((4, 4)),
        ("dx", 0, 0): (getter, "x", (slice(0, 4), slice(0, 4))),
        ("alias", 0, 0): ("dx", 0, 0),
        ("dx2", 0): (getitem, ("alias", 0, 0), (slice(None), 0)),
    }
    keys = [("dx2", 0)]
    dsk2 = optimize(dsk, keys)
    assert len(dsk2) == 3
    fused_key = (dsk2.keys() - {"x", ("dx2", 0)}).pop()
    assert _check_get_task_eq(dsk2[fused_key], (getter, "x", (slice(0, 4), 0)))


def test_dont_fuse_fancy_indexing_in_getter_nofancy(getitem, getter_nofancy):
    dsk = {
        "a": (
            getitem,
            (getter_nofancy, "x", (slice(10, 20, None), slice(100, 200, None))),
            ([1, 3], slice(50, 60, None)),
        )
    }
    _assert_getter_dsk_eq(optimize_slices(dsk), dsk)

    dsk = {"a": (getitem, (getter_nofancy, "x", [1, 2, 3]), 0)}
    _assert_getter_dsk_eq(optimize_slices(dsk), dsk)


@pytest.mark.parametrize("chunks", [10, 5, 3])
def test_fuse_getter_with_asarray(chunks):
    x = np.ones(10) * 1234567890
    y = da.ones(10, chunks=chunks)
    z = x + y
    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())
    assert any(v is x for v in dsk.values())
    for v in dsk.values():
        s = str(v)
        assert s.count("getitem") + s.count("getter") <= 1
        if v is not x:
            assert "1234567890" not in s
    n_getters = len([v for v in dsk.values() if _is_getter_task(v)])
    if y.npartitions > 1:
        assert n_getters == y.npartitions
    else:
        assert n_getters == 0

    assert_eq(z, x + 1)


def test_remove_no_op_slices_for_getitem(getitem):
    null = slice(0, None)
    opts = [
        ((getitem, "x", null, False, False), "x"),
        ((getitem, (getitem, "x", null, False, False), null), "x"),
        ((getitem, (getitem, "x", (null, null), False, False), ()), "x"),
    ]
    for orig, final in opts:
        _assert_getter_dsk_eq(optimize_slices({"a": orig}), {"a": final})


@pytest.mark.parametrize("which", ["getter", "getter_nofancy"])
def test_dont_remove_no_op_slices_for_getter_or_getter_nofancy(
    which, getitem, getter, getter_nofancy
):
    # Test that no-op slices are *not* removed when using getter or
    # getter_nofancy. This ensures that `get` calls are always made in all
    # tasks created by `from_array`, even after optimization

    # Pytest doesn't make it easy to parameterize over parameterized fixtures
    if which == "getter":
        get = getter
    else:
        get = getter_nofancy

    null = slice(0, None)
    opts = [
        (
            (get, "x", null, False, False),
            (get, "x", null, False, False),
        ),
        (
            (getitem, (get, "x", null, False, False), null),
            (get, "x", null, False, False),
        ),
        (
            (getitem, (get, "x", (null, null), False, False), ()),
            (get, "x", (null, null), False, False),
        ),
    ]
    for orig, final in opts:
        _assert_getter_dsk_eq(optimize_slices({"a": orig}), {"a": final})


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
