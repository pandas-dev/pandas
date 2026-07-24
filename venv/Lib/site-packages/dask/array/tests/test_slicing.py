from __future__ import annotations

import itertools
import warnings

import pytest

from dask._task_spec import Alias, Task, TaskRef
from dask.delayed import delayed

np = pytest.importorskip("numpy")

import dask
import dask.array as da
from dask.array.chunk import getitem
from dask.array.slicing import (
    SlicingNoop,
    _sanitize_index_element,
    _slice_1d,
    make_block_sorted_slices,
    new_blockdim,
    normalize_index,
    sanitize_index,
    shuffle_slice,
    slice_array,
    take,
)
from dask.array.utils import assert_eq, same_keys


def test_slice_1d():
    expected = {0: slice(10, 25, 1), 1: slice(None, None, None), 2: slice(0, 1, 1)}
    result = _slice_1d(100, [25] * 4, slice(10, 51, None))
    assert expected == result

    # x[100:12:-3]
    expected = {
        0: slice(-2, -8, -3),
        1: slice(-1, -21, -3),
        2: slice(-3, -21, -3),
        3: slice(-2, -21, -3),
        4: slice(-1, -21, -3),
    }
    result = _slice_1d(100, [20] * 5, slice(100, 12, -3))
    assert expected == result

    # x[102::-3]
    expected = {
        0: slice(-2, -21, -3),
        1: slice(-1, -21, -3),
        2: slice(-3, -21, -3),
        3: slice(-2, -21, -3),
        4: slice(-1, -21, -3),
    }
    result = _slice_1d(100, [20] * 5, slice(102, None, -3))
    assert expected == result

    # x[::-4]
    expected = {
        0: slice(-1, -21, -4),
        1: slice(-1, -21, -4),
        2: slice(-1, -21, -4),
        3: slice(-1, -21, -4),
        4: slice(-1, -21, -4),
    }
    result = _slice_1d(100, [20] * 5, slice(None, None, -4))
    assert expected == result

    # x[::-7]
    expected = {
        0: slice(-5, -21, -7),
        1: slice(-4, -21, -7),
        2: slice(-3, -21, -7),
        3: slice(-2, -21, -7),
        4: slice(-1, -21, -7),
    }
    result = _slice_1d(100, [20] * 5, slice(None, None, -7))
    assert expected == result

    # x=range(115)
    # x[::-7]
    expected = {
        0: slice(-7, -24, -7),
        1: slice(-2, -24, -7),
        2: slice(-4, -24, -7),
        3: slice(-6, -24, -7),
        4: slice(-1, -24, -7),
    }
    result = _slice_1d(115, [23] * 5, slice(None, None, -7))
    assert expected == result

    # x[79::-3]
    expected = {
        0: slice(-1, -21, -3),
        1: slice(-3, -21, -3),
        2: slice(-2, -21, -3),
        3: slice(-1, -21, -3),
    }
    result = _slice_1d(100, [20] * 5, slice(79, None, -3))
    assert expected == result

    # x[-1:-8:-1]
    expected = {4: slice(-1, -8, -1)}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(-1, 92, -1))
    assert expected == result

    # x[20:0:-1]
    expected = {0: slice(-1, -20, -1), 1: slice(-20, -21, -1)}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(20, 0, -1))
    assert expected == result

    # x[:0]
    expected = {}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(0))
    assert result

    # x=range(99)
    expected = {
        0: slice(-3, -21, -3),
        1: slice(-2, -21, -3),
        2: slice(-1, -21, -3),
        3: slice(-2, -20, -3),
        4: slice(-1, -21, -3),
    }
    # This array has non-uniformly sized blocks
    result = _slice_1d(99, [20, 20, 20, 19, 20], slice(100, None, -3))
    assert expected == result

    # x=range(104)
    # x[::-3]
    expected = {
        0: slice(-1, -21, -3),
        1: slice(-3, -24, -3),
        2: slice(-3, -28, -3),
        3: slice(-1, -14, -3),
        4: slice(-1, -22, -3),
    }
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(None, None, -3))
    assert expected == result

    # x=range(104)
    # x[:27:-3]
    expected = {
        1: slice(-3, -16, -3),
        2: slice(-3, -28, -3),
        3: slice(-1, -14, -3),
        4: slice(-1, -22, -3),
    }
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(None, 27, -3))
    assert expected == result

    # x=range(104)
    # x[100:27:-3]
    expected = {
        1: slice(-3, -16, -3),
        2: slice(-3, -28, -3),
        3: slice(-1, -14, -3),
        4: slice(-4, -22, -3),
    }
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(100, 27, -3))
    assert expected == result

    # x=range(1000000000000)
    # x[1000:]
    expected = {0: slice(1000, 1000000000, 1)}
    expected.update({ii: slice(None, None, None) for ii in range(1, 1000)})
    # This array is large
    result = _slice_1d(1000000000000, [1000000000] * 1000, slice(1000, None, None))
    assert expected == result


def test_slice_singleton_value_on_boundary():
    assert _slice_1d(15, [5, 5, 5], 10) == {2: 0}
    assert _slice_1d(30, (5, 5, 5, 5, 5, 5), 10) == {2: 0}


def test_slice_array_1d():
    # x[24::2]
    expected = {
        ("y", 0): Task(("y", 0), getitem, TaskRef(("x", 0)), (slice(24, 25, 2),)),
        ("y", 1): Task(("y", 1), getitem, TaskRef(("x", 1)), (slice(1, 25, 2),)),
        ("y", 2): Task(("y", 2), getitem, TaskRef(("x", 2)), (slice(0, 25, 2),)),
        ("y", 3): Task(("y", 3), getitem, TaskRef(("x", 3)), (slice(1, 25, 2),)),
    }
    result, chunks = slice_array("y", "x", [[25] * 4], [slice(24, None, 2)])

    assert expected == result

    # x[26::2]
    expected = {
        ("y", 0): Task(("y", 0), getitem, TaskRef(("x", 1)), (slice(1, 25, 2),)),
        ("y", 1): Task(("y", 1), getitem, TaskRef(("x", 2)), (slice(0, 25, 2),)),
        ("y", 2): Task(("y", 2), getitem, TaskRef(("x", 3)), (slice(1, 25, 2),)),
    }

    result, chunks = slice_array("y", "x", [[25] * 4], [slice(26, None, 2)])
    assert expected == result

    # x[24::2]
    expected = {
        ("y", 0): Task(("y", 0), getitem, TaskRef(("x", 0)), (slice(24, 25, 2),)),
        ("y", 1): Task(("y", 1), getitem, TaskRef(("x", 1)), (slice(1, 25, 2),)),
        ("y", 2): Task(("y", 2), getitem, TaskRef(("x", 2)), (slice(0, 25, 2),)),
        ("y", 3): Task(("y", 3), getitem, TaskRef(("x", 3)), (slice(1, 25, 2),)),
    }
    result, chunks = slice_array("y", "x", [(25,) * 4], (slice(24, None, 2),))

    assert expected == result

    # x[26::2]
    expected = {
        ("y", 0): Task(("y", 0), getitem, TaskRef(("x", 1)), (slice(1, 25, 2),)),
        ("y", 1): Task(("y", 1), getitem, TaskRef(("x", 2)), (slice(0, 25, 2),)),
        ("y", 2): Task(("y", 2), getitem, TaskRef(("x", 3)), (slice(1, 25, 2),)),
    }

    result, chunks = slice_array("y", "x", [(25,) * 4], (slice(26, None, 2),))
    assert expected == result


def test_slice_array_2d():
    # 2d slices: x[13::2,10::1]
    expected = {
        ("y", 0, 0): Task(
            ("y", 0, 0),
            getitem,
            TaskRef(("x", 0, 0)),
            (slice(13, 20, 2), slice(10, 20, 1)),
        ),
        ("y", 0, 1): Task(
            ("y", 0, 1), getitem, TaskRef(("x", 0, 1)), (slice(13, 20, 2), slice(None))
        ),
        ("y", 0, 2): Task(
            ("y", 0, 2), getitem, TaskRef(("x", 0, 2)), (slice(13, 20, 2), slice(None))
        ),
    }

    result, chunks = slice_array(
        "y",
        "x",
        [[20], [20, 20, 5]],
        [slice(13, None, 2), slice(10, None, 1)],
    )

    assert expected == result

    # 2d slices with one dimension: x[5,10::1]
    expected = {
        ("y", 0): Task(("y", 0), getitem, TaskRef(("x", 0, 0)), (5, slice(10, 20, 1))),
        ("y", 1): Task(("y", 1), getitem, TaskRef(("x", 0, 1)), (5, slice(None))),
        ("y", 2): Task(("y", 2), getitem, TaskRef(("x", 0, 2)), (5, slice(None))),
    }

    result, chunks = slice_array("y", "x", ([20], [20, 20, 5]), [5, slice(10, None, 1)])

    assert expected == result


def test_mixed_index():
    da_array = da.ones((1, 1, 31, 40))
    new = da_array[(np.array([0]), 0, slice(None), slice(None))]
    assert isinstance(new, da.Array)
    assert_eq(new, da_array[0])


def test_slice_optimizations():
    # bar[:]
    with pytest.raises(SlicingNoop):
        slice_array("foo", "bar", [[100]], (slice(None, None, None),))

    # bar[:,:,:]
    with pytest.raises(SlicingNoop):
        slice_array(
            "foo",
            "bar",
            [(100, 1000, 10000)],
            (slice(None, None, None), slice(None, None, None), slice(None, None, None)),
        )


def test_slicing_with_singleton_indices():
    result, chunks = slice_array("y", "x", ([5, 5], [5, 5]), (slice(0, 5), 8))

    expected = {
        ("y", 0): Task(
            ("y", 0), getitem, TaskRef(("x", 0, 1)), (slice(None, None, None), 3)
        )
    }

    assert expected == result


def test_slicing_with_newaxis():
    result, chunks = slice_array(
        "y",
        "x",
        ([5, 5], [5, 5]),
        (slice(0, 3), None, slice(None, None, None)),
    )

    expected = {
        ("y", 0, 0, 0): Task(
            ("y", 0, 0, 0),
            getitem,
            TaskRef(("x", 0, 0)),
            (slice(0, 3, 1), None, slice(None, None, None)),
        ),
        ("y", 0, 0, 1): Task(
            ("y", 0, 0, 1),
            getitem,
            TaskRef(("x", 0, 1)),
            (slice(0, 3, 1), None, slice(None, None, None)),
        ),
    }

    assert expected == result
    assert chunks == ((3,), (1,), (5, 5))


def test_take():
    chunks, dsk = take("y-y", "x", [(20, 20, 20, 20)], [5, 1, 47, 3], axis=0)
    assert len(dsk) == 6
    assert chunks == ((4,),)

    chunks, dsk = take("y-y", "x", [(20, 20, 20, 20), (20, 20)], [5, 1, 47, 3], axis=0)
    assert len(dsk) == 9
    assert chunks == ((4,), (20, 20))


def test_take_sorted():
    chunks, dsk = take("y-y", "x", [(20, 20, 20, 20)], [1, 3, 5, 47], axis=0)
    assert len(dsk) == 6
    assert chunks == ((4,),)
    chunks, dsk = take("y", "x", [(20, 20, 20, 20)], np.arange(0, 80), axis=0)
    assert len(dsk) == 4
    assert chunks == ((20, 20, 20, 20),)


def test_slicing_chunks():
    result, chunks = slice_array("y", "x", ([5, 5], [5, 5]), (1, np.array([2, 0, 3])))
    assert chunks == ((3,),)

    result, chunks = slice_array(
        "y", "x", ([5, 5], [5, 5]), (slice(0, 7), np.array([2, 0, 3]))
    )
    assert chunks == ((5, 2), (3,))

    result, chunks = slice_array("y", "x", ([5, 5], [5, 5]), (slice(0, 7), 1))
    assert chunks == ((5, 2),)


def test_slicing_with_numpy_arrays():
    a, bd1 = slice_array(
        "y-y",
        "x",
        ((3, 3, 3, 1), (3, 3, 3, 1)),
        (np.array([1, 2, 9]), slice(None, None, None)),
    )
    b, bd2 = slice_array(
        "y-y",
        "x",
        ((3, 3, 3, 1), (3, 3, 3, 1)),
        (np.array([1, 2, 9]), slice(None, None, None)),
    )

    assert bd1 == bd2
    np.testing.assert_equal(a, b)

    i = [False, True, True, False, False, False, False, False, False, True]
    index = (i, slice(None, None, None))
    index = normalize_index(index, (10, 10))
    c, bd3 = slice_array("y-y", "x", ((3, 3, 3, 1), (3, 3, 3, 1)), index)
    assert bd1 == bd3
    np.testing.assert_equal(a, c)


def test_slicing_and_chunks():
    o = da.ones((24, 16), chunks=((4, 8, 8, 4), (2, 6, 6, 2)))
    t = o[4:-4, 2:-2]
    assert t.chunks == ((8, 8), (6, 6))


def test_slicing_and_unknown_chunks():
    a = da.ones((10, 5), chunks=5)
    a._chunks = ((np.nan, np.nan), (5,))
    with pytest.raises(ValueError, match="Array chunk size or shape is unknown"):
        a[[0, 5]].compute()


def test_slicing_identities():
    a = da.ones((24, 16), chunks=((4, 8, 8, 4), (2, 6, 6, 2)))

    assert a is a[slice(None)]
    assert a is a[:]
    assert a is a[::]
    assert a is a[...]
    assert a is a[0:]
    assert a is a[0::]
    assert a is a[::1]
    assert a is a[0 : len(a)]
    assert a is a[0::1]
    assert a is a[0 : len(a) : 1]


def test_slice_stop_0():
    # from gh-125
    a = da.ones(10, chunks=(10,))[:0].compute()
    b = np.ones(10)[:0]
    assert_eq(a, b)


def test_slice_list_then_None():
    x = da.zeros(shape=(5, 5), chunks=(3, 3))
    y = x[[2, 1]][None]

    assert_eq(y, np.zeros((1, 2, 5)))


class ReturnItem:
    def __getitem__(self, key):
        return key


@pytest.mark.skip(reason="really long test")
def test_slicing_exhaustively():
    x = np.random.default_rng().random(6, 7, 8)
    a = da.from_array(x, chunks=(3, 3, 3))
    I = ReturnItem()

    # independent indexing along different axes
    indexers = [0, -2, I[:], I[:5], [0, 1], [0, 1, 2], [4, 2], I[::-1], None, I[:0], []]
    for i in indexers:
        assert_eq(x[i], a[i])
        for j in indexers:
            assert_eq(x[i][:, j], a[i][:, j])
            assert_eq(x[:, i][j], a[:, i][j])
            for k in indexers:
                assert_eq(x[..., i][:, j][k], a[..., i][:, j][k])

    # repeated indexing along the first axis
    first_indexers = [I[:], I[:5], np.arange(5), [3, 1, 4, 5, 0], np.arange(6) < 6]
    second_indexers = [0, -1, 3, I[:], I[:3], I[2:-1], [2, 4], [], I[:0]]
    for i in first_indexers:
        for j in second_indexers:
            assert_eq(x[i][j], a[i][j])


def test_slicing_with_negative_step_flops_keys():
    x = da.arange(10, chunks=5)
    y = x[:1:-1]
    assert (x.name, 1) in y.dask[(y.name, 0)].dependencies
    assert (x.name, 0) in y.dask[(y.name, 1)].dependencies

    assert_eq(y, np.arange(10)[:1:-1])

    assert y.chunks == ((5, 3),)

    assert y.dask[(y.name, 0)] == Task(
        (y.name, 0), getitem, TaskRef((x.name, 1)), (slice(-1, -6, -1),)
    )
    assert y.dask[(y.name, 1)] == Task(
        (y.name, 1), getitem, TaskRef((x.name, 0)), (slice(-1, -4, -1),)
    )


def test_empty_slice():
    x = da.ones((5, 5), chunks=(2, 2), dtype="i4")
    y = x[:0]

    assert_eq(y, np.ones((5, 5), dtype="i4")[:0])


def test_multiple_list_slicing():
    x = np.random.default_rng().random((6, 7, 8))
    a = da.from_array(x, chunks=(3, 3, 3))
    assert_eq(x[:, [0, 1, 2]][[0, 1]], a[:, [0, 1, 2]][[0, 1]])


def test_boolean_list_slicing():
    with pytest.raises(IndexError):
        da.asarray(range(2))[[True]]
    with pytest.raises(IndexError):
        da.asarray(range(2))[[False, False, False]]
    x = np.arange(5)
    ind = [True, False, False, False, True]
    assert_eq(da.asarray(x)[ind], x[ind])
    # https://github.com/dask/dask/issues/3706
    ind = [True]
    assert_eq(da.asarray([0])[ind], np.arange(1)[ind])


def test_boolean_numpy_array_slicing():
    with pytest.raises(IndexError):
        da.asarray(range(2))[np.array([True])]
    with pytest.raises(IndexError):
        da.asarray(range(2))[np.array([False, False, False])]
    x = np.arange(5)
    ind = np.array([True, False, False, False, True])
    assert_eq(da.asarray(x)[ind], x[ind])
    # https://github.com/dask/dask/issues/3706
    ind = np.array([True])
    assert_eq(da.asarray([0])[ind], np.arange(1)[ind])


def test_empty_list():
    x = np.ones((5, 5, 5), dtype="i4")
    dx = da.from_array(x, chunks=2)

    assert_eq(dx[[], :3, :2], x[[], :3, :2])
    assert_eq(dx[:3, [], :2], x[:3, [], :2])
    assert_eq(dx[:3, :2, []], x[:3, :2, []])


def test_uneven_chunks():
    assert da.ones(20, chunks=5)[::2].chunks == ((3, 2, 3, 2),)


def test_new_blockdim():
    assert new_blockdim(20, [5, 5, 5, 5], slice(0, None, 2)) == [3, 2, 3, 2]


def test_slicing_consistent_names():
    x = np.arange(100).reshape((10, 10))
    a = da.from_array(x, chunks=(5, 5))
    assert same_keys(a[0], a[0])
    assert same_keys(a[:, [1, 2, 3]], a[:, [1, 2, 3]])
    assert same_keys(a[:, 5:2:-1], a[:, 5:2:-1])
    assert same_keys(a[0, ...], a[0, ...])
    assert same_keys(a[...], a[...])
    assert same_keys(a[[1, 3, 5]], a[[1, 3, 5]])
    assert same_keys(a[-11:11], a[:])
    assert same_keys(a[-11:-9], a[:1])
    assert same_keys(a[-1], a[9])
    assert same_keys(a[0::-1], a[0:-11:-1])


def test_slicing_consistent_names_after_normalization():
    x = da.zeros(10, chunks=(5,))
    assert same_keys(x[0:], x[:10])
    assert same_keys(x[0:], x[0:10])
    assert same_keys(x[0:], x[0:10:1])
    assert same_keys(x[:], x[0:10:1])


def test_sanitize_index_element():
    with pytest.raises(TypeError):
        _sanitize_index_element("Hello!")


def test_sanitize_index():
    pd = pytest.importorskip("pandas")
    with pytest.raises(TypeError):
        sanitize_index("Hello!")

    np.testing.assert_equal(sanitize_index(pd.Series([1, 2, 3])), [1, 2, 3])
    np.testing.assert_equal(sanitize_index((1, 2, 3)), [1, 2, 3])


def test_uneven_blockdims():
    blockdims = ((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30), (100,))
    index = (slice(240, 270), slice(None))
    dsk_out, bd_out = slice_array("in", "out", blockdims, index)
    sol = {
        ("in", 0, 0): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 7, 0)),
            (slice(28, 31, 1), slice(None)),
        ),
        ("in", 1, 0): Task(
            ("in", 1, 0),
            getitem,
            TaskRef(("out", 8, 0)),
            (slice(0, 27, 1), slice(None)),
        ),
    }
    assert dsk_out == sol
    assert bd_out == ((3, 27), (100,))

    blockdims = ((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30),) * 2
    index = (slice(240, 270), slice(180, 230))
    dsk_out, bd_out = slice_array("in", "out", blockdims, index)
    sol = {
        ("in", 0, 0): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 7, 5)),
            (slice(28, 31, 1), slice(29, 30, 1)),
        ),
        ("in", 0, 1): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 7, 6)),
            (slice(28, 31, 1), slice(None)),
        ),
        ("in", 0, 2): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 7, 7)),
            (slice(28, 31, 1), slice(0, 18, 1)),
        ),
        ("in", 1, 0): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 8, 5)),
            (slice(0, 27, 1), slice(29, 30, 1)),
        ),
        ("in", 1, 1): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 8, 6)),
            (slice(0, 27, 1), slice(None)),
        ),
        ("in", 1, 2): Task(
            ("in", 0, 0),
            getitem,
            TaskRef(("out", 8, 7)),
            (slice(0, 27, 1), slice(0, 18, 1)),
        ),
    }
    assert dsk_out == sol
    assert bd_out == ((3, 27), (1, 31, 18))


def test_oob_check():
    x = da.ones(5, chunks=(2,))
    with pytest.raises(IndexError):
        x[6]
    with pytest.raises(IndexError):
        x[[6]]
    with pytest.raises(IndexError):
        x[-10]
    with pytest.raises(IndexError):
        x[[-10]]
    with pytest.raises(IndexError):
        x[0, 0]


@pytest.mark.parametrize("idx_chunks", [None, 3, 2, 1])
@pytest.mark.parametrize("x_chunks", [None, (3, 5), (2, 3), (1, 2), (1, 1)])
def test_index_with_int_dask_array(x_chunks, idx_chunks):
    # test data is crafted to stress use cases:
    # - pick from different chunks of x out of order
    # - a chunk of x contains no matches
    # - only one chunk of x
    x = np.array(
        [[10, 20, 30, 40, 50], [60, 70, 80, 90, 100], [110, 120, 130, 140, 150]]
    )
    idx = np.array([3, 0, 1])
    expect = np.array([[40, 10, 20], [90, 60, 70], [140, 110, 120]])

    if x_chunks is not None:
        x = da.from_array(x, chunks=x_chunks)
    if idx_chunks is not None:
        idx = da.from_array(idx, chunks=idx_chunks)

    assert_eq(x[:, idx], expect)
    assert_eq(x.T[idx, :], expect.T)


@pytest.mark.parametrize("chunks", [1, 2, 3])
def test_index_with_int_dask_array_0d(chunks):
    # Slice by 0-dimensional array
    x = da.from_array([[10, 20, 30], [40, 50, 60]], chunks=chunks)
    idx0 = da.from_array(1, chunks=1)
    assert_eq(x[idx0, :], x[1, :])
    assert_eq(x[:, idx0], x[:, 1])


@pytest.mark.parametrize("chunks", [1, 2, 3, 4, 5])
def test_index_with_int_dask_array_nanchunks(chunks):
    # Slice by array with nan-sized chunks
    a = da.arange(-2, 3, chunks=chunks)
    assert_eq(a[a.nonzero()], np.array([-2, -1, 1, 2]))
    # Edge case: the nan-sized chunks resolve to size 0
    a = da.zeros(5, chunks=chunks)
    assert_eq(a[a.nonzero()], np.array([]))


@pytest.mark.parametrize("chunks", [2, 4])
def test_index_with_int_dask_array_negindex(chunks):
    a = da.arange(4, chunks=chunks)
    idx = da.from_array([-1, -4], chunks=1)
    assert_eq(a[idx], np.array([3, 0]))


@pytest.mark.parametrize("chunks", [2, 4])
def test_index_with_int_dask_array_indexerror(chunks):
    a = da.arange(4, chunks=chunks)
    idx = da.from_array([4], chunks=1)
    with pytest.raises(IndexError):
        a[idx].compute()
    idx = da.from_array([-5], chunks=1)
    with pytest.raises(IndexError):
        a[idx].compute()


@pytest.mark.parametrize(
    "dtype", ["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]
)
def test_index_with_int_dask_array_dtypes(dtype):
    a = da.from_array([10, 20, 30, 40], chunks=-1)
    idx = da.from_array(np.array([1, 2]).astype(dtype), chunks=1)
    assert_eq(a[idx], np.array([20, 30]))


def test_index_with_int_dask_array_nocompute():
    """Test that when the indices are a dask array
    they are not accidentally computed
    """

    def crash():
        raise NotImplementedError()

    x = da.arange(5, chunks=-1)
    idx = da.Array({("x", 0): (crash,)}, name="x", chunks=((2,),), dtype=np.int64)
    result = x[idx]
    with pytest.raises(NotImplementedError):
        result.compute()


def test_index_with_bool_dask_array():
    x = np.arange(36).reshape((6, 6))
    d = da.from_array(x, chunks=(3, 3))
    ind = np.asarray([True, True, False, True, False, False], dtype=bool)
    ind = da.from_array(ind, chunks=2)
    for index in [ind, (slice(1, 9, 2), ind), (ind, slice(2, 8, 1))]:
        x_index = dask.compute(index)[0]
        assert_eq(x[x_index], d[index])


def test_index_with_bool_dask_array_2():
    rng = np.random.default_rng()
    x = rng.random((10, 10, 10))
    ind = rng.random(10) > 0.5

    d = da.from_array(x, chunks=(3, 4, 5))
    dind = da.from_array(ind, chunks=4)

    index = [slice(1, 9, 1), slice(None)]

    for i in range(x.ndim):
        index2 = index[:]
        index2.insert(i, dind)

        index3 = index[:]
        index3.insert(i, ind)

        assert_eq(x[tuple(index3)], d[tuple(index2)])


@pytest.mark.xfail
def test_cull():
    x = da.ones(1000, chunks=(10,))

    for slc in [1, slice(0, 30), slice(0, None, 100)]:
        y = x[slc]
        assert len(y.dask) < len(x.dask)


@pytest.mark.parametrize("shape", [(2,), (2, 3), (2, 3, 5)])
@pytest.mark.parametrize(
    "index", [(Ellipsis,), (None, Ellipsis), (Ellipsis, None), (None, Ellipsis, None)]
)
def test_slicing_with_Nones(shape, index):
    x = np.random.default_rng().random(shape)
    d = da.from_array(x, chunks=shape)

    assert_eq(x[index], d[index])


indexers = [Ellipsis, slice(2), 0, 1, -2, -1, slice(-2, None), None]


"""
# We comment this out because it is 4096 tests
@pytest.mark.parametrize('a', indexers)
@pytest.mark.parametrize('b', indexers)
@pytest.mark.parametrize('c', indexers)
@pytest.mark.parametrize('d', indexers)
def test_slicing_none_int_ellipses(a, b, c, d):
    if (a, b, c, d).count(Ellipsis) > 1:
        return
    shape = (2,3,5,7,11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)

    xx = x[a, b, c, d]
    yy = y[a, b, c, d]
    assert_eq(xx, yy)
"""


def test_slicing_integer_no_warnings():
    # https://github.com/dask/dask/pull/2457/
    X = da.random.default_rng().random(size=(100, 2), chunks=(2, 2))
    idx = np.array([0, 0, 1, 1])
    with warnings.catch_warnings(record=True) as record:
        X[idx].compute()
    assert not record


@pytest.mark.slow
def test_slicing_none_int_ellipes():
    shape = (2, 3, 5, 7, 11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)
    for ind in itertools.product(indexers, indexers, indexers, indexers):
        if ind.count(Ellipsis) > 1:
            continue

        assert_eq(x[ind], y[ind])


def test_None_overlap_int():
    a, b, c, d = (0, slice(None, 2, None), None, Ellipsis)
    shape = (2, 3, 5, 7, 11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)

    xx = x[a, b, c, d]
    yy = y[a, b, c, d]
    assert_eq(xx, yy)


def test_negative_n_slicing():
    assert_eq(da.ones(2, chunks=2)[-2], np.ones(2)[-2])


def test_negative_list_slicing():
    x = np.arange(5)
    dx = da.from_array(x, chunks=2)
    assert_eq(dx[[0, -5]], x[[0, -5]])
    assert_eq(dx[[4, -1]], x[[4, -1]])


def test_permit_oob_slices():
    x = np.arange(5)
    dx = da.from_array(x, chunks=2)

    assert_eq(x[-102:], dx[-102:])
    assert_eq(x[102:], dx[102:])
    assert_eq(x[:102], dx[:102])
    assert_eq(x[:-102], dx[:-102])


def test_normalize_index():
    assert normalize_index((Ellipsis, None), (10,)) == (slice(None), None)
    assert normalize_index(5, (np.nan,)) == (5,)
    assert normalize_index(-5, (np.nan,)) == (-5,)
    (result,) = normalize_index([-5, -2, 1], (np.nan,))
    assert result.tolist() == [-5, -2, 1]
    assert normalize_index(slice(-5, -2), (np.nan,)) == (slice(-5, -2),)


def test_take_semi_sorted():
    x = da.ones(10, chunks=(5,))
    index = np.arange(15) % 10

    y = x[index]
    assert y.chunks == ((5, 5, 5),)


def test_getitem_avoids_large_chunks():
    with dask.config.set({"array.chunk-size": "0.1Mb"}):
        a = np.arange(2 * 128 * 128, dtype="int64").reshape(2, 128, 128)
        indexer = [0] + [1] * 11
        arr = da.from_array(a, chunks=(1, 8, 8))
        result = arr[
            indexer
        ]  # small chunks within the chunk-size limit should NOT raise PerformanceWarning
        expected = a[indexer]
        assert_eq(result, expected)

        arr = da.from_array(a, chunks=(1, 128, 128))  # large chunks
        expected = a[indexer]

        result = arr[indexer]
        assert_eq(result, expected)
        assert result.chunks == ((1,) * 12, (128,), (128,))

        # Users can silence the warning
        with dask.config.set({"array.slicing.split-large-chunks": False}):
            with warnings.catch_warnings(record=True) as record:
                result = arr[indexer]
            assert_eq(result, expected)
            assert not record

        # Users can silence the warning
        with dask.config.set({"array.slicing.split-large-chunks": True}):
            with warnings.catch_warnings(record=True) as record:
                result = arr[indexer]
            assert_eq(result, expected)
            assert not record

            assert result.chunks == ((1,) * 12, (128,), (128,))


def test_getitem_avoids_large_chunks_missing():
    # We cannot apply the "avoid large chunks" optimization when
    # the chunks have unknown sizes.
    with dask.config.set({"array.chunk-size": "0.1Mb"}):
        a = np.arange(4 * 500 * 500).reshape(4, 500, 500)
        arr = da.from_array(a, chunks=(1, 500, 500))
        arr._chunks = ((1, 1, 1, 1), (np.nan,), (np.nan,))
        indexer = [0, 1] + [2] * 100 + [3]
        expected = a[indexer]
        result = arr[indexer]
        assert_eq(result, expected)


def test_take_avoids_large_chunks():
    # unit test for https://github.com/dask/dask/issues/6270
    with dask.config.set({"array.slicing.split-large-chunks": True}):
        chunks = ((1, 1, 1, 1), (500,), (500,))
        index = np.array([0, 1] + [2] * 101 + [3])
        chunks2, dsk = take("a", "b", chunks, index)
        assert chunks2 == ((1,) * 104, (500,), (500,))
        assert len(dsk) == 105

        index = np.array([0] * 101 + [1, 2, 3])
        chunks2, dsk = take("a", "b", chunks, index)
        assert chunks2 == ((1,) * 104, (500,), (500,))
        assert len(dsk) == 105

        index = np.array([0, 1, 2] + [3] * 101)
        chunks2, dsk = take("a", "b", chunks, index)
        assert chunks2 == ((1,) * 104, (500,), (500,))
        assert len(dsk) == 105

        chunks = ((500,), (1, 1, 1, 1), (500,))
        index = np.array([0, 1, 2] + [3] * 101)
        chunks2, dsk = take("a", "b", chunks, index, axis=1)
        assert chunks2 == ((500,), (1,) * 104, (500,))
        assert len(dsk) == 105


def test_pathological_unsorted_slicing():
    x = da.ones(100, chunks=10)

    # [0, 10, 20, ... 90, 1, 11, 21, ... 91, ...]
    index = np.arange(100).reshape(10, 10).ravel(order="F")

    assert_eq(x[index], x.compute()[index])


@pytest.mark.parametrize("params", [(2, 2, 1), (5, 3, 2)])
def test_setitem_with_different_chunks_preserves_shape(params):
    """Reproducer for https://github.com/dask/dask/issues/3730.

    Mutating based on an array with different chunks can cause new chunks to be
    used.  We need to ensure those new chunk sizes are applied to the mutated
    array, otherwise the array won't generate the correct keys.
    """
    array_size, chunk_size1, chunk_size2 = params
    x = da.zeros(array_size, chunks=chunk_size1)
    mask = da.zeros(array_size, chunks=chunk_size2)
    x[mask] = 1
    result = x.compute()
    assert x.shape == result.shape


def test_gh3579():
    assert_eq(np.arange(10)[0::-1], da.arange(10, chunks=3)[0::-1])
    assert_eq(np.arange(10)[::-1], da.arange(10, chunks=3)[::-1])


def test_make_blockwise_sorted_slice():
    x = da.arange(8, chunks=4)
    index = np.array([6, 0, 4, 2, 7, 1, 5, 3])

    a, b = make_block_sorted_slices(index, x.chunks)

    index2 = np.array([0, 2, 4, 6, 1, 3, 5, 7])
    index3 = np.array([3, 0, 2, 1, 7, 4, 6, 5])
    np.testing.assert_array_equal(a, index2)
    np.testing.assert_array_equal(b, index3)


@pytest.mark.filterwarnings("ignore:Slicing:dask.array.core.PerformanceWarning")
@pytest.mark.parametrize(
    "size, chunks", [((100, 2), (50, 2)), ((100, 2), (37, 1)), ((100,), (55,))]
)
def test_shuffle_slice(size, chunks):
    x = da.random.default_rng().integers(0, 1000, size=size, chunks=chunks)
    index = np.arange(len(x))
    np.random.default_rng().shuffle(index)

    a = x[index]
    b = shuffle_slice(x, index)
    assert_eq(a, b)

    index = np.arange(1, len(x)).tolist()
    index.append(0)
    index = np.array(index)
    a = x[index]
    b = shuffle_slice(x, index)
    assert_eq(a, b)


def test_unknown_chunks_length_one():
    a = np.arange(256, dtype=int)
    arr = da.from_array(a, chunks=(256,))
    # np.flatnonzero dispatches
    result = np.flatnonzero(arr)
    assert_eq(result[[0, -1]], np.flatnonzero(a)[[0, -1]])

    result = da.flatnonzero(arr)
    assert_eq(result[[0, -1]], np.flatnonzero(a)[[0, -1]])

    a = a.reshape(16, 16)
    arr = da.from_array(a, chunks=(8, 16))
    arr._chunks = ((8, 8), (np.nan,))
    result = arr[:, [0, -1]]
    expected = a[:, [0, -1]]
    assert_eq(result, expected)

    arr = da.from_array(a, chunks=(8, 8))
    arr._chunks = ((8, 8), (np.nan, np.nan))
    with pytest.raises(ValueError, match="Array chunk size or shape"):
        arr[:, [0, -1]]


@pytest.mark.parametrize("lock", [True, False])
@pytest.mark.parametrize("asarray", [True, False])
@pytest.mark.parametrize("fancy", [True, False])
def test_gh4043(lock, asarray, fancy):
    a1 = da.from_array(np.zeros(3), chunks=1, asarray=asarray, lock=lock, fancy=fancy)
    a2 = da.from_array(np.ones(3), chunks=1, asarray=asarray, lock=lock, fancy=fancy)
    al = da.stack([a1, a2])
    assert_eq(al, al)


def test_slice_array_3d_with_bool_numpy_array():
    # https://github.com/dask/dask/issues/6089
    array = da.arange(0, 24).reshape((4, 3, 2))
    mask = np.arange(0, 24).reshape((4, 3, 2)) > 12

    actual = array[mask].compute()
    expected = np.arange(13, 24)
    assert_eq(actual, expected)


def test_slice_masked_arrays():
    arr = np.ma.array(range(8), mask=[0, 0, 1, 0, 0, 1, 0, 1])
    darr = da.from_array(arr, chunks=(4, 4))
    assert_eq(darr[[2, 6]], arr[[2, 6]])


def test_slice_array_null_dimension():
    array = da.from_array(np.zeros((3, 0)))
    expected = np.zeros((3, 0))[[0]]
    assert_eq(array[[0]], expected)


def test_take_sorted_indexer():
    arr = da.ones((250, 100), chunks=((50, 100, 33, 67), 100))
    indexer = list(range(0, 250))
    result = arr[indexer, :]
    assert_eq(arr, result)
    assert {
        **dict(arr.dask),
        **{
            k: Alias(k, k2)
            for k, k2 in zip(
                [k for k in dict(result.dask) if "getitem" in k[0]],
                dict(arr.dask).keys(),
            )
        },
    } == dict(result.dask)


def test_all_none_slices_just_mappings():
    arr = da.ones((10, 10), chunks=(1, 5))
    result = arr[slice(None, 6), slice(None)]
    dsk = dict(result.dask)
    assert len([k for k in dsk if "getitem" in k[0]]) == 12
    # check that we are just mapping the keys
    assert all(isinstance(v, Alias) for k, v in dsk.items() if "getitem" in k[0])
    assert_eq(result, np.ones((6, 10)))


def test_minimal_dtype_doesnt_overflow():
    x = np.arange(1980)
    dx = dask.array.from_array(x, chunks=[248])
    ib = np.zeros(1980, dtype=bool)
    ib[1560:1860] = True
    assert_eq(dx[ib], x[ib])


def test_vindex_with_dask_array():
    arr = np.array([0.2, 0.4, 0.6])
    darr = da.from_array(arr, chunks=-1)

    indexer = np.random.randint(0, 3, 8).reshape(4, 2).astype(int)
    dindexer = da.from_array(indexer, chunks=(2, 2))
    assert_eq(darr.vindex[dindexer], arr[indexer])

    msg = "vindex does not support indexing"

    with pytest.raises(IndexError, match=msg):
        darr.rechunk((1, 1)).vindex[dindexer]

    with pytest.raises(IndexError, match=msg):
        darr.reshape((3, 1)).vindex[dindexer]

    with pytest.raises(IndexError, match=msg):
        darr.vindex[(dindexer, None)]


def test_positional_indexer_newaxis():
    arr = da.array([0, 1, 2])
    new = arr[[True, True, False], np.newaxis]
    assert_eq(new, arr.compute()[[True, True, False], np.newaxis])


@pytest.mark.parametrize(
    "shapes",
    [
        (10, 10),
        (np.nan, np.nan),
        pytest.param(
            (10, np.nan), marks=pytest.mark.xfail(reason="Not implemented", strict=True)
        ),
        pytest.param(
            (np.nan, 10), marks=pytest.mark.xfail(reason="Not implemented", strict=True)
        ),
    ],
)
def test_boolean_mask_with_unknown_shape(
    shapes: tuple[float | int, float | int],
) -> None:
    x_shape, mask_shape = shapes
    arr = delayed(np.ones(10))
    x = da.concatenate(
        [
            da.from_delayed(arr, shape=(x_shape,), dtype=float),
            da.from_delayed(arr, shape=(x_shape,), dtype=float),
        ]
    )
    mask = da.concatenate(
        [
            da.from_delayed(arr, shape=(mask_shape,), dtype=bool),
            da.from_delayed(arr, shape=(mask_shape,), dtype=bool),
        ]
    )
    x[mask] = 2

    expected = np.full(20, 2.0)
    assert_eq(x, expected)
