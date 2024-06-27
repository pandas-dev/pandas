from __future__ import annotations

import pytest

pytest.importorskip("numpy")

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

import dask.array as da
from dask.array.lib.stride_tricks import sliding_window_view
from dask.array.overlap import (
    boundaries,
    constant,
    ensure_minimum_chunksize,
    nearest,
    overlap,
    overlap_internal,
    periodic,
    reflect,
    trim_internal,
)
from dask.array.utils import assert_eq, same_keys


def test_overlap_internal():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    g = overlap_internal(d, {0: 2, 1: 1})
    result = g.compute(scheduler="sync")
    assert g.chunks == ((6, 6), (5, 5))

    expected = np.array(
        [
            [0, 1, 2, 3, 4, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 12, 11, 12, 13, 14, 15],
            [16, 17, 18, 19, 20, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 28, 27, 28, 29, 30, 31],
            [32, 33, 34, 35, 36, 35, 36, 37, 38, 39],
            [40, 41, 42, 43, 44, 43, 44, 45, 46, 47],
            [16, 17, 18, 19, 20, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 28, 27, 28, 29, 30, 31],
            [32, 33, 34, 35, 36, 35, 36, 37, 38, 39],
            [40, 41, 42, 43, 44, 43, 44, 45, 46, 47],
            [48, 49, 50, 51, 52, 51, 52, 53, 54, 55],
            [56, 57, 58, 59, 60, 59, 60, 61, 62, 63],
        ]
    )

    assert_eq(result, expected)
    assert same_keys(overlap_internal(d, {0: 2, 1: 1}), g)


def test_overlap_internal_asymmetric():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    result = overlap_internal(d, {0: (2, 0), 1: (1, 0)})
    assert result.chunks == ((4, 6), (4, 5))

    expected = np.array(
        [
            [0, 1, 2, 3, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 11, 12, 13, 14, 15],
            [16, 17, 18, 19, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 27, 28, 29, 30, 31],
            [16, 17, 18, 19, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 27, 28, 29, 30, 31],
            [32, 33, 34, 35, 35, 36, 37, 38, 39],
            [40, 41, 42, 43, 43, 44, 45, 46, 47],
            [48, 49, 50, 51, 51, 52, 53, 54, 55],
            [56, 57, 58, 59, 59, 60, 61, 62, 63],
        ]
    )
    assert_eq(result, expected)
    assert same_keys(overlap_internal(d, {0: (2, 0), 1: (1, 0)}), result)


def test_overlap_internal_asymmetric_small():
    x = np.arange(32).reshape((2, 16))
    d = da.from_array(x, chunks=(2, 4))

    result = overlap_internal(d, {0: (0, 0), 1: (1, 1)})
    assert result.chunks == ((2,), (5, 6, 6, 5))

    expected = np.array(
        [
            [0, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 7, 8, 9, 10, 11, 12, 11, 12, 13, 14, 15],
            [
                16,
                17,
                18,
                19,
                20,
                19,
                20,
                21,
                22,
                23,
                24,
                23,
                24,
                25,
                26,
                27,
                28,
                27,
                28,
                29,
                30,
                31,
            ],
        ]
    )

    assert_eq(result, expected)
    assert same_keys(overlap_internal(d, {0: (0, 0), 1: (1, 1)}), result)


def test_trim_internal():
    d = da.ones((40, 60), chunks=(10, 10))
    e = trim_internal(d, axes={0: 1, 1: 2}, boundary="reflect")

    assert e.chunks == ((8, 8, 8, 8), (6, 6, 6, 6, 6, 6))


def test_periodic():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = periodic(d, axis=0, depth=2)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], d[-1, :])
    assert_eq(e[0, :], d[-2, :])


def test_reflect():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    e = reflect(d, axis=0, depth=2)
    expected = np.array([1, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8])
    assert_eq(e, expected)

    e = reflect(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected)


def test_nearest():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    e = nearest(d, axis=0, depth=2)
    expected = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9])
    assert_eq(e, expected)

    e = nearest(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected)


def test_constant():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = constant(d, axis=0, depth=2, value=10)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], np.ones(8, dtype=x.dtype) * 10)
    assert_eq(e[-1, :], np.ones(8, dtype=x.dtype) * 10)


def test_boundaries():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = boundaries(d, {0: 2, 1: 1}, {0: 0, 1: "periodic"})

    expected = np.array(
        [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [7, 0, 1, 2, 3, 4, 5, 6, 7, 0],
            [15, 8, 9, 10, 11, 12, 13, 14, 15, 8],
            [23, 16, 17, 18, 19, 20, 21, 22, 23, 16],
            [31, 24, 25, 26, 27, 28, 29, 30, 31, 24],
            [39, 32, 33, 34, 35, 36, 37, 38, 39, 32],
            [47, 40, 41, 42, 43, 44, 45, 46, 47, 40],
            [55, 48, 49, 50, 51, 52, 53, 54, 55, 48],
            [63, 56, 57, 58, 59, 60, 61, 62, 63, 56],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    )
    assert_eq(e, expected)


def test_overlap():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))
    g = overlap(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: "reflect"})
    assert g.chunks == ((8, 8), (6, 6))
    expected = np.array(
        [
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [0, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7, 7],
            [8, 8, 9, 10, 11, 12, 11, 12, 13, 14, 15, 15],
            [16, 16, 17, 18, 19, 20, 19, 20, 21, 22, 23, 23],
            [24, 24, 25, 26, 27, 28, 27, 28, 29, 30, 31, 31],
            [32, 32, 33, 34, 35, 36, 35, 36, 37, 38, 39, 39],
            [40, 40, 41, 42, 43, 44, 43, 44, 45, 46, 47, 47],
            [16, 16, 17, 18, 19, 20, 19, 20, 21, 22, 23, 23],
            [24, 24, 25, 26, 27, 28, 27, 28, 29, 30, 31, 31],
            [32, 32, 33, 34, 35, 36, 35, 36, 37, 38, 39, 39],
            [40, 40, 41, 42, 43, 44, 43, 44, 45, 46, 47, 47],
            [48, 48, 49, 50, 51, 52, 51, 52, 53, 54, 55, 55],
            [56, 56, 57, 58, 59, 60, 59, 60, 61, 62, 63, 63],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
        ]
    )
    assert_eq(g, expected)
    assert same_keys(g, overlap(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: "reflect"}))

    u_depth = np.uint16([2, 1])
    u_depth = {k: v for k, v in enumerate(u_depth)}
    g = overlap(d, depth=u_depth, boundary={0: 100, 1: "reflect"})
    assert g.chunks == ((8, 8), (6, 6))
    assert_eq(g, expected)
    assert same_keys(g, overlap(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: "reflect"}))

    g = overlap(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: "none"})
    expected = np.array(
        [
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [0, 1, 2, 3, 4, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 12, 11, 12, 13, 14, 15],
            [16, 17, 18, 19, 20, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 28, 27, 28, 29, 30, 31],
            [32, 33, 34, 35, 36, 35, 36, 37, 38, 39],
            [40, 41, 42, 43, 44, 43, 44, 45, 46, 47],
            [16, 17, 18, 19, 20, 19, 20, 21, 22, 23],
            [24, 25, 26, 27, 28, 27, 28, 29, 30, 31],
            [32, 33, 34, 35, 36, 35, 36, 37, 38, 39],
            [40, 41, 42, 43, 44, 43, 44, 45, 46, 47],
            [48, 49, 50, 51, 52, 51, 52, 53, 54, 55],
            [56, 57, 58, 59, 60, 59, 60, 61, 62, 63],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
        ]
    )
    assert_eq(g, expected)
    assert g.chunks == ((8, 8), (5, 5))

    u_depth = np.uint16([2, 1])
    u_depth = {k: v for k, v in enumerate(u_depth)}
    g = overlap(d, depth=u_depth, boundary={0: 100, 1: "none"})
    assert_eq(g, expected)
    assert g.chunks == ((8, 8), (5, 5))


def test_overlap_allow_rechunk_kwarg():
    # The smallest array chunk is too small to fit overlap depth
    arr = da.arange(6, chunks=5)
    da.overlap.overlap(arr, 2, "reflect", allow_rechunk=True)
    arr.map_overlap(lambda x: x, 2, "reflect", allow_rechunk=True)
    with pytest.raises(ValueError):
        da.overlap.overlap(arr, 2, "reflect", allow_rechunk=False)
    with pytest.raises(ValueError):
        arr.map_overlap(lambda x: x, 2, "reflect", allow_rechunk=False)
    # No rechunking required
    arr = da.arange(6, chunks=4)
    da.overlap.overlap(arr, 2, "reflect", allow_rechunk=False)


def test_asymmetric_overlap_boundary_exception():
    x = da.arange(10, chunks=5)
    with pytest.raises(NotImplementedError):
        x.map_overlap(
            lambda x: x + len(x), depth={0: (0, 2)}, boundary="reflect", dtype=x.dtype
        )


def test_map_overlap():
    x = da.arange(10, chunks=5)
    y = x.map_overlap(lambda x: x + len(x), depth=2, dtype=x.dtype, boundary="reflect")
    assert_eq(y, np.arange(10) + 5 + 2 + 2)

    x = da.arange(10, chunks=5)
    y = x.map_overlap(
        lambda x: x + len(x), depth=np.int64(2), dtype=x.dtype, boundary="reflect"
    )
    assert all([(type(s) is int) for s in y.shape])
    assert_eq(y, np.arange(10) + 5 + 2 + 2)

    x = da.ones((10, 10), chunks=(3, 4))
    z = x.map_overlap(lambda x: x, depth={1: 5}, boundary="reflect")
    assert z.chunks[0] == x.chunks[0]  # don't rechunk the first dimension

    x = np.arange(16).reshape((4, 4))
    d = da.from_array(x, chunks=(2, 2))
    exp1 = d.map_overlap(
        lambda x: x + x.size, depth=1, dtype=d.dtype, boundary="reflect"
    )
    exp2 = d.map_overlap(
        lambda x: x + x.size,
        depth={0: 1, 1: 1},
        boundary={0: "reflect", 1: "none"},
        dtype=d.dtype,
    )
    exp3 = d.map_overlap(
        lambda x: x + x.size, depth={1: 1}, boundary={1: "reflect"}, dtype=d.dtype
    )
    exp4 = d.map_overlap(
        lambda x: x + x.size,
        depth={1: (1, 0)},
        boundary={0: "none", 1: "none"},
        dtype=d.dtype,
    )
    assert_eq(exp1, x + 16)
    assert_eq(exp2, x + 12)
    assert_eq(exp3, x + 8)
    assert_eq(
        exp4,
        np.block(
            [[x[0:2, 0:2] + 4, x[0:2, 2:4] + 6], [x[2:4, 0:2] + 4, x[2:4, 2:4] + 6]]
        ),
    )


def test_map_overlap_escapes_to_map_blocks_when_depth_is_zero():
    x = da.arange(10, chunks=5)
    y = x.map_overlap(lambda x: x + 1, depth=0, boundary="none")
    assert len(y.dask) == 2 * x.numblocks[0]  # depth=0 --> map_blocks
    assert_eq(y, np.arange(10) + 1)


@pytest.mark.parametrize(
    "boundary", [None, "reflect", "periodic", "nearest", "none", 0]
)
def test_map_overlap_no_depth(boundary):
    x = da.arange(10, chunks=5)
    y = x.map_overlap(lambda i: i, depth=0, boundary=boundary, dtype=x.dtype)
    assert_eq(y, x)


def test_map_overlap_multiarray():
    # Same ndim, same numblocks, same chunks
    x = da.arange(10, chunks=5)
    y = da.arange(10, chunks=5)
    z = da.map_overlap(lambda x, y: x + y, x, y, depth=1, boundary="none")
    assert_eq(z, 2 * np.arange(10))

    # Same ndim, same numblocks, different chunks
    x = da.arange(10, chunks=(2, 3, 5))
    y = da.arange(10, chunks=(5, 3, 2))
    z = da.map_overlap(lambda x, y: x + y, x, y, depth=1, boundary="none")
    assert z.chunks == ((2, 3, 3, 2),)
    assert_eq(z, 2 * np.arange(10))

    # Same ndim, different numblocks, different chunks
    x = da.arange(10, chunks=(10,))
    y = da.arange(10, chunks=(4, 4, 2))
    z = da.map_overlap(lambda x, y: x + y, x, y, depth=1, boundary="none")
    assert z.chunks == ((4, 4, 2),)
    assert_eq(z, 2 * np.arange(10))

    # Different ndim, different numblocks, different chunks
    x = da.arange(10, chunks=(10,))
    y = da.arange(10).reshape(1, 10).rechunk((1, (4, 4, 2)))
    z = da.map_overlap(lambda x, y: x + y, x, y, depth=1, boundary="none")
    assert z.chunks == ((1,), (4, 4, 2))
    assert z.shape == (1, 10)
    assert_eq(z, 2 * np.arange(10)[np.newaxis])

    # Note: checks on arange equality in all of the above help ensure that
    # trimming is applied appropriately to result chunks (i.e. results
    # are not somehow shifted)


def test_map_overlap_multiarray_defaults():
    # Check that by default, chunk alignment and arrays of varying dimensionality
    # are supported by with no effect on result shape
    # (i.e. defaults are pass-through to map_blocks)
    x = da.ones((10,), chunks=10)
    y = da.ones((1, 10), chunks=5)
    z = da.map_overlap(lambda x, y: x + y, x, y, boundary="none")
    # func should be called twice and get (5,) and (1, 5) arrays of ones each time
    assert_eq(z.shape, (1, 10))
    assert_eq(z.sum(), 20.0)


def test_map_overlap_multiarray_different_depths():
    x = da.ones(5, dtype="int")
    y = da.ones(5, dtype="int")

    def run(depth):
        return da.map_overlap(
            lambda x, y: x.sum() + y.sum(),
            x,
            y,
            depth=depth,
            chunks=(0,),
            trim=False,
            boundary="reflect",
        ).compute()

    # Check that the number of elements added
    # to arrays in overlap works as expected
    # when depths differ for each array
    assert run([0, 0]) == 10
    assert run([0, 1]) == 12
    assert run([1, 1]) == 14
    assert run([1, 2]) == 16
    assert run([0, 5]) == 20
    assert run([5, 5]) == 30

    # Ensure that depth > chunk size results in error
    with pytest.raises(ValueError):
        run([0, 6])


def test_map_overlap_multiarray_uneven_numblocks_exception():
    x = da.arange(10, chunks=(10,))
    y = da.arange(10, chunks=(5, 5))
    with pytest.raises(ValueError):
        # Fail with chunk alignment explicitly disabled
        da.map_overlap(
            lambda x, y: x + y, x, y, align_arrays=False, boundary="none"
        ).compute()


def test_map_overlap_multiarray_block_broadcast():
    def func(x, y):
        # Return result with expected padding
        z = x.size + y.size
        return np.ones((3, 3)) * z

    # Chunks in trailing dimension will be unified to two chunks of size 6
    # and block broadcast will allow chunks from x to repeat
    x = da.ones((12,), chunks=12)  # numblocks = (1,) -> (2, 2) after broadcast
    y = da.ones((16, 12), chunks=(8, 6))  # numblocks = (2, 2)
    z = da.map_overlap(
        func, x, y, chunks=(3, 3), depth=1, trim=True, boundary="reflect"
    )
    assert_eq(z, z)
    assert z.shape == (2, 2)
    # func call will receive (8,) and (10, 8) arrays for each of 4 blocks
    assert_eq(z.sum(), 4.0 * (10 * 8 + 8))


def test_map_overlap_multiarray_variadic():
    # Test overlapping row slices from 3D arrays
    xs = [
        # Dim 0 will unify to chunks of size 4 for all:
        da.ones((12, 1, 1), chunks=((12,), 1, 1)),
        da.ones((12, 8, 1), chunks=((8, 4), 8, 1)),
        da.ones((12, 8, 4), chunks=((4, 8), 8, 4)),
    ]

    def func(*args):
        return np.array([sum(x.size for x in args)])

    x = da.map_overlap(
        func,
        *xs,
        chunks=(1,),
        depth=1,
        trim=False,
        drop_axis=[1, 2],
        boundary="reflect",
    )

    # Each func call should get 4 rows from each array padded by 1 in each dimension
    size_per_slice = sum(np.pad(x[:4], 1, mode="constant").size for x in xs)
    assert x.shape == (3,)
    assert all(x.compute() == size_per_slice)


@pytest.mark.parametrize(
    "drop_axis",
    (
        (0,),
        (1,),
        (2,),
        (0, 1),
        (1, 2),
        (2, 0),
        1,
        (-3,),
        (-2,),
        (-1,),
        (-3, -2),
        (-2, -1),
        (-1, -3),
        -2,
    ),
)
def test_map_overlap_trim_using_drop_axis_and_different_depths(drop_axis):
    x = da.random.default_rng().standard_normal((5, 10, 8), chunks=(2, 5, 4))

    def _mean(x):
        return x.mean(axis=drop_axis)

    expected = _mean(x)

    # unique boundary and depth value per axis
    boundary = (0, "reflect", "nearest")
    depth = (1, 3, 2)
    # to match expected result, dropped axes must have depth 0
    _drop_axis = (drop_axis,) if np.isscalar(drop_axis) else drop_axis
    _drop_axis = [d % x.ndim for d in _drop_axis]
    depth = tuple(0 if i in _drop_axis else d for i, d in enumerate(depth))

    y = da.map_overlap(
        _mean, x, depth=depth, boundary=boundary, drop_axis=drop_axis, dtype=float
    ).compute()
    assert_array_almost_equal(expected, y)


def test_map_overlap_assumes_shape_matches_first_array_if_trim_is_false():
    # https://github.com/dask/dask/issues/6681
    x1 = da.ones((10,), chunks=(5, 5))
    x2 = x1.rechunk(10)

    def oversum(x):
        return x[2:-2]

    z1 = da.map_overlap(oversum, x1, depth=2, trim=False, boundary="none")
    assert z1.shape == (10,)

    z2 = da.map_overlap(oversum, x2, depth=2, trim=False, boundary="none")
    assert z2.shape == (10,)


def test_map_overlap_deprecated_signature():
    def func(x):
        return np.array(x.sum())

    x = da.ones(3)

    # Old positional signature: func, depth, boundary, trim
    with pytest.warns(FutureWarning):
        y = da.map_overlap(x, func, 0, "reflect", True)
        assert y.compute() == 3
        assert y.shape == (3,)

    with pytest.warns(FutureWarning):
        y = da.map_overlap(x, func, 1, "reflect", True)
        assert y.compute() == 5
        assert y.shape == (3,)

    with pytest.warns(FutureWarning):
        y = da.map_overlap(x, func, 1, "reflect", False)
        assert y.compute() == 5
        assert y.shape == (3,)


def test_nearest_overlap():
    a = np.arange(144).reshape(12, 12).astype(float)

    darr = da.from_array(a, chunks=(6, 6))
    garr = overlap(darr, depth={0: 5, 1: 5}, boundary={0: "nearest", 1: "nearest"})
    tarr = trim_internal(garr, {0: 5, 1: 5}, boundary="nearest")
    assert_array_almost_equal(tarr, a)


@pytest.mark.parametrize(
    "depth",
    [
        {0: 0, 1: 0},  # depth all zeros
        {0: 4, 1: 0},  # depth with some zeros
        {0: 5, 1: 5},  # depth equal to boundary length
        {0: 8, 1: 7},  # depth greater than boundary length
    ],
)
def test_different_depths_and_boundary_combinations(depth):
    expected = np.arange(100).reshape(10, 10)
    darr = da.from_array(expected, chunks=(5, 2))

    reflected = overlap(darr, depth=depth, boundary="reflect")
    nearest = overlap(darr, depth=depth, boundary="nearest")
    periodic = overlap(darr, depth=depth, boundary="periodic")
    constant = overlap(darr, depth=depth, boundary=42)

    result = trim_internal(reflected, depth, boundary="reflect")
    assert_array_equal(result, expected)

    result = trim_internal(nearest, depth, boundary="nearest")
    assert_array_equal(result, expected)

    result = trim_internal(periodic, depth, boundary="periodic")
    assert_array_equal(result, expected)

    result = trim_internal(constant, depth, boundary=42)
    assert_array_equal(result, expected)


def test_one_chunk_along_axis():
    a = np.arange(2 * 9).reshape(2, 9)
    darr = da.from_array(a, chunks=((2,), (2, 2, 2, 3)))
    g = overlap(darr, depth=0, boundary=0)
    assert a.shape == g.shape


def test_constant_boundaries():
    a = np.arange(1 * 9).reshape(1, 9)
    darr = da.from_array(a, chunks=((1,), (2, 2, 2, 3)))
    b = boundaries(darr, {0: 0, 1: 0}, {0: 0, 1: 0})
    assert b.chunks == darr.chunks


@pytest.mark.parametrize(
    "chunks",
    [
        ((5, 5, 2), (5, 5, 2)),
        ((3, 3, 3, 3), (11, 1)),
    ],
)
def test_depth_greater_than_smallest_chunk_combines_chunks(chunks):
    a = np.arange(144).reshape(12, 12)
    darr = da.from_array(a, chunks=chunks)

    depth = {0: 4, 1: 2}
    output = overlap(darr, depth=depth, boundary=1)

    assert all(c >= depth[0] * 2 for c in output.chunks[0])
    assert all(c >= depth[1] * 2 for c in output.chunks[1])


def test_depth_greater_than_dim():
    a = np.arange(144).reshape(12, 12)
    darr = da.from_array(a, chunks=(3, 5))

    depth = {0: 13, 1: 4}
    with pytest.raises(ValueError, match="The overlapping depth"):
        overlap(darr, depth=depth, boundary=1)


def test_none_boundaries():
    x = da.from_array(np.arange(16).reshape(4, 4), chunks=(2, 2))
    exp = boundaries(x, 2, {0: "none", 1: 33})
    res = np.array(
        [
            [33, 33, 0, 1, 2, 3, 33, 33],
            [33, 33, 4, 5, 6, 7, 33, 33],
            [33, 33, 8, 9, 10, 11, 33, 33],
            [33, 33, 12, 13, 14, 15, 33, 33],
        ]
    )
    assert_eq(exp, res)


def test_overlap_small():
    x = da.ones((10, 10), chunks=(5, 5))

    y = x.map_overlap(lambda x: x, depth=1, boundary="none")
    assert len(y.dask) < 200

    y = x.map_overlap(lambda x: x, depth=1, boundary="none")
    assert len(y.dask) < 100


def test_no_shared_keys_with_different_depths():
    rng = da.random.default_rng(0)
    a = rng.random((9, 9), chunks=(3, 3))

    def check(x):
        assert x.shape == (3, 3)
        return x

    r = [
        a.map_overlap(
            lambda a: a + 1,
            dtype=a.dtype,
            depth={j: int(i == j) for j in range(a.ndim)},
            boundary="none",
        ).map_blocks(check, dtype=a.dtype)
        for i in range(a.ndim)
    ]

    assert set(r[0].dask) & set(r[1].dask) == set(a.dask)
    da.compute(*r, scheduler="single-threaded")


def test_overlap_few_dimensions_small():
    x = da.ones((20, 20), chunks=(10, 10))

    a = x.map_overlap(lambda x: x, depth={0: 1}, boundary="none")
    assert_eq(x, a)
    assert any(isinstance(k[1], float) for k in a.dask)
    assert all(isinstance(k[2], int) for k in a.dask)

    b = x.map_overlap(lambda x: x, depth={1: 1}, boundary="none")
    assert_eq(x, b)
    assert all(isinstance(k[1], int) for k in b.dask)
    assert any(isinstance(k[2], float) for k in b.dask)

    c = x.map_overlap(lambda x: x, depth={0: 1, 1: 1}, boundary="none")
    assert_eq(x, c)
    assert any(isinstance(k[1], float) for k in c.dask)
    assert any(isinstance(k[2], float) for k in c.dask)


def test_overlap_few_dimensions():
    x = da.ones((100, 100), chunks=(10, 10))

    a = x.map_overlap(lambda x: x, depth={0: 1}, boundary="none")
    b = x.map_overlap(lambda x: x, depth={1: 1}, boundary="none")
    c = x.map_overlap(lambda x: x, depth={0: 1, 1: 1}, boundary="none")

    assert len(a.dask) == len(b.dask)
    assert len(a.dask) < len(c.dask)

    assert len(c.dask) < 10 * len(a.dask)


@pytest.mark.parametrize("boundary", ["reflect", "periodic", "nearest", "none"])
def test_trim_boundary(boundary):
    x = da.from_array(np.arange(24).reshape(4, 6), chunks=(2, 3))
    x_overlaped = da.overlap.overlap(x, 2, boundary={0: "reflect", 1: boundary})
    x_trimmed = da.overlap.trim_overlap(
        x_overlaped, 2, boundary={0: "reflect", 1: boundary}
    )
    assert np.all(x == x_trimmed)

    x_overlaped = da.overlap.overlap(x, 2, boundary={1: boundary})
    x_trimmed = da.overlap.trim_overlap(x_overlaped, 2, boundary={1: boundary})
    assert np.all(x == x_trimmed)

    x_overlaped = da.overlap.overlap(x, 2, boundary=boundary)
    x_trimmed = da.overlap.trim_overlap(x_overlaped, 2, boundary=boundary)
    assert np.all(x == x_trimmed)


def test_map_overlap_rechunks_array_if_needed():
    # https://github.com/dask/dask/issues/6597
    expected = np.arange(11)
    x = da.from_array(expected, chunks=5)
    y = x.map_overlap(lambda x: x, depth=2, boundary=0)
    assert all(c >= 2 for c in y.chunks[0])
    assert_eq(y, expected)


def test_map_overlap_rechunks_array_along_multiple_dims_if_needed():
    # https://github.com/dask/dask/issues/6688
    rand = da.random.default_rng().random((860, 1024, 1024), chunks=(1, 1024, 1024))
    filtered = rand.map_overlap(
        lambda arr: arr,
        depth=(2, 2, 2),
        boundary="reflect",
    )
    assert all(all(c >= 2 for c in chunks) for chunks in filtered.chunks)


@pytest.mark.parametrize(
    "chunks,expected",
    [
        [(10,), (10,)],
        [(10, 10), (10, 10)],
        [
            (10, 10, 1),
            (10, 11),
        ],
        [(20, 20, 20, 1), (20, 20, 11, 10)],
        [(20, 20, 10, 1), (20, 20, 11)],
        [(2, 20, 2, 20), (14, 10, 20)],
        [(1, 1, 1, 1, 7), (11,)],
        [(20, 20, 2, 20, 20, 2), (20, 12, 10, 20, 12, 10)],
    ],
)
def test_ensure_minimum_chunksize(chunks, expected):
    actual = ensure_minimum_chunksize(10, chunks)
    assert actual == expected


def test_ensure_minimum_chunksize_raises_error():
    chunks = (5, 2, 1, 1)
    with pytest.raises(ValueError, match="overlapping depth 10 is larger than"):
        ensure_minimum_chunksize(10, chunks)


@pytest.mark.parametrize(
    "shape, chunks, window_shape, axis",
    [
        ((6, 7, 8), (6, (2, 2, 2, 1), 4), (3, 2), (1, 2)),  # chunks vary along axis
        ((40, 30, 2), 5, (3,), (0,)),  # window < chunk
        ((21,), 3, (7,), (0,)),  # window > chunk
        ((9,), 3, 3, 0),  # window == chunk, axis is integer
        ((9,), 3, 3, -1),  # axis=-1
        ((9,), 3, 3, None),  # axis=None
        ((9, 8), 3, (2, 4), None),  # axis=None
        ((9,), 3, (3, 3, 3), (0, 0, 0)),  # axis is repeated
        ((9,), 3, (3, 3), (0, -1)),  # axis is repeated, with -1
        ((9,), 3, [3, 3], [0, -1]),  # list instead of tuple
    ],
)
def test_sliding_window_view(shape, chunks, window_shape, axis):
    arr = da.from_array(np.arange(np.prod(shape)).reshape(shape), chunks=chunks)
    actual = sliding_window_view(arr, window_shape, axis)
    expected = np.lib.stride_tricks.sliding_window_view(
        arr.compute(), window_shape, axis
    )
    assert_eq(expected, actual)


@pytest.mark.parametrize(
    "window_shape, axis",
    [
        ((10,), 0),  # window > axis shape
        ((2,), 3),  # axis > ndim
        (-1, 0),  # window shape is negative
        (2, (0, 1)),  # len(window shape) < len(axis)
        (2, None),  # len(window shape) < len(axis)
        (0, None),  # window_shape = 0
    ],
)
def test_sliding_window_errors(window_shape, axis):
    arr = da.zeros((4, 3))
    with pytest.raises(ValueError):
        sliding_window_view(arr, window_shape, axis)
