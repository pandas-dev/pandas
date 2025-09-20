from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.utils import assert_eq, same_keys

cupy = pytest.importorskip("cupy")
cupy_version = Version(cupy.__version__)


def test_overlap_internal():
    x = cupy.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4), asarray=False)

    g = da.overlap.overlap_internal(d, {0: 2, 1: 1})
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

    assert_eq(g, expected, check_type=False)
    assert same_keys(da.overlap.overlap_internal(d, {0: 2, 1: 1}), g)


def test_trim_internal():
    x = cupy.ones((40, 60))
    d = da.from_array(x, chunks=(10, 10), asarray=False)
    e = da.overlap.trim_internal(d, axes={0: 1, 1: 2}, boundary="reflect")

    assert e.chunks == ((8, 8, 8, 8), (6, 6, 6, 6, 6, 6))


def test_periodic():
    x = cupy.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4), asarray=False)

    e = da.overlap.periodic(d, axis=0, depth=2)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], d[-1, :])
    assert_eq(e[0, :], d[-2, :])


def test_reflect():
    x = cupy.arange(10)
    d = da.from_array(x, chunks=(5, 5), asarray=False)

    e = da.overlap.reflect(d, axis=0, depth=2)
    expected = np.array([1, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8])
    assert_eq(e, expected, check_type=False)

    e = da.overlap.reflect(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected, check_type=False)


def test_nearest():
    x = cupy.arange(10)
    d = da.from_array(x, chunks=(5, 5), asarray=False)

    e = da.overlap.nearest(d, axis=0, depth=2)
    expected = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9])
    assert_eq(e, expected, check_type=False)

    e = da.overlap.nearest(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected, check_type=False)


@pytest.mark.skipif(
    cupy_version < Version("6.4.0"),
    reason="Requires CuPy 6.4.0+ (with https://github.com/cupy/cupy/pull/2418)",
)
def test_constant():
    x = cupy.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4), asarray=False)

    e = da.overlap.constant(d, axis=0, depth=2, value=10)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], np.ones(8, dtype=x.dtype) * 10, check_type=False)
    assert_eq(e[-1, :], np.ones(8, dtype=x.dtype) * 10, check_type=False)


@pytest.mark.skipif(
    cupy_version < Version("6.4.0"),
    reason="Requires CuPy 6.4.0+ (with https://github.com/cupy/cupy/pull/2418)",
)
def test_boundaries():
    x = cupy.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4), asarray=False)

    e = da.overlap.boundaries(d, {0: 2, 1: 1}, {0: 0, 1: "periodic"})

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
    assert_eq(e, expected, check_type=False)
