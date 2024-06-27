from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.utils import assert_eq, same_keys

cupy = pytest.importorskip("cupy")
cupy_version = Version(cupy.__version__)


@pytest.mark.skipif(
    cupy_version < Version("6.4.0"),
    reason="Requires CuPy 6.4.0+ (with https://github.com/cupy/cupy/pull/2418)",
)
def test_bincount():
    x = cupy.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2, asarray=False)
    e = da.bincount(d, minlength=6)
    assert_eq(e, np.bincount(x, minlength=6))
    assert same_keys(da.bincount(d, minlength=6), e)

    assert da.bincount(d, minlength=6).name != da.bincount(d, minlength=7).name
    assert da.bincount(d, minlength=6).name == da.bincount(d, minlength=6).name


def test_compress():
    carr = cupy.random.default_rng().integers(0, 3, size=(10, 10))

    darr = da.from_array(carr, chunks=(20, 5))

    c = cupy.asarray([True])
    res = da.compress(c, darr, axis=0)

    # cupy.compress is not implemented but dask implementation does not
    # rely on np.compress -- move original data back to host and
    # compare da.compress with np.compress
    assert_eq(np.compress(c.tolist(), carr.tolist(), axis=0), res, check_type=False)


@pytest.mark.parametrize(
    "shape, axis",
    [[(10, 15, 20), 0], [(10, 15, 20), 1], [(10, 15, 20), 2], [(10, 15, 20), -1]],
)
@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff(shape, n, axis):
    x = cupy.random.default_rng().integers(0, 10, shape)
    a = da.from_array(x, chunks=(len(shape) * (5,)))

    assert_eq(da.diff(a, n, axis), cupy.diff(x, n, axis))


@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff_prepend(n):
    x = cupy.arange(5) + 1
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, prepend=0), cupy.diff(x, n, prepend=0))
    assert_eq(da.diff(a, n, prepend=[0]), cupy.diff(x, n, prepend=[0]))
    assert_eq(da.diff(a, n, prepend=[-1, 0]), cupy.diff(x, n, prepend=[-1, 0]))

    x = cupy.arange(16).reshape(4, 4)
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, axis=1, prepend=0), cupy.diff(x, n, axis=1, prepend=0))
    assert_eq(
        da.diff(a, n, axis=1, prepend=[[0], [0], [0], [0]]),
        cupy.diff(x, n, axis=1, prepend=[[0], [0], [0], [0]]),
    )
    assert_eq(da.diff(a, n, axis=0, prepend=0), cupy.diff(x, n, axis=0, prepend=0))
    assert_eq(
        da.diff(a, n, axis=0, prepend=[[0, 0, 0, 0]]),
        cupy.diff(x, n, axis=0, prepend=[[0, 0, 0, 0]]),
    )

    if n > 0:
        # When order is 0 the result is the icupyut array, it doesn't raise
        # an error
        with pytest.raises(ValueError):
            da.diff(a, n, prepend=cupy.zeros((3, 3)))


@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff_append(n):
    x = cupy.arange(5) + 1
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, append=0), cupy.diff(x, n, append=0))
    assert_eq(da.diff(a, n, append=[0]), cupy.diff(x, n, append=[0]))
    assert_eq(da.diff(a, n, append=[-1, 0]), cupy.diff(x, n, append=[-1, 0]))

    x = cupy.arange(16).reshape(4, 4)
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, axis=1, append=0), cupy.diff(x, n, axis=1, append=0))
    assert_eq(
        da.diff(a, n, axis=1, append=[[0], [0], [0], [0]]),
        cupy.diff(x, n, axis=1, append=[[0], [0], [0], [0]]),
    )
    assert_eq(da.diff(a, n, axis=0, append=0), cupy.diff(x, n, axis=0, append=0))
    assert_eq(
        da.diff(a, n, axis=0, append=[[0, 0, 0, 0]]),
        cupy.diff(x, n, axis=0, append=[[0, 0, 0, 0]]),
    )

    if n > 0:
        with pytest.raises(ValueError):
            # When order is 0 the result is the icupyut array, it doesn't raise
            # an error
            da.diff(a, n, append=cupy.zeros((3, 3)))


@pytest.mark.parametrize("bins_type", [np, cupy])
def test_digitize(bins_type):
    x = cupy.array([2, 4, 5, 6, 1])
    bins = bins_type.array([1, 2, 3, 4, 5])
    for chunks in [2, 4]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            bins_cupy = cupy.array(bins)
            assert_eq(
                da.digitize(d, bins, right=right),
                np.digitize(x, bins_cupy, right=right),
                check_type=False,
            )

    x = cupy.random.default_rng().random(size=(100, 100))
    bins = bins_type.random.default_rng().random(size=13)
    bins.sort()
    for chunks in [(10, 10), (10, 20), (13, 17), (87, 54)]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            bins_cupy = cupy.array(bins)
            assert_eq(
                da.digitize(d, bins, right=right),
                np.digitize(x, bins_cupy, right=right),
            )


@pytest.mark.skipif(
    cupy_version < Version("6.4.0"),
    reason="Requires CuPy 6.4.0+ (with https://github.com/cupy/cupy/pull/2418)",
)
def test_tril_triu():
    A = cupy.random.default_rng().standard_normal((20, 20))
    for chk in [5, 4]:
        dA = da.from_array(A, (chk, chk), asarray=False)

        assert_eq(da.triu(dA), np.triu(A))
        assert_eq(da.tril(dA), np.tril(A))

        for k in [-25, -20, -9, -1, 1, 8, 19, 21]:
            assert_eq(da.triu(dA, k), np.triu(A, k))
            assert_eq(da.tril(dA, k), np.tril(A, k))


@pytest.mark.skipif(
    cupy_version < Version("6.4.0"),
    reason="Requires CuPy 6.4.0+ (with https://github.com/cupy/cupy/pull/2418)",
)
def test_tril_triu_non_square_arrays():
    A = cupy.random.default_rng().integers(0, 11, (30, 35))
    dA = da.from_array(A, chunks=(5, 5), asarray=False)
    assert_eq(da.triu(dA), np.triu(A))
    assert_eq(da.tril(dA), np.tril(A))


@pytest.mark.parametrize("return_index", [False, True])
@pytest.mark.parametrize("return_inverse", [False, True])
@pytest.mark.parametrize("return_counts", [False, True])
def test_unique_kwargs(return_index, return_inverse, return_counts):
    kwargs = dict(
        return_index=return_index,
        return_inverse=return_inverse,
        return_counts=return_counts,
    )

    a = cupy.array([1, 2, 4, 4, 5, 2])
    d = da.from_array(a, chunks=(3,))

    def _test_unique_kwargs():
        r_a = np.unique(a, **kwargs)
        r_d = da.unique(d, **kwargs)

        if not any([return_index, return_inverse, return_counts]):
            assert isinstance(r_a, cupy.ndarray)
            assert isinstance(r_d, da.Array)

            r_a = (r_a,)
            r_d = (r_d,)

        assert len(r_a) == len(r_d)

        if return_inverse:
            i = 1 + int(return_index)
            assert (d.size,) == r_d[i].shape

        for e_r_a, e_r_d in zip(r_a, r_d):
            assert_eq(e_r_d, e_r_a)

    # `return_index`, `return_inverse` and `return_counts` are currently
    # unsupported on CuPy-backed Dask arrays.
    if any(kwargs.values()):
        with pytest.raises(ValueError):
            _test_unique_kwargs()


@pytest.mark.parametrize("seed", [23, 796])
@pytest.mark.parametrize("low, high", [[0, 10]])
@pytest.mark.parametrize(
    "shape, chunks",
    [[(10,), (5,)], [(10,), (3,)], [(4, 5), (3, 2)], [(20, 20), (4, 5)]],
)
def test_unique_rand(seed, low, high, shape, chunks):
    rng = cupy.random.default_rng(seed)

    a = rng.integers(low, high, size=shape)
    d = da.from_array(a, chunks=chunks)

    r_a = np.unique(a)
    r_d = da.unique(d)
    assert_eq(r_d, r_a)
