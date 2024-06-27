from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

pytestmark = pytest.mark.gpu

import dask
import dask.array as da
from dask.array.utils import assert_eq
from dask.sizeof import sizeof

cupy = pytest.importorskip("cupy")
cupy_version = Version(cupy.__version__)


functions = [
    lambda x: x,
    lambda x: da.expm1(x),
    lambda x: 2 * x,
    lambda x: x / 2,
    lambda x: x**2,
    lambda x: x + x,
    lambda x: x * x,
    lambda x: x[0],
    lambda x: x[:, 1],
    lambda x: x[:1, None, 1:3],
    lambda x: x.T,
    lambda x: da.transpose(x, (1, 2, 0)),
    lambda x: x.sum(),
    lambda x: da.empty_like(x),
    lambda x: da.ones_like(x),
    lambda x: da.zeros_like(x),
    lambda x: da.full_like(x, 5),
    pytest.param(
        lambda x: x.mean(),
        marks=pytest.mark.skipif(
            cupy_version < Version("6.4.0"),
            reason="Requires CuPy 6.4.0+ "
            "(with https://github.com/cupy/cupy/pull/2418)",
        ),
    ),
    pytest.param(lambda x: x.moment(order=0)),
    lambda x: x.moment(order=2),
    pytest.param(
        lambda x: x.std(),
        marks=pytest.mark.skipif(
            cupy_version < Version("6.4.0"),
            reason="Requires CuPy 6.4.0+ "
            "(with https://github.com/cupy/cupy/pull/2418)",
        ),
    ),
    pytest.param(
        lambda x: x.var(),
        marks=pytest.mark.skipif(
            cupy_version < Version("6.4.0"),
            reason="Requires CuPy 6.4.0+ "
            "(with https://github.com/cupy/cupy/pull/2418)",
        ),
    ),
    pytest.param(
        lambda x: x.dot(np.arange(x.shape[-1])),
        marks=pytest.mark.xfail(reason="cupy.dot(numpy) fails"),
    ),
    pytest.param(
        lambda x: x.dot(np.eye(x.shape[-1])),
        marks=pytest.mark.xfail(reason="cupy.dot(numpy) fails"),
    ),
    pytest.param(
        lambda x: da.tensordot(x, np.ones(x.shape[:2]), axes=[(0, 1), (0, 1)]),
        marks=pytest.mark.xfail(reason="cupy.dot(numpy) fails"),
    ),
    lambda x: x.sum(axis=0),
    lambda x: x.max(axis=0),
    lambda x: x.sum(axis=(1, 2)),
    lambda x: x.astype(np.complex128),
    lambda x: x.map_blocks(lambda x: x * 2),
    pytest.param(lambda x: x.round(1)),
    lambda x: x.reshape((x.shape[0] * x.shape[1], x.shape[2])),
    # Rechunking here is required, see https://github.com/dask/dask/issues/2561
    lambda x: (x.rechunk(x.shape)).reshape((x.shape[1], x.shape[0], x.shape[2])),
    lambda x: x.reshape((x.shape[0], x.shape[1], x.shape[2] / 2, x.shape[2] / 2)),
    lambda x: abs(x),
    lambda x: x > 0.5,
    lambda x: x.rechunk((4, 4, 4)),
    lambda x: x.rechunk((2, 2, 1)),
    pytest.param(lambda x: da.einsum("ijk,ijk", x, x)),
    lambda x: np.isneginf(x),
    lambda x: np.isposinf(x),
    lambda x: np.isreal(x),
    lambda x: np.iscomplex(x),
    lambda x: np.real(x),
    lambda x: np.imag(x),
    lambda x: np.exp(x),
    lambda x: np.fix(x),
    lambda x: np.i0(x.reshape((24,))),
    lambda x: np.sinc(x),
    lambda x: np.nan_to_num(x),
    lambda x: np.max(x),
    lambda x: np.min(x),
    lambda x: np.prod(x),
    lambda x: np.any(x),
    lambda x: np.all(x),
    lambda x: np.nansum(x),
    lambda x: np.nanprod(x),
    lambda x: np.nanmin(x),
    lambda x: np.nanmax(x),
    pytest.param(
        lambda x: np.angle(x),
        marks=pytest.mark.skipif(
            not dask.utils.has_keyword(cupy.angle, "deg"),
            reason="Requires `deg` argument in `cupy.angle()` introduced in "
            "https://github.com/cupy/cupy/pull/6905",
        ),
    ),
    pytest.param(
        lambda x: np.angle(x, True),
        marks=pytest.mark.skipif(
            not dask.utils.has_keyword(cupy.angle, "deg"),
            reason="Requires `deg` argument in `cupy.angle()` introduced in "
            "https://github.com/cupy/cupy/pull/6905",
        ),
    ),
]


@pytest.mark.parametrize("func", functions)
def test_basic(func):
    c = cupy.random.default_rng().random((2, 3, 4))
    n = c.get()
    dc = da.from_array(c, chunks=(1, 2, 2), asarray=False)
    dn = da.from_array(n, chunks=(1, 2, 2))

    ddc = func(dc)
    ddn = func(dn)

    assert type(ddc._meta) is cupy.ndarray

    if next(iter(ddc.dask.keys()))[0].startswith("empty"):
        # We can't verify for data correctness when testing empty_like
        assert type(ddc._meta) is type(ddc.compute())
    else:
        assert_eq(ddc, ddc)  # Check that _meta and computed arrays match types
        assert_eq(ddc, ddn, check_type=False)


@pytest.mark.parametrize("dtype", ["f4", "f8"])
def test_sizeof(dtype):
    c = cupy.random.default_rng().random((2, 3, 4), dtype=dtype)

    assert sizeof(c) == c.nbytes


@pytest.mark.parametrize(
    "arr", [np.arange(5), cupy.arange(5), da.arange(5), da.from_array(cupy.arange(5))]
)
@pytest.mark.parametrize(
    "like", [np.arange(5), cupy.arange(5), da.arange(5), da.from_array(cupy.arange(5))]
)
def test_asanyarray(arr, like):
    if isinstance(like, np.ndarray) and isinstance(
        da.utils.meta_from_array(arr), cupy.ndarray
    ):
        with pytest.raises(TypeError):
            a = da.utils.asanyarray_safe(arr, like=like)
    else:
        a = da.utils.asanyarray_safe(arr, like=like)
        assert type(a) is type(like)


def test_vindex():
    x_np = np.arange(56).reshape((7, 8))
    x_cp = cupy.arange(56).reshape((7, 8))

    d_np = da.from_array(x_np, chunks=(3, 4))
    d_cp = da.from_array(x_cp, chunks=(3, 4))

    res_np = da.core._vindex(d_np, [0, 1, 6, 0], [0, 1, 0, 7])
    res_cp = da.core._vindex(d_cp, [0, 1, 6, 0], [0, 1, 0, 7])

    assert type(res_cp._meta) == cupy.ndarray
    assert_eq(
        res_cp, res_cp, check_type=False
    )  # Check that _meta and computed arrays match types

    assert_eq(res_np, res_cp, check_type=False)


def test_view():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(cupy.array(x), chunks=(2, 3))

    result = d.view()
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, x.view(), check_type=False)

    result = d.view("i4")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, x.view("i4"), check_type=False)

    result = d.view("i2")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, x.view("i2"), check_type=False)
    assert all(isinstance(s, int) for s in d.shape)

    x = np.arange(8, dtype="i1")
    d = da.from_array(cupy.array(x), chunks=(4,))
    result = d.view("i4")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(x.view("i4"), d.view("i4"), check_type=False)

    with pytest.raises(ValueError):
        x = np.arange(8, dtype="i1")
        d = da.from_array(cupy.array(x), chunks=(3,))
        d.view("i4")

    with pytest.raises(ValueError):
        d.view("i4", order="asdf")


def test_view_fortran():
    x = np.asfortranarray(np.arange(64).reshape((8, 8)))
    d = da.from_array(cupy.asfortranarray(cupy.array(x)), chunks=(2, 3))

    result = d.view("i4", order="F")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, x.T.view("i4").T, check_type=False)

    result = d.view("i2", order="F")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, x.T.view("i2").T, check_type=False)


def test_getter():
    result = da.core.getter(cupy.arange(5), (None, slice(None, None)))

    assert type(result) == cupy.ndarray
    assert_eq(result, np.arange(5)[None, :], check_type=False)


def test_store_kwargs():
    d = da.from_array(cupy.ones((10, 10)), chunks=(2, 2))
    a = d + 1

    called = [False]

    def get_func(*args, **kwargs):
        assert kwargs.pop("foo") == "test kwarg"
        r = dask.get(*args, **kwargs)
        called[0] = True
        return r

    called[0] = False
    at = cupy.zeros(shape=(10, 10))
    da.core.store([a], [at], scheduler=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = cupy.zeros(shape=(10, 10))
    a.store(at, scheduler=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = cupy.zeros(shape=(10, 10))
    da.core.store([a], [at], scheduler=get_func, return_stored=True, foo="test kwarg")
    assert called[0]


def test_setitem_1d():
    x = cupy.arange(10)
    dx = da.from_array(x.copy(), chunks=(5,))

    x[x > 6] = -1
    x[x % 2 == 0] = -2

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2

    assert_eq(x, dx)


def test_setitem_2d():
    x = cupy.arange(24).reshape((4, 6))
    dx = da.from_array(x.copy(), chunks=(2, 2))

    x[x > 6] = -1
    x[x % 2 == 0] = -2

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2

    assert_eq(x, dx)


def test_setitem_extended_API_0d():
    # 0-d array
    x = cupy.array(9)
    dx = da.from_array(x.copy())

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
        pytest.param(
            slice(8, None, 2),
            [-30],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
        [slice(1, None, -2), -4],
        pytest.param(
            slice(1, None, -2),
            [-40],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
        [slice(3, None, 2), -5],
        [slice(-3, None, -2), -6],
        [slice(1, None, -2), -4],
        [slice(3, None, 2), -5],
        pytest.param(
            slice(3, None, 2),
            [10, 11, 12, 13],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
        pytest.param(
            slice(-4, None, -2),
            [14, 15, 16, 17],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
    ],
)
def test_setitem_extended_API_1d(index, value):
    # 1-d array
    x = cupy.arange(10)
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
        pytest.param(
            (slice(None), 2),
            range(6),
            marks=pytest.mark.skip(
                reason="Assigning `range` to CuPy array is not supported"
            ),
        ),
        pytest.param(
            3,
            range(10),
            marks=pytest.mark.skip(
                reason="Assigning `range` to CuPy array is not supported"
            ),
        ),
        [(slice(None), [3, 5, 6]), [-30, -31, -32]],
        [([-1, 0, 1], 2), [-30, -31, -32]],
        pytest.param(
            (slice(None, 2), slice(None, 3)),
            [-50, -51, -52],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
        [(slice(None), [6, 1, 3]), [-60, -61, -62]],
        pytest.param(
            (slice(1, 3), slice(1, 4)),
            [[-70, -71, -72]],
            marks=pytest.mark.skip(reason="Unsupported assigning `list` to CuPy array"),
        ),
        pytest.param(
            (slice(None), [9, 8, 8]),
            [[-80, -81, 91]],
            marks=pytest.mark.flaky(reruns=10),
        ),
        [([True, False, False, False, True, False], 2), -1],
        [(3, [True, True, False, True, True, False, True, False, True, True]), -1],
        [(np.array([False, False, True, True, False, False]), slice(5, 7)), -1],
        [(cupy.array([False, False, True, True, False, False]), slice(5, 7)), -1],
        pytest.param(
            (
                4,
                da.from_array(
                    [False, False, True, True, False, False, True, False, False, True]
                ),
            ),
            -1,
            marks=pytest.mark.skip(
                reason="Unsupported assigning Dask Array to CuPy array"
            ),
        ),
        [slice(5, None, 2), -99],
        pytest.param(
            slice(5, None, 2),
            range(1, 11),
            marks=pytest.mark.skip(
                reason="Assigning `range` to CuPy array is not supported"
            ),
        ),
        [slice(1, None, -2), -98],
        pytest.param(
            slice(1, None, -2),
            range(11, 21),
            marks=pytest.mark.skip(
                reason="Assigning `range` to CuPy array is not supported"
            ),
        ),
    ],
)
def test_setitem_extended_API_2d(index, value):
    # 2-d array
    x = cupy.arange(60).reshape((6, 10))
    dx = da.from_array(x, chunks=(2, 3))
    dx[index] = value
    x[index] = value
    assert_eq(x, dx.compute())


def test_setitem_extended_API_2d_rhs_func_of_lhs():
    # Cases:
    # * RHS and/or indices are a function of the LHS
    # * Indices have unknown chunk sizes
    # * RHS has extra leading size 1 dimensions compared to LHS
    x = cupy.arange(60).reshape((6, 10))
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

    # Both tests below fail in CuPy due to leading singular dimensions
    if False:
        # RHS has extra leading size 1 dimensions compared to LHS
        dx = da.from_array(x.copy(), chunks=(2, 3))
        v = x.reshape((1, 1) + x.shape)
        x[...] = v
        dx[...] = v
        assert_eq(x, dx.compute())

        index = da.where(da.arange(3, chunks=(1,)) < 2)[0]
        v = -cupy.arange(12).reshape(1, 1, 6, 2)
        x[:, [0, 1]] = v
        dx[:, index] = v
        assert_eq(x, dx.compute())


def test_setitem_on_read_only_blocks():
    # Outputs of broadcast_trick-style functions contain read-only
    # arrays
    dx = da.empty_like(cupy.array(()), shape=(4, 6), dtype=float, chunks=(2, 2))
    dx[0] = 99

    assert_eq(dx[0, 0], 99.0)

    dx[0:2] = 88

    assert_eq(dx[0, 0], 88.0)


def test_setitem_errs():
    x = da.ones_like(cupy.array(()), shape=(4, 4), chunks=(2, 2))

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

    # Too many/not enough booleans
    with pytest.raises(IndexError):
        x[[True, True, False]] = 5

    with pytest.raises(IndexError):
        x[[False, True, True, True, False]] = 5

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
        x[cupy.array(True)] = 5

    with pytest.raises(NotImplementedError):
        x[0, da.from_array(True)] = 5

    # Scalar arrays
    y = da.from_array(cupy.array(1))
    with pytest.raises(IndexError):
        y[:] = 2

    # RHS has non-brodacastable extra leading dimensions
    x = cupy.arange(12).reshape((3, 4))
    dx = da.from_array(x, chunks=(2, 2))
    with pytest.raises(ValueError):
        dx[...] = cupy.arange(24).reshape((2, 1, 3, 4))

    # RHS has extra leading size 1 dimensions compared to LHS
    x = cupy.arange(12).reshape((3, 4))
    dx = da.from_array(x, chunks=(2, 3))


@pytest.mark.parametrize("xp", [np, da])
@pytest.mark.parametrize("orig_arr", [np.array, da.array])
@pytest.mark.parametrize("array_func", ["array", "asarray", "asanyarray"])
def test_array_like(xp, orig_arr, array_func):
    cp_func = getattr(cupy, array_func)
    xp_func = getattr(xp, array_func)

    cp_a = cp_func([1, 2, 3])
    xp_a = xp_func(orig_arr([1, 2, 3]), like=da.from_array(cupy.array(())))
    assert isinstance(xp_a, da.Array)
    assert isinstance(xp_a._meta, cupy.ndarray)
    assert_eq(xp_a, cp_a)
