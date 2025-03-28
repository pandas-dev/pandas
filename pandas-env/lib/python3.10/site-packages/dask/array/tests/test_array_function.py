from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

import dask.array as da
from dask.array.tests.test_dispatch import EncapsulateNDArray, WrappedArray
from dask.array.utils import assert_eq

if da._array_expr_enabled():
    pytest.skip("parametrize using unsupported functions", allow_module_level=True)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: np.append(x, x),
        lambda x: np.concatenate([x, x, x]),
        lambda x: np.cov(x, x),
        lambda x: np.dot(x, x),
        lambda x: np.dstack((x, x)),
        lambda x: np.flip(x, axis=0),
        lambda x: np.hstack((x, x)),
        lambda x: np.matmul(x, x),
        lambda x: np.mean(x),
        lambda x: np.stack([x, x]),
        lambda x: np.block([x, x]),
        lambda x: np.sum(x),
        lambda x: np.var(x),
        lambda x: np.vstack((x, x)),
        lambda x: np.linalg.norm(x),
        lambda x: np.min(x),
        lambda x: np.amin(x),
        lambda x: np.round(x),
        lambda x: np.insert(x, 0, 3, axis=0),
        lambda x: np.delete(x, 0, axis=0),
        lambda x: np.select(
            [x < 0.3, x < 0.6, x > 0.7], [x * 2, x, x / 2], default=0.65
        ),
    ],
)
def test_array_function_dask(func):
    x = np.random.default_rng().random((100, 100))
    y = da.from_array(x, chunks=(50, 50))
    res_x = func(x)
    res_y = func(y)

    assert isinstance(res_y, da.Array)
    assert_eq(res_y, res_x)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: np.dstack(x),
        lambda x: np.hstack(x),
        lambda x: np.vstack(x),
    ],
)
def test_stack_functions_require_sequence_of_arrays(func):
    x = np.random.default_rng().random((100, 100))
    y = da.from_array(x, chunks=(50, 50))

    with pytest.raises(
        NotImplementedError, match="expects a sequence of arrays as the first argument"
    ):
        func(y)


@pytest.mark.parametrize("func", [np.fft.fft, np.fft.fft2])
def test_array_function_fft(func):
    x = np.random.default_rng().random((100, 100))
    y = da.from_array(x, chunks=(100, 100))
    res_x = func(x)
    res_y = func(y)

    if func.__module__ != "mkl_fft._numpy_fft":
        assert isinstance(res_y, da.Array)
    assert_eq(res_y, res_x)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: np.min_scalar_type(x),
        lambda x: np.linalg.det(x),
        lambda x: np.linalg.eigvals(x),
    ],
)
def test_array_notimpl_function_dask(func):
    x = np.random.default_rng().random((100, 100))
    y = da.from_array(x, chunks=(50, 50))

    with pytest.warns(
        FutureWarning, match="The `.*` function is not implemented by Dask"
    ):
        func(y)


@pytest.mark.parametrize(
    "func", [lambda x: np.real(x), lambda x: np.imag(x), lambda x: np.transpose(x)]
)
def test_array_function_sparse(func):
    sparse = pytest.importorskip("sparse")
    if Version(sparse.__version__) == Version("0.15.2"):
        pytest.skip(reason="https://github.com/pydata/sparse/issues/682")
    x = da.random.default_rng().random((500, 500), chunks=(100, 100))
    x[x < 0.9] = 0

    y = x.map_blocks(sparse.COO)

    assert_eq(func(x), func(y))


def test_array_function_sparse_tensordot():
    sparse = pytest.importorskip("sparse")
    rng = np.random.default_rng()
    x = rng.random((2, 3, 4))
    x[x < 0.9] = 0
    y = rng.random((4, 3, 2))
    y[y < 0.9] = 0

    xx = sparse.COO(x)
    yy = sparse.COO(y)

    assert_eq(
        np.tensordot(x, y, axes=(2, 0)), np.tensordot(xx, yy, axes=(2, 0)).todense()
    )


@pytest.mark.parametrize("chunks", [(100, 100), (500, 100)])
def test_array_function_cupy_svd(chunks):
    cupy = pytest.importorskip("cupy")
    x = cupy.random.default_rng().random((500, 100))

    y = da.from_array(x, chunks=chunks, asarray=False)

    u_base, s_base, v_base = da.linalg.svd(y)
    u, s, v = np.linalg.svd(y)

    assert_eq(u, u_base)
    assert_eq(s, s_base)
    assert_eq(v, v_base)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: np.concatenate([x, x, x]),
        lambda x: np.cov(x, x),
        lambda x: np.dot(x, x),
        lambda x: np.dstack((x, x)),
        lambda x: np.flip(x, axis=0),
        lambda x: np.hstack((x, x)),
        lambda x: np.matmul(x, x),
        lambda x: np.mean(x),
        lambda x: np.stack([x, x]),
        lambda x: np.sum(x),
        lambda x: np.var(x),
        lambda x: np.vstack((x, x)),
        lambda x: np.linalg.norm(x),
    ],
)
def test_unregistered_func(func):
    # Wrap a procol-based encapsulated ndarray
    x = EncapsulateNDArray(np.random.default_rng().random((100, 100)))

    # See if Dask holds the array fine
    y = da.from_array(x, chunks=(50, 50))

    # Check if it's an equivalent array
    assert_eq(x, y, check_meta=False, check_type=False)

    # Perform two NumPy functions, one on the
    # Encapsulated array
    xx = func(x)

    # And one on the Dask array holding these
    # encapsulated arrays
    yy = func(y)

    # Check that they are equivalent arrays.
    assert_eq(xx, yy, check_meta=False, check_type=False)


def test_non_existent_func():
    # Regression test for __array_function__ becoming default in numpy 1.17
    # dask has no sort function, so ensure that this still calls np.sort
    x = da.from_array(np.array([1, 2, 4, 3]), chunks=(2,))
    with pytest.warns(
        FutureWarning, match="The `numpy.sort` function is not implemented by Dask"
    ):
        assert list(np.sort(x)) == [1, 2, 3, 4]


@pytest.mark.parametrize(
    "func",
    [
        np.equal,
        np.matmul,
        np.dot,
        lambda x, y: np.stack([x, y]),
    ],
)
@pytest.mark.parametrize(
    "arr_upcast, arr_downcast",
    [
        (
            WrappedArray(np.random.default_rng().random((10, 10))),
            da.random.default_rng().random((10, 10), chunks=(5, 5)),
        ),
        (
            da.random.default_rng().random((10, 10), chunks=(5, 5)),
            EncapsulateNDArray(np.random.default_rng().random((10, 10))),
        ),
        (
            WrappedArray(np.random.default_rng().random((10, 10))),
            EncapsulateNDArray(np.random.default_rng().random((10, 10))),
        ),
    ],
)
def test_binary_function_type_precedence(func, arr_upcast, arr_downcast):
    """Test proper dispatch on binary NumPy functions"""
    assert (
        type(func(arr_upcast, arr_downcast))
        == type(func(arr_downcast, arr_upcast))
        == type(arr_upcast)
    )


@pytest.mark.parametrize("like_ns", (np, da))
@pytest.mark.parametrize("func", [da.array, da.asarray, da.asanyarray])
def test_like(func, like_ns):
    like = like_ns.array(0, dtype=np.int16)

    assert_eq(func(1, like=like), da.array(1))
    assert_eq(func(1.1, like=like), da.array(1.1))
    assert_eq(func(np.int8(1), like=like, dtype=np.int32), da.array(1, dtype=np.int32))
    assert_eq(func(1, like=like, dtype=np.int32), da.array(1, dtype=np.int32))
    assert_eq(func(1.1, like=like, dtype=np.int32), da.array(1, dtype=np.int32))

    assert_eq(func(da.array(1, dtype=np.int8), like=like), da.array(1, dtype=np.int8))
    assert_eq(func(da.array(1.1), like=like), da.array(1.1))
    assert_eq(
        func(da.array(1, dtype=np.int8), like=like, dtype=np.int32),
        da.array(1, dtype=np.int32),
    )


@pytest.mark.parametrize("like_ns", (np, da))
def test_tri_like(like_ns):
    like = like_ns.array(0)
    assert_eq(da.tri(1, like=like), da.tri(1))
    assert_eq(da.tri(1, like=like, dtype=np.float32), da.tri(1, dtype=np.float32))


@pytest.mark.parametrize("func", [np.array, np.asarray, np.asanyarray])
def test_like_with_numpy_func(func):
    assert_eq(func(1, like=da.array(1)), func(1))


@pytest.mark.parametrize("func", [np.array, np.asarray, np.asanyarray])
def test_like_with_numpy_func_and_dtype(func):
    assert_eq(func(1, dtype=float, like=da.array(1)), func(1, dtype=float))
