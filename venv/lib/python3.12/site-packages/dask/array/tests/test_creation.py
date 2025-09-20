from __future__ import annotations

import pickle

import pytest

from dask._task_spec import Alias
from dask.base import collections_to_expr

pytest.importorskip("numpy")

import numpy as np
import pytest
from tlz import concat

import dask
import dask.array as da
from dask.array.core import normalize_chunks
from dask.array.numpy_compat import NUMPY_GE_210, AxisError
from dask.array.utils import assert_eq, same_keys


@pytest.mark.parametrize(
    "backend",
    ["numpy", pytest.param("cupy", marks=pytest.mark.gpu)],
)
@pytest.mark.parametrize(
    "funcname",
    [
        "empty_like",
        "empty",
        "ones_like",
        "ones",
        "zeros_like",
        "zeros",
        "full_like",
        "full",
    ],
)
@pytest.mark.parametrize("cast_shape", [tuple, list, np.asarray])
@pytest.mark.parametrize("cast_chunks", [tuple, list, np.asarray])
@pytest.mark.parametrize("shape, chunks", [((10, 10), (4, 4))])
@pytest.mark.parametrize("name", [None, "my-name"])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("dtype", ["i4"])
def test_arr_like(
    funcname, shape, cast_shape, dtype, cast_chunks, chunks, name, order, backend
):
    backend_lib = pytest.importorskip(backend)
    with dask.config.set({"array.backend": backend}):
        np_func = getattr(backend_lib, funcname)
        da_func = getattr(da, funcname)
        shape = cast_shape(shape)
        chunks = cast_chunks(chunks)

        if "full" in funcname:
            old_np_func = np_func
            old_da_func = da_func

            np_func = lambda *a, **k: old_np_func(*a, fill_value=5, **k)
            da_func = lambda *a, **k: old_da_func(*a, fill_value=5, **k)

        dtype = np.dtype(dtype)

        if "like" in funcname:
            a = backend_lib.random.randint(0, 10, shape).astype(dtype)

            np_r = np_func(a, order=order)
            da_r = da_func(a, order=order, chunks=chunks, name=name)
        else:
            np_r = np_func(shape, order=order, dtype=dtype)
            da_r = da_func(shape, order=order, dtype=dtype, chunks=chunks, name=name)

        assert np_r.shape == da_r.shape
        assert np_r.dtype == da_r.dtype

        # Make sure we are using the desired backend
        assert isinstance(da_r._meta, backend_lib.ndarray)
        assert isinstance(da_r.compute(), backend_lib.ndarray)

        if "empty" not in funcname:
            assert_eq(np_r, da_r)

        if name is None:
            assert funcname.split("_")[0] in da_r.name
        else:
            assert da_r.name == name

        if "order" == "F":
            assert np.isfortran(da_r.compute())
        else:
            assert not np.isfortran(da_r.compute())


@pytest.mark.parametrize(
    "funcname, kwargs",
    [
        ("empty_like", {}),
        ("ones_like", {}),
        ("zeros_like", {}),
        ("full_like", {"fill_value": 5}),
    ],
)
@pytest.mark.parametrize(
    "shape, chunks, out_shape",
    [
        ((10, 10), (4, 4), None),
        ((10, 10), (4, 4), (20, 3)),
        ((10, 10), (4), (20)),
        ((10, 10, 10), (4, 2), (5, 5)),
        ((2, 3, 5, 7), None, (3, 5, 7)),
        ((2, 3, 5, 7), (2, 5, 3), (3, 5, 7)),
        ((2, 3, 5, 7), (2, 5, 3, "auto", 3), (11,) + (2, 3, 5, 7)),
        ((2, 3, 5, 7), "auto", (3, 5, 7)),
    ],
)
@pytest.mark.parametrize("dtype", ["i4"])
def test_arr_like_shape(funcname, kwargs, shape, dtype, chunks, out_shape):
    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)
    a = np.random.randint(0, 10, shape).astype(dtype)
    np_r = np_func(a, shape=out_shape, **kwargs)
    da_r = da_func(a, chunks=chunks, shape=out_shape, **kwargs)

    assert np_r.shape == da_r.shape
    assert np_r.dtype == da_r.dtype

    if "empty" not in funcname:
        assert_eq(np_r, da_r)


@pytest.mark.parametrize("endpoint", [True, False])
def test_linspace(endpoint):
    darr = da.linspace(6, 49, endpoint=endpoint, chunks=5)
    nparr = np.linspace(6, 49, endpoint=endpoint)
    assert_eq(darr, nparr)

    darr = da.linspace(1.4, 4.9, endpoint=endpoint, chunks=5, num=13)
    nparr = np.linspace(1.4, 4.9, endpoint=endpoint, num=13)
    assert_eq(darr, nparr)

    darr = da.linspace(6, 49, endpoint=endpoint, chunks=5, dtype=float)
    nparr = np.linspace(6, 49, endpoint=endpoint, dtype=float)
    assert_eq(darr, nparr)

    darr, dstep = da.linspace(6, 49, endpoint=endpoint, chunks=5, retstep=True)
    nparr, npstep = np.linspace(6, 49, endpoint=endpoint, retstep=True)
    assert np.allclose(dstep, npstep)
    assert_eq(darr, nparr)

    darr = da.linspace(1.4, 4.9, endpoint=endpoint, chunks=5, num=13, dtype=int)
    nparr = np.linspace(1.4, 4.9, num=13, endpoint=endpoint, dtype=int)
    assert_eq(darr, nparr)
    assert sorted(
        da.linspace(1.4, 4.9, endpoint=endpoint, chunks=5, num=13).dask
    ) == sorted(da.linspace(1.4, 4.9, endpoint=endpoint, chunks=5, num=13).dask)
    assert sorted(
        da.linspace(6, 49, endpoint=endpoint, chunks=5, dtype=float).dask
    ) == sorted(da.linspace(6, 49, endpoint=endpoint, chunks=5, dtype=float).dask)

    x = da.array([0.2, 6.4, 3.0, 1.6])
    nparr = np.linspace(0, 2, 8, endpoint=endpoint)
    darr = da.linspace(da.argmin(x), da.argmax(x) + 1, 8, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(0, 0, 0, endpoint=endpoint)
    darr = da.linspace(0, 0, 0, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(1, 1, 0, endpoint=endpoint)
    darr = da.linspace(1, 1, 0, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(1, 5, 0, endpoint=endpoint)
    darr = da.linspace(1, 5, 0, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(0, 0, 1, endpoint=endpoint)
    darr = da.linspace(0, 0, 1, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(1, 1, 1, endpoint=endpoint)
    darr = da.linspace(1, 1, 1, endpoint=endpoint)
    assert_eq(darr, nparr)

    nparr = np.linspace(1, 5, 1, endpoint=endpoint)
    darr = da.linspace(1, 5, 1, endpoint=endpoint)
    assert_eq(darr, nparr)


def test_arange():
    darr = da.arange(77, chunks=13)
    nparr = np.arange(77)
    assert_eq(darr, nparr)

    darr = da.arange(2, 13, chunks=5)
    nparr = np.arange(2, 13)
    assert_eq(darr, nparr)

    darr = da.arange(4, 21, 9, chunks=13)
    nparr = np.arange(4, 21, 9)
    assert_eq(darr, nparr)

    # negative steps
    darr = da.arange(53, 5, -3, chunks=5)
    nparr = np.arange(53, 5, -3)
    assert_eq(darr, nparr)

    darr = da.arange(77, chunks=13, dtype=float)
    nparr = np.arange(77, dtype=float)
    assert_eq(darr, nparr)

    darr = da.arange(2, 13, chunks=5, dtype=int)
    nparr = np.arange(2, 13, dtype=int)
    assert_eq(darr, nparr)
    assert sorted(da.arange(2, 13, chunks=5).dask) == sorted(
        da.arange(2, 13, chunks=5).dask
    )
    assert sorted(da.arange(77, chunks=13, dtype=float).dask) == sorted(
        da.arange(77, chunks=13, dtype=float).dask
    )

    # 0 size output
    darr = da.arange(0, 1, -0.5, chunks=20)
    nparr = np.arange(0, 1, -0.5)
    assert_eq(darr, nparr)

    darr = da.arange(0, -1, 0.5, chunks=20)
    nparr = np.arange(0, -1, 0.5)
    assert_eq(darr, nparr)

    # stop and/or step as kwargs
    darr = da.arange(stop=10)
    nparr = np.arange(stop=10)
    assert_eq(darr, nparr)

    darr = da.arange(10, step=2)
    nparr = np.arange(10, step=2)
    assert_eq(darr, nparr)

    darr = da.arange(stop=10, step=2)
    nparr = np.arange(stop=10, step=2)
    assert_eq(darr, nparr)

    darr = da.arange(3, stop=10, step=2)
    nparr = np.arange(3, stop=10, step=2)
    assert_eq(darr, nparr)

    with pytest.raises(TypeError, match="requires stop"):
        da.arange()

    # Unexpected or missing kwargs
    with pytest.raises(TypeError, match="whatsthis"):
        da.arange(10, chunks=-1, whatsthis=1)

    assert da.arange(10).chunks == ((10,),)


arange_dtypes = [
    np.uint8,
    np.uint64,
    np.int8,
    np.int64,
    np.float32,
    np.float64,
]


# FIXME hypothesis would be much better suited for this
@pytest.mark.parametrize("start_type", arange_dtypes + [int, float])
@pytest.mark.parametrize("stop_type", arange_dtypes + [int, float])
@pytest.mark.parametrize("step_type", arange_dtypes + [int, float])
def test_arange_dtype_infer(start_type, stop_type, step_type):
    start = start_type(3)
    stop = stop_type(13)
    step = step_type(2)
    a_np = np.arange(start, stop, step)
    a_da = da.arange(start, stop, step)
    assert_eq(a_np, a_da)


@pytest.mark.parametrize("dtype", arange_dtypes)
def test_arange_dtype_force(dtype):
    assert da.arange(10, dtype=dtype).dtype == dtype
    assert da.arange(np.float32(10), dtype=dtype).dtype == dtype
    assert da.arange(np.int64(10), dtype=dtype).dtype == dtype


@pytest.mark.skipif(np.array(0).dtype == np.int32, reason="64-bit platforms only")
@pytest.mark.parametrize(
    "start,stop,step",
    [
        (2**63 - 10_000, 2**63 - 1, 1),
        (2**63 - 1, 2**63 - 10_000, -1),
        (0, 2**63 - 1, 2**63 - 10_000),
        (0.0, 2**63 - 1, 2**63 - 10_000),
        (0.0, -9_131_138_316_486_228_481, -92_233_720_368_547_759),
        (-72_057_594_037_927_945, -72_057_594_037_927_938, 1.0),
        (-72_057_594_037_927_945, -72_057_594_037_927_938, 1.5),
    ],
)
@pytest.mark.parametrize("chunks", ["auto", 1])
def test_arange_very_large_args(start, stop, step, chunks):
    """Test args that are very close to 2**63
    https://github.com/dask/dask/issues/11706
    """
    a_np = np.arange(start, stop, step)
    a_da = da.arange(start, stop, step, chunks=chunks)
    assert_eq(a_np, a_da)


@pytest.mark.xfail(
    reason="Casting floats to ints is not supported since edge "
    "behavior is not specified or guaranteed by NumPy."
)
def test_arange_cast_float_int_step():
    darr = da.arange(3.3, -9.1, -0.25, chunks=3, dtype="i8")
    nparr = np.arange(3.3, -9.1, -0.25, dtype="i8")
    assert_eq(darr, nparr)


def test_arange_float_step():
    darr = da.arange(2.0, 13.0, 0.3, chunks=4)
    nparr = np.arange(2.0, 13.0, 0.3)
    assert_eq(darr, nparr)

    darr = da.arange(7.7, 1.5, -0.8, chunks=3)
    nparr = np.arange(7.7, 1.5, -0.8)
    assert_eq(darr, nparr)

    darr = da.arange(0, 1, 0.01, chunks=20)
    nparr = np.arange(0, 1, 0.01)
    assert_eq(darr, nparr)

    darr = da.arange(0, 1, 0.03, chunks=20)
    nparr = np.arange(0, 1, 0.03)
    assert_eq(darr, nparr)


def test_indices_wrong_chunks():
    with pytest.raises(ValueError):
        da.indices((1,), chunks=tuple())


def test_indices_dimensions_chunks():
    chunks = ((1, 4, 2, 3), (5, 5))
    darr = da.indices((10, 10), chunks=chunks)
    assert darr.chunks == ((1, 1),) + chunks

    with dask.config.set({"array.chunk-size": "50 MiB"}):
        shape = (10000, 10000)
        expected = normalize_chunks("auto", shape=shape, dtype=int)
        result = da.indices(shape, chunks="auto")
        # indices prepends a dimension
        actual = result.chunks[1:]
        assert expected == actual


def test_empty_indices():
    darr = da.indices(tuple(), chunks=tuple())
    nparr = np.indices(tuple())
    assert darr.shape == nparr.shape
    assert darr.dtype == nparr.dtype
    assert_eq(darr, nparr)

    darr = da.indices(tuple(), float, chunks=tuple())
    nparr = np.indices(tuple(), float)
    assert darr.shape == nparr.shape
    assert darr.dtype == nparr.dtype
    assert_eq(darr, nparr)

    darr = da.indices((0,), float, chunks=(1,))
    nparr = np.indices((0,), float)
    assert darr.shape == nparr.shape
    assert darr.dtype == nparr.dtype
    assert_eq(darr, nparr)

    darr = da.indices((0, 1, 2), float, chunks=(1, 1, 2))
    nparr = np.indices((0, 1, 2), float)
    assert darr.shape == nparr.shape
    assert darr.dtype == nparr.dtype
    assert_eq(darr, nparr)


def test_indices():
    darr = da.indices((1,), chunks=(1,))
    nparr = np.indices((1,))
    assert_eq(darr, nparr)

    darr = da.indices((1,), float, chunks=(1,))
    nparr = np.indices((1,), float)
    assert_eq(darr, nparr)

    darr = da.indices((2, 1), chunks=(2, 1))
    nparr = np.indices((2, 1))
    assert_eq(darr, nparr)

    darr = da.indices((2, 3), chunks=(1, 2))
    nparr = np.indices((2, 3))
    assert_eq(darr, nparr)


@pytest.mark.parametrize(
    "shapes, chunks",
    [
        ([()], [()]),
        ([(0,)], [(0,)]),
        ([(2,), (3,)], [(1,), (2,)]),
        ([(2,), (3,), (4,)], [(1,), (2,), (3,)]),
        ([(2,), (3,), (4,), (5,)], [(1,), (2,), (3,), (4,)]),
        ([(2, 3), (4,)], [(1, 2), (3,)]),
    ],
)
@pytest.mark.parametrize("indexing", ["ij", "xy"])
@pytest.mark.parametrize("sparse", [False, True])
def test_meshgrid(shapes, chunks, indexing, sparse):
    xi_a = []
    xi_d = []
    xi_dc = []
    for each_shape, each_chunk in zip(shapes, chunks):
        xi_a.append(np.random.random(each_shape))
        xi_d_e = da.from_array(xi_a[-1], chunks=each_chunk)
        xi_d.append(xi_d_e)
        xi_d_ef = xi_d_e.flatten()
        xi_dc.append(xi_d_ef.chunks[0])
    do = list(range(len(xi_dc)))
    if indexing == "xy" and len(xi_dc) > 1:
        do[0], do[1] = do[1], do[0]
        xi_dc[0], xi_dc[1] = xi_dc[1], xi_dc[0]
    xi_dc = tuple(xi_dc)

    r_a = np.meshgrid(*xi_a, indexing=indexing, sparse=sparse)
    r_d = da.meshgrid(*xi_d, indexing=indexing, sparse=sparse)

    assert type(r_d) is type(r_a)
    assert len(r_a) == len(r_d)

    for e_r_a, e_r_d, i in zip(r_a, r_d, do):
        assert_eq(e_r_a, e_r_d)
        if sparse:
            assert e_r_d.chunks[i] == xi_dc[i]
        else:
            assert e_r_d.chunks == xi_dc


def test_meshgrid_inputcoercion():
    a = [1, 2, 3]
    b = np.array([4, 5, 6, 7])
    x, y = np.meshgrid(a, b, indexing="ij")
    z = x * y

    x_d, y_d = da.meshgrid(a, b, indexing="ij")
    z_d = x_d * y_d

    assert z_d.shape == (len(a), len(b))
    assert_eq(z, z_d)


@pytest.mark.parametrize(
    "N, M, k, dtype, chunks",
    [
        (3, None, 0, float, "auto"),
        (4, None, 0, float, "auto"),
        (3, 4, 0, bool, "auto"),
        (3, None, 1, int, "auto"),
        (3, None, -1, int, "auto"),
        (3, None, 2, int, 1),
        (6, 8, -2, int, (3, 4)),
        (6, 8, 0, int, (3, "auto")),
    ],
)
def test_tri(N, M, k, dtype, chunks):
    assert_eq(da.tri(N, M, k, dtype, chunks), np.tri(N, M, k, dtype))


def test_eye():
    assert_eq(da.eye(9, chunks=3), np.eye(9))
    assert_eq(da.eye(9), np.eye(9))
    assert_eq(da.eye(10, chunks=3), np.eye(10))
    assert_eq(da.eye(9, chunks=3, M=11), np.eye(9, M=11))
    assert_eq(da.eye(11, chunks=3, M=9), np.eye(11, M=9))
    assert_eq(da.eye(7, chunks=3, M=11), np.eye(7, M=11))
    assert_eq(da.eye(11, chunks=3, M=7), np.eye(11, M=7))
    assert_eq(da.eye(9, chunks=3, k=2), np.eye(9, k=2))
    assert_eq(da.eye(9, chunks=3, k=-2), np.eye(9, k=-2))
    assert_eq(da.eye(7, chunks=3, M=11, k=5), np.eye(7, M=11, k=5))
    assert_eq(da.eye(11, chunks=3, M=7, k=-6), np.eye(11, M=7, k=-6))
    assert_eq(da.eye(6, chunks=3, M=9, k=7), np.eye(6, M=9, k=7))
    assert_eq(da.eye(12, chunks=3, M=6, k=-3), np.eye(12, M=6, k=-3))

    assert_eq(da.eye(9, chunks=3, dtype=int), np.eye(9, dtype=int))
    assert_eq(da.eye(10, chunks=3, dtype=int), np.eye(10, dtype=int))
    assert_eq(da.eye(10, chunks=-1, dtype=int), np.eye(10, dtype=int))
    assert_eq(da.eye(9, chunks=3, dtype=None), np.eye(9, dtype=None))

    with dask.config.set({"array.chunk-size": "50 MiB"}):
        x = da.eye(10000, "auto")
        assert 4 < x.npartitions < 32


@pytest.mark.parametrize("k", [0, 3, -3, 8])
def test_diag_bad_input(k):
    # when input numpy array is neither 1d nor 2d:
    v = np.arange(2 * 3 * 4).reshape((2, 3, 4))
    with pytest.raises(ValueError, match="Array must be 1d or 2d only"):
        da.diag(v, k)

    # when input dask array is neither 1d nor 2d:
    v = da.arange(2 * 3 * 4).reshape((2, 3, 4))
    with pytest.raises(ValueError, match="Array must be 1d or 2d only"):
        da.diag(v, k)

    # when input is not an array:
    v = 1
    with pytest.raises(TypeError, match="v must be a dask array or numpy array"):
        da.diag(v, k)


@pytest.mark.parametrize("k", [0, 3, -3, 8])
def test_diag_2d_array_creation(k):
    # when input 1d-array is a numpy array:
    v = np.arange(11)
    assert_eq(da.diag(v, k), np.diag(v, k))

    # when input 1d-array is a dask array:
    v = da.arange(11, chunks=3)
    darr = da.diag(v, k)
    nparr = np.diag(v, k)
    assert_eq(darr, nparr)
    assert sorted(da.diag(v, k).dask) == sorted(da.diag(v, k).dask)

    v = v + v + 3
    darr = da.diag(v, k)
    nparr = np.diag(v, k)
    assert_eq(darr, nparr)

    v = da.arange(11, chunks=11)
    darr = da.diag(v, k)
    nparr = np.diag(v, k)
    assert_eq(darr, nparr)
    assert sorted(da.diag(v, k).dask) == sorted(da.diag(v, k).dask)


@pytest.mark.parametrize("k", [0, 3, -3, 8])
def test_diag_extraction(k):
    # when input 2d-array is a square numpy array:
    x = np.arange(64).reshape((8, 8))
    assert_eq(da.diag(x, k), np.diag(x, k))
    # when input 2d-array is a square dask array:
    d = da.from_array(x, chunks=(4, 4))
    assert_eq(da.diag(d, k), np.diag(x, k))
    # heterogeneous chunks:
    d = da.from_array(x, chunks=((3, 2, 3), (4, 1, 2, 1)))
    assert_eq(da.diag(d, k), np.diag(x, k))

    # when input 2d-array is a rectangular numpy array:
    y = np.arange(5 * 8).reshape((5, 8))
    assert_eq(da.diag(y, k), np.diag(y, k))
    # when input 2d-array is a rectangular dask array:
    d = da.from_array(y, chunks=(4, 4))
    assert_eq(da.diag(d, k), np.diag(y, k))
    # heterogeneous chunks:
    d = da.from_array(y, chunks=((3, 2), (4, 1, 2, 1)))
    assert_eq(da.diag(d, k), np.diag(y, k))


def test_creation_data_producers():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))
    dsk = collections_to_expr([d]).__dask_graph__()
    assert all(v.data_producer for v in dsk.values())

    # blockwise fusion
    x = d.astype("float64")
    dsk = collections_to_expr([x]).__dask_graph__()
    assert sum(v.data_producer for v in dsk.values()) == 4
    assert sum(isinstance(v, Alias) for v in dsk.values()) == 4
    assert len(dsk) == 8

    # linear fusion
    x = d[slice(0, 6), None].astype("float64")
    dsk = collections_to_expr([x]).__dask_graph__()
    assert sum(v.data_producer for v in dsk.values()) == 4
    assert sum(isinstance(v, Alias) for v in dsk.values()) == 4
    assert len(dsk) == 8

    # no fusion
    x = d[[1, 3, 5, 7, 6, 4, 2, 0]].astype("float64")
    dsk = collections_to_expr([x]).__dask_graph__()
    assert sum(v.data_producer and "array-" in k[0] for k, v in dsk.items()) == 4
    assert sum(v.data_producer for v in dsk.values()) == 8  # getitem data nodes
    assert len(dsk) == 24


def test_diagonal():
    v = np.arange(11)
    with pytest.raises(ValueError):
        da.diagonal(v)

    v = np.arange(4).reshape((2, 2))
    with pytest.raises(ValueError):
        da.diagonal(v, axis1=0, axis2=0)

    with pytest.raises(AxisError):
        da.diagonal(v, axis1=-4)

    with pytest.raises(AxisError):
        da.diagonal(v, axis2=-4)

    v = np.arange(4 * 5 * 6).reshape((4, 5, 6))
    v = da.from_array(v, chunks=2)
    assert_eq(da.diagonal(v), np.diagonal(v))
    # Empty diagonal.
    assert_eq(da.diagonal(v, offset=10), np.diagonal(v, offset=10))
    assert_eq(da.diagonal(v, offset=-10), np.diagonal(v, offset=-10))

    with pytest.raises(ValueError):
        da.diagonal(v, axis1=-2)

    # Negative axis.
    assert_eq(da.diagonal(v, axis1=-1), np.diagonal(v, axis1=-1))
    assert_eq(da.diagonal(v, offset=1, axis1=-1), np.diagonal(v, offset=1, axis1=-1))

    # Heterogeneous chunks.
    v = np.arange(2 * 3 * 4 * 5 * 6).reshape((2, 3, 4, 5, 6))
    v = da.from_array(v, chunks=(1, (1, 2), (1, 2, 1), (2, 1, 2), (5, 1)))

    assert_eq(da.diagonal(v), np.diagonal(v))
    assert_eq(
        da.diagonal(v, offset=2, axis1=3, axis2=1),
        np.diagonal(v, offset=2, axis1=3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-2, axis1=3, axis2=1),
        np.diagonal(v, offset=-2, axis1=3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-2, axis1=3, axis2=4),
        np.diagonal(v, offset=-2, axis1=3, axis2=4),
    )

    assert_eq(da.diagonal(v, 1), np.diagonal(v, 1))
    assert_eq(da.diagonal(v, -1), np.diagonal(v, -1))
    # Positional arguments
    assert_eq(da.diagonal(v, 1, 2, 1), np.diagonal(v, 1, 2, 1))

    v = np.arange(2 * 3 * 4 * 5 * 6).reshape((2, 3, 4, 5, 6))
    assert_eq(da.diagonal(v, axis1=1, axis2=3), np.diagonal(v, axis1=1, axis2=3))
    assert_eq(
        da.diagonal(v, offset=1, axis1=1, axis2=3),
        np.diagonal(v, offset=1, axis1=1, axis2=3),
    )

    assert_eq(
        da.diagonal(v, offset=1, axis1=3, axis2=1),
        np.diagonal(v, offset=1, axis1=3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-5, axis1=3, axis2=1),
        np.diagonal(v, offset=-5, axis1=3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-6, axis1=3, axis2=1),
        np.diagonal(v, offset=-6, axis1=3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-6, axis1=-3, axis2=1),
        np.diagonal(v, offset=-6, axis1=-3, axis2=1),
    )

    assert_eq(
        da.diagonal(v, offset=-6, axis1=-3, axis2=1),
        np.diagonal(v, offset=-6, axis1=-3, axis2=1),
    )

    v = da.from_array(v, chunks=2)
    assert_eq(
        da.diagonal(v, offset=1, axis1=3, axis2=1),
        np.diagonal(v, offset=1, axis1=3, axis2=1),
    )
    assert_eq(
        da.diagonal(v, offset=-1, axis1=3, axis2=1),
        np.diagonal(v, offset=-1, axis1=3, axis2=1),
    )

    v = np.arange(384).reshape((8, 8, 6))
    assert_eq(da.diagonal(v, offset=-1, axis1=2), np.diagonal(v, offset=-1, axis1=2))

    v = da.from_array(v, chunks=(4, 4, 2))
    assert_eq(da.diagonal(v, offset=-1, axis1=2), np.diagonal(v, offset=-1, axis1=2))


@pytest.mark.parametrize("dtype", [None, "f8", "i8"])
@pytest.mark.parametrize(
    "func, kwargs",
    [
        (lambda x, y: x + y, {}),
        (lambda x, y, c=1: x + c * y, {}),
        (lambda x, y, c=1: x + c * y, {"c": 3}),
    ],
)
def test_fromfunction(func, dtype, kwargs):
    a = np.fromfunction(func, shape=(5, 5), dtype=dtype, **kwargs)
    d = da.fromfunction(func, shape=(5, 5), chunks=(2, 2), dtype=dtype, **kwargs)

    assert_eq(d, a)

    d2 = da.fromfunction(func, shape=(5, 5), chunks=(2, 2), dtype=dtype, **kwargs)

    assert same_keys(d, d2)


def test_repeat():
    x = np.random.random((10, 11, 13))
    d = da.from_array(x, chunks=(4, 5, 3))

    repeats = [0, 1, 2, 5]
    axes = [-3, -2, -1, 0, 1, 2]

    for r in repeats:
        for a in axes:
            assert_eq(x.repeat(r, axis=a), d.repeat(r, axis=a))

    assert_eq(d.repeat(2, 0), da.repeat(d, 2, 0))

    with pytest.raises(NotImplementedError):
        da.repeat(d, np.arange(10))

    with pytest.raises(NotImplementedError):
        da.repeat(d, 2, None)

    with pytest.raises(NotImplementedError):
        da.repeat(d, 2)

    for invalid_axis in [3, -4]:
        with pytest.raises(ValueError):
            da.repeat(d, 2, axis=invalid_axis)

    x = np.arange(5)
    d = da.arange(5, chunks=(2,))

    assert_eq(x.repeat(3), d.repeat(3))

    for r in [1, 2, 3, 4]:
        assert all(concat(d.repeat(r).chunks))


@pytest.mark.parametrize("reps", [2, (2, 2), (1, 2), (2, 1), (2, 3, 4, 0)])
def test_tile_basic(reps):
    a = da.asarray([0, 1, 2])
    b = [[1, 2], [3, 4]]

    assert_eq(np.tile(a.compute(), reps), da.tile(a, reps))
    assert_eq(np.tile(b, reps), da.tile(b, reps))


@pytest.mark.parametrize("shape, chunks", [((10,), (1,)), ((10, 11, 13), (4, 5, 3))])
@pytest.mark.parametrize("reps", [0, 1, 2, 3, 5, (1,), (1, 2)])
def test_tile_chunks(shape, chunks, reps):
    x = np.random.random(shape)
    d = da.from_array(x, chunks=chunks)

    assert_eq(np.tile(x, reps), da.tile(d, reps))


@pytest.mark.parametrize("shape, chunks", [((10,), (1,)), ((10, 11, 13), (4, 5, 3))])
@pytest.mark.parametrize("reps", [-1, -5])
def test_tile_neg_reps(shape, chunks, reps):
    x = np.random.random(shape)
    d = da.from_array(x, chunks=chunks)

    with pytest.raises(ValueError):
        da.tile(d, reps)


@pytest.mark.parametrize("shape, chunks", [((10,), (1,)), ((10, 11, 13), (4, 5, 3))])
@pytest.mark.parametrize("reps", [0, (0,), (2, 0), (0, 3, 0, 4)])
def test_tile_zero_reps(shape, chunks, reps):
    x = np.random.random(shape)
    d = da.from_array(x, chunks=chunks)

    assert_eq(np.tile(x, reps), da.tile(d, reps))


@pytest.mark.parametrize("shape, chunks", [((1, 1, 0), (1, 1, 0)), ((2, 0), (1, 0))])
@pytest.mark.parametrize("reps", [2, (3, 2, 5)])
def test_tile_empty_array(shape, chunks, reps):
    x = np.empty(shape)
    d = da.from_array(x, chunks=chunks)

    assert_eq(np.tile(x, reps), da.tile(d, reps))


@pytest.mark.parametrize(
    "shape", [(3,), (2, 3), (3, 4, 3), (3, 2, 3), (4, 3, 2, 4), (2, 2)]
)
@pytest.mark.parametrize("reps", [(2,), (1, 2), (2, 1), (2, 2), (2, 3, 2), (3, 2)])
def test_tile_np_kroncompare_examples(shape, reps):
    x = np.random.random(shape)
    d = da.asarray(x)

    assert_eq(np.tile(x, reps), da.tile(d, reps))


@pytest.mark.parametrize(
    "shape, chunks, pad_width, mode, kwargs",
    [
        ((10, 11), (4, 5), 0, "constant", {"constant_values": 2}),
        ((10, 11), (4, 5), 0, "edge", {}),
        ((10, 11), (4, 5), 0, "linear_ramp", {"end_values": 2}),
        ((10, 11), (4, 5), 0, "reflect", {}),
        ((10, 11), (4, 5), 0, "symmetric", {}),
        ((10, 11), (4, 5), 0, "wrap", {}),
        ((10, 11), (4, 5), 0, "empty", {}),
    ],
)
def test_pad_0_width(shape, chunks, pad_width, mode, kwargs):
    np_a = np.random.random(shape)
    da_a = da.from_array(np_a, chunks=chunks)

    np_r = np.pad(np_a, pad_width, mode, **kwargs)
    da_r = da.pad(da_a, pad_width, mode, **kwargs)

    assert da_r is da_a

    assert_eq(np_r, da_r)


@pytest.mark.parametrize(
    "shape, chunks, pad_width, mode, kwargs",
    [
        ((10,), (3,), 1, "constant", {}),
        ((10,), (3,), 2, "constant", {"constant_values": -1}),
        ((10,), (3,), 2, "constant", {"constant_values": np.array(-1)}),
        ((10,), (3,), ((2, 3)), "constant", {"constant_values": (-1, -2)}),
        (
            (10, 11),
            (4, 5),
            ((1, 4), (2, 3)),
            "constant",
            {"constant_values": ((-1, -2), (2, 1))},
        ),
        ((10,), (3,), 3, "edge", {}),
        ((10,), (3,), 3, "linear_ramp", {}),
        ((10,), (3,), 3, "linear_ramp", {"end_values": 0}),
        (
            (10, 11),
            (4, 5),
            ((1, 4), (2, 3)),
            "linear_ramp",
            {"end_values": ((-1, -2), (4, 3))},
        ),
        ((10, 11), (4, 5), ((1, 4), (2, 3)), "reflect", {}),
        ((10, 11), (4, 5), ((1, 4), (2, 3)), "symmetric", {}),
        ((10, 11), (4, 5), ((1, 4), (2, 3)), "wrap", {}),
        ((10,), (3,), ((2, 3)), "maximum", {"stat_length": (1, 2)}),
        ((10, 11), (4, 5), ((1, 4), (2, 3)), "mean", {"stat_length": ((3, 4), (2, 1))}),
        ((10,), (3,), ((2, 3)), "minimum", {"stat_length": (2, 3)}),
        ((10,), (3,), 1, "empty", {}),
    ],
)
def test_pad(shape, chunks, pad_width, mode, kwargs):
    np_a = np.random.random(shape)
    da_a = da.from_array(np_a, chunks=chunks)

    np_r = np.pad(np_a, pad_width, mode, **kwargs)
    da_r = da.pad(da_a, pad_width, mode, **kwargs)

    if mode == "empty":
        # empty pads lead to undefined values which may be different
        assert_eq(np_r[pad_width:-pad_width], da_r[pad_width:-pad_width])
    else:
        assert_eq(np_r, da_r)


@pytest.mark.parametrize(
    ["np_a", "pad_value"],
    (
        (np.arange(4, dtype="int64"), np.int64(1)),
        (np.arange(4, dtype="float64"), np.float64(0)),
        (
            np.array(
                ["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04"],
                dtype="datetime64[ns]",
            ),
            np.datetime64("1972-01-01"),
        ),
        (np.array([True, False, True, True], dtype=np.bool_), np.bool_(False)),
        (np.array(["ab", "bc", "de", "ef"], dtype=np.str_), np.str_("00")),
        (np.arange(4, dtype="int64"), np.array(1, dtype="int64")),
        (np.arange(4, dtype="float64"), np.array(0, dtype="float64")),
        (
            np.array(
                ["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04"],
                dtype="datetime64[ns]",
            ),
            np.array("1972-01-01", dtype="datetime64[ns]"),
        ),
        (
            np.array([True, False, True, True], dtype=np.bool_),
            np.array(False, dtype=np.bool_),
        ),
        (
            np.array(["ab", "bc", "de", "ef"], dtype=np.str_),
            np.array("00", dtype=np.str_),
        ),
    ),
)
def test_pad_constant_values(np_a, pad_value):
    pad_width = (1, 1)
    da_a = da.from_array(np_a, chunks=(2,))

    np_r = np.pad(np_a, pad_width, mode="constant", constant_values=pad_value)
    da_r = da.pad(da_a, pad_width, mode="constant", constant_values=pad_value)

    assert_eq(np_r, da_r)


@pytest.mark.parametrize("dtype", [np.uint8, np.int16, np.float32, bool])
@pytest.mark.parametrize(
    "pad_widths", [2, (2,), (2, 3), ((2, 3),), ((3, 1), (0, 0), (2, 0))]
)
@pytest.mark.parametrize(
    "mode",
    [
        "constant",
        "edge",
        "linear_ramp",
        "maximum",
        "mean",
        "minimum",
        pytest.param(
            "reflect",
            marks=pytest.mark.skip(
                reason="Bug when pad_width is larger than dimension: https://github.com/dask/dask/issues/5303"
            ),
        ),
        pytest.param(
            "symmetric",
            marks=pytest.mark.skip(
                reason="Bug when pad_width is larger than dimension: https://github.com/dask/dask/issues/5303"
            ),
        ),
        pytest.param(
            "wrap",
            marks=pytest.mark.skip(
                reason="Bug when pad_width is larger than dimension: https://github.com/dask/dask/issues/5303"
            ),
        ),
        pytest.param(
            "median",
            marks=pytest.mark.skip(reason="Not implemented"),
        ),
        pytest.param(
            "empty",
            marks=pytest.mark.skip(
                reason="Empty leads to undefined values, which may be different"
            ),
        ),
    ],
)
def test_pad_3d_data(dtype, pad_widths, mode):
    np_a = np.arange(2 * 3 * 4).reshape(2, 3, 4).astype(dtype)
    da_a = da.from_array(np_a, chunks="auto")

    np_r = np.pad(np_a, pad_widths, mode=mode)
    da_r = da.pad(da_a, pad_widths, mode=mode)

    assert_eq(np_r, da_r)


@pytest.mark.parametrize("kwargs", [{}, {"scaler": 2}])
def test_pad_udf(kwargs):
    def udf_pad(vector, pad_width, iaxis, inner_kwargs):
        assert kwargs == inner_kwargs
        scaler = inner_kwargs.get("scaler", 1)
        vector[: pad_width[0]] = -scaler * pad_width[0]
        vector[-pad_width[1] :] = scaler * pad_width[1]
        return vector

    shape = (10, 11)
    chunks = (4, 5)
    pad_width = ((1, 2), (2, 3))

    np_a = np.random.random(shape)
    da_a = da.from_array(np_a, chunks=chunks)

    np_r = np.pad(np_a, pad_width, udf_pad, **kwargs)
    da_r = da.pad(da_a, pad_width, udf_pad, **kwargs)

    assert_eq(np_r, da_r)


def test_pad_constant_chunksizes():
    array = dask.array.ones((10, 10), chunks=(1, 1))
    result = dask.array.pad(
        array, ((0, 16 - 10), (0, 0)), mode="constant", constant_values=0
    )
    assert tuple(map(max, result.chunks)) == (1, 1)
    assert_eq(
        result,
        np.pad(
            array.compute(),
            mode="constant",
            constant_values=0,
            pad_width=((0, 16 - 10), (0, 0)),
        ),
    )


def test_auto_chunks():
    with dask.config.set({"array.chunk-size": "50 MiB"}):
        x = da.ones((10000, 10000))
        assert 4 < x.npartitions < 32


def test_string_auto_chunk():
    with pytest.raises(ValueError):
        da.full((10000, 10000), "auto_chunk", chunks="auto")


def test_diagonal_zero_chunks():
    x = da.ones((8, 8), chunks=(4, 4))
    dd = da.ones((8, 8), chunks=(4, 4))
    d = da.diagonal(dd)

    expected = np.ones((8,))
    assert_eq(d, expected)
    assert_eq(d + d, 2 * expected)
    A = d + x
    assert_eq(A, np.full((8, 8), 2.0))


@pytest.mark.parametrize("fn", ["zeros_like", "ones_like"])
@pytest.mark.parametrize("shape_chunks", [((50, 4), (10, 2)), ((50,), (10,))])
@pytest.mark.parametrize("dtype", ["u4", np.float32, None, np.int64])
def test_nan_zeros_ones_like(fn, shape_chunks, dtype):
    dafn = getattr(da, fn)
    npfn = getattr(np, fn)
    shape, chunks = shape_chunks
    x1 = da.random.standard_normal(size=shape, chunks=chunks)
    y1 = x1[x1 < 0.5]
    x2 = x1.compute()
    y2 = x2[x2 < 0.5]
    assert_eq(
        dafn(y1, dtype=dtype),
        npfn(y2, dtype=dtype),
    )


def test_from_array_getitem_fused():
    arr = np.arange(100).reshape(10, 10)
    darr = da.from_array(arr, chunks=(5, 5))
    result = darr[slice(1, 5), :][slice(1, 3), :]
    dsk = collections_to_expr([result]).__dask_graph__()
    # Ensure that slices are merged properly
    key = [k for k in dsk if "array-getitem" in k[0]][0]
    key_2 = [
        k
        for k, v in dsk[key].args[0].items()
        if "getitem" in k[0] and not isinstance(v, Alias)
    ][0]
    assert dsk[key].args[0][key_2].args[1] == ((slice(2, 4), slice(0, None)))
    assert_eq(result, arr[slice(1, 5), :][slice(1, 3), :])


@pytest.mark.parametrize("shape_chunks", [((50, 4), (10, 2)), ((50,), (10,))])
@pytest.mark.parametrize("dtype", ["u4", np.float32, None, np.int64])
def test_nan_empty_like(shape_chunks, dtype):
    shape, chunks = shape_chunks
    x1 = da.random.standard_normal(size=shape, chunks=chunks)
    y1 = x1[x1 < 0.5]
    x2 = x1.compute()
    y2 = x2[x2 < 0.5]
    a_da = da.empty_like(y1, dtype=dtype).compute()
    a_np = np.empty_like(y2, dtype=dtype)
    assert a_da.shape == a_np.shape
    assert a_da.dtype == a_np.dtype


@pytest.mark.parametrize("val", [0, 0.0, 99, -1])
@pytest.mark.parametrize("shape_chunks", [((50, 4), (10, 2)), ((50,), (10,))])
@pytest.mark.parametrize("dtype", ["u4", np.float32, None, np.int64])
def test_nan_full_like(val, shape_chunks, dtype):
    if NUMPY_GE_210 and val == -1 and dtype == "u4":
        pytest.xfail("can't insert negative numbers into unsigned integer")
    shape, chunks = shape_chunks
    x1 = da.random.standard_normal(size=shape, chunks=chunks)
    y1 = x1[x1 < 0.5]
    x2 = x1.compute()
    y2 = x2[x2 < 0.5]
    assert_eq(
        da.full_like(y1, val, dtype=dtype),
        np.full_like(y2, val, dtype=dtype),
    )


@pytest.mark.parametrize(
    "func", [da.array, da.asarray, da.asanyarray, da.arange, da.tri]
)
def test_like_forgets_graph(func):
    """Test that array creation functions with like=x do not
    internally store the graph of x
    """
    x = da.arange(3).map_blocks(lambda x: x)
    with pytest.raises(Exception, match="local object"):
        pickle.dumps(x)

    a = func(1, like=x)
    pickle.dumps(a)
