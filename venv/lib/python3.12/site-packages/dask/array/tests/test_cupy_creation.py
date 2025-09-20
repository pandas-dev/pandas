from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import dask.array as da
from dask import config
from dask.array.numpy_compat import AxisError
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")


def test_diag():
    v = cupy.arange(11)
    dv = da.from_array(v, chunks=(4,), asarray=False)
    assert type(dv._meta) == cupy.ndarray
    assert_eq(dv, dv)  # Check that _meta and computed arrays match types
    assert_eq(da.diag(dv), cupy.diag(v))

    v = v + v + 3
    dv = dv + dv + 3
    darr = da.diag(dv)
    cupyarr = cupy.diag(v)
    assert type(darr._meta) == cupy.ndarray
    assert_eq(darr, darr)  # Check that _meta and computed arrays match types
    assert_eq(darr, cupyarr)

    x = cupy.arange(64).reshape((8, 8))
    dx = da.from_array(x, chunks=(4, 4), asarray=False)
    assert type(dx._meta) == cupy.ndarray
    assert_eq(dx, dx)  # Check that _meta and computed arrays match types
    assert_eq(da.diag(dx), cupy.diag(x))


def test_diagonal():
    v = cupy.arange(11)
    with pytest.raises(ValueError):
        da.diagonal(v)

    v = cupy.arange(4).reshape((2, 2))
    with pytest.raises(ValueError):
        da.diagonal(v, axis1=0, axis2=0)

    with pytest.raises(AxisError):
        da.diagonal(v, axis1=-4)

    with pytest.raises(AxisError):
        da.diagonal(v, axis2=-4)

    v = cupy.arange(4 * 5 * 6).reshape((4, 5, 6))
    v = da.from_array(v, chunks=2, asarray=False)
    assert_eq(da.diagonal(v), np.diagonal(v))
    # Empty diagonal.
    assert_eq(da.diagonal(v, offset=10), np.diagonal(v, offset=10))
    assert_eq(da.diagonal(v, offset=-10), np.diagonal(v, offset=-10))
    assert isinstance(da.diagonal(v).compute(), cupy.ndarray)

    with pytest.raises(ValueError):
        da.diagonal(v, axis1=-2)

    # Negative axis.
    assert_eq(da.diagonal(v, axis1=-1), np.diagonal(v, axis1=-1))
    assert_eq(da.diagonal(v, offset=1, axis1=-1), np.diagonal(v, offset=1, axis1=-1))

    # Heterogeneous chunks.
    v = cupy.arange(2 * 3 * 4 * 5 * 6).reshape((2, 3, 4, 5, 6))
    v = da.from_array(
        v, chunks=(1, (1, 2), (1, 2, 1), (2, 1, 2), (5, 1)), asarray=False
    )

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


@pytest.mark.parametrize(
    "shape, chunks, pad_width, mode, kwargs",
    [
        ((10,), (3,), 1, "constant", {}),
        ((10,), (3,), 2, "constant", {"constant_values": -1}),
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
    np_a = np.random.default_rng().random(shape)
    da_a = da.from_array(cupy.array(np_a), chunks=chunks)

    np_r = np.pad(np_a, pad_width, mode, **kwargs)
    da_r = da.pad(da_a, pad_width, mode, **kwargs)

    assert isinstance(da_r._meta, cupy.ndarray)
    assert isinstance(da_r.compute(), cupy.ndarray)

    if mode == "empty":
        # empty pads lead to undefined values which may be different
        assert_eq(
            np_r[pad_width:-pad_width], da_r[pad_width:-pad_width], check_type=False
        )
    else:
        assert_eq(np_r, da_r, check_type=False)


@pytest.mark.parametrize("xp", [np, da])
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
def test_tri_like(xp, N, M, k, dtype, chunks):
    args = [N, M, k, dtype]

    cp_a = cupy.tri(*args)

    if xp is da:
        args.append(chunks)
    xp_a = xp.tri(*args, like=da.from_array(cupy.array(())))

    assert_eq(xp_a, cp_a)


def test_to_backend_cupy():
    # Test that `Array.to_backend` works as expected
    with config.set({"array.backend": "numpy"}):
        # Start with cupy-backed array
        x = da.from_array(cupy.arange(11), chunks=(4,))
        assert isinstance(x._meta, cupy.ndarray)

        # Calling default `to_backend` should move
        # backend to `numpy`
        x_new = x.to_backend()
        assert isinstance(x_new._meta, np.ndarray)

        # Calling `to_backend("cudf")` should always
        # move the data back to `cupy`
        x_new = x.to_backend("cupy")
        assert isinstance(x_new._meta, cupy.ndarray)

        # Change global "array.backend" config to `cupy`
        with config.set({"array.backend": "cupy"}):
            # Calling `to_backend("numpy")` should
            # always move the data to `numpy`
            x_new = x.to_backend("numpy")
            assert isinstance(x_new._meta, np.ndarray)

            # Calling default `to_backend()` should
            # satisfy the global config (now `cupy`)
            x_new = x.to_backend()
            assert isinstance(x_new._meta, cupy.ndarray)

    assert_eq(x, x.to_backend("numpy"), check_type=False)
    assert_eq(x, x.to_backend("cupy"))
