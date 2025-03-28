from __future__ import annotations

import contextlib
import itertools
import pickle
import sys
import warnings
from numbers import Number

import pytest

import dask
from dask.delayed import delayed

np = pytest.importorskip("numpy")

import dask.array as da
from dask.array.numpy_compat import NUMPY_GE_200, AxisError
from dask.array.utils import assert_eq, same_keys

if da._array_expr_enabled():
    pytest.skip("parametrize using unsupported functions", allow_module_level=True)


def test_array():
    x = np.ones(5, dtype="i4")
    d = da.ones(5, chunks=3, dtype="i4")
    assert_eq(da.array(d, ndmin=3, dtype="i8"), np.array(x, ndmin=3, dtype="i8"))

    # regression #1847 this shall not raise an exception.
    x = da.ones((100, 3), chunks=10)
    y = da.array(x)
    assert isinstance(y, da.Array)


def test_array_return_type():
    # Regression test for https://github.com/dask/dask/issues/5426
    x = [0, 1, 2, 3]
    dx = da.array(x)
    assert isinstance(dx, da.Array)
    assert_eq(x, dx)


def test_derived_docstrings():
    assert "This docstring was copied from numpy.array" in da.routines.array.__doc__
    assert "Create an array." in da.routines.array.__doc__


@pytest.mark.parametrize("funcname", ["atleast_1d", "atleast_2d", "atleast_3d"])
def test_atleast_nd_no_args(funcname):
    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r_n = np_func()
    da_r_n = da_func()

    assert np_r_n == da_r_n


@pytest.mark.parametrize("funcname", ["atleast_1d", "atleast_2d", "atleast_3d"])
@pytest.mark.parametrize(
    "shape, chunks",
    [
        (tuple(), tuple()),
        ((4,), (2,)),
        ((4, 6), (2, 3)),
        ((4, 6, 8), (2, 3, 4)),
        ((4, 6, 8, 10), (2, 3, 4, 5)),
    ],
)
def test_atleast_nd_one_arg(funcname, shape, chunks):
    np_a = np.random.default_rng().random(shape)
    da_a = da.from_array(np_a, chunks=chunks)

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r = np_func(np_a)
    da_r = da_func(da_a)

    assert_eq(np_r, da_r)


@pytest.mark.parametrize("funcname", ["atleast_1d", "atleast_2d", "atleast_3d"])
@pytest.mark.parametrize(
    "shape1, shape2",
    list(
        itertools.combinations_with_replacement(
            [tuple(), (4,), (4, 6), (4, 6, 8), (4, 6, 8, 10)], 2
        )
    ),
)
def test_atleast_nd_two_args(funcname, shape1, shape2):
    np_a_1 = np.random.default_rng().random(shape1)
    da_a_1 = da.from_array(np_a_1, chunks=tuple(c // 2 for c in shape1))

    np_a_2 = np.random.default_rng().random(shape2)
    da_a_2 = da.from_array(np_a_2, chunks=tuple(c // 2 for c in shape2))

    np_a_n = [np_a_1, np_a_2]
    da_a_n = [da_a_1, da_a_2]

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r_n = np_func(*np_a_n)
    da_r_n = da_func(*da_a_n)

    assert type(np_r_n) is type(da_r_n)

    assert len(np_r_n) == len(da_r_n)

    for np_r, da_r in zip(np_r_n, da_r_n):
        assert_eq(np_r, da_r)


def test_transpose():
    x = np.arange(240).reshape((4, 6, 10))
    d = da.from_array(x, (2, 3, 4))

    assert_eq(d.transpose((2, 0, 1)), x.transpose((2, 0, 1)))
    assert same_keys(d.transpose((2, 0, 1)), d.transpose((2, 0, 1)))

    assert_eq(d.transpose(2, 0, 1), x.transpose(2, 0, 1))
    assert same_keys(d.transpose(2, 0, 1), d.transpose(2, 0, 1))

    with pytest.raises(ValueError):
        d.transpose(1, 2)

    with pytest.raises(ValueError):
        d.transpose((1, 2))


def test_transpose_negative_axes():
    x = np.ones((2, 3, 4, 5))
    y = da.ones((2, 3, 4, 5), chunks=3)

    assert_eq(x.transpose([-1, -2, 0, 1]), y.transpose([-1, -2, 0, 1]))


def test_transpose_skip_when_possible():
    x = da.ones((2, 3, 4), chunks=3)
    assert x.transpose((0, 1, 2)) is x
    assert x.transpose((-3, -2, -1)) is x


def test_swapaxes():
    x = np.random.default_rng().normal(0, 10, size=(10, 12, 7))
    d = da.from_array(x, chunks=(4, 5, 2))

    assert_eq(np.swapaxes(x, 0, 1), da.swapaxes(d, 0, 1))
    assert_eq(np.swapaxes(x, 2, 1), da.swapaxes(d, 2, 1))
    assert_eq(x.swapaxes(2, 1), d.swapaxes(2, 1))
    assert_eq(x.swapaxes(0, 0), d.swapaxes(0, 0))
    assert_eq(x.swapaxes(1, 2), d.swapaxes(1, 2))
    assert_eq(x.swapaxes(0, -1), d.swapaxes(0, -1))
    assert_eq(x.swapaxes(-1, 1), d.swapaxes(-1, 1))

    assert d.swapaxes(0, 1).name == d.swapaxes(0, 1).name
    assert d.swapaxes(0, 1).name != d.swapaxes(1, 0).name


@pytest.mark.parametrize("funcname", ["moveaxis", "rollaxis"])
@pytest.mark.parametrize("shape", [(), (5,), (3, 5, 7, 3)])
def test_moveaxis_rollaxis(funcname, shape):
    x = np.random.default_rng().random(shape)
    d = da.from_array(x, chunks=(len(shape) * (2,)))
    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)
    for axis1 in range(-x.ndim, x.ndim):
        assert isinstance(da_func(d, 0, axis1), da.Array)
        for axis2 in range(-x.ndim, x.ndim):
            assert_eq(np_func(x, axis1, axis2), da_func(d, axis1, axis2))


def test_moveaxis_rollaxis_keyword():
    x = np.random.default_rng().random((10, 12, 7))
    d = da.from_array(x, chunks=(4, 5, 2))
    assert_eq(
        np.moveaxis(x, destination=1, source=0), da.moveaxis(d, destination=1, source=0)
    )
    assert_eq(np.rollaxis(x, 2), da.rollaxis(d, 2))
    assert isinstance(da.rollaxis(d, 1), da.Array)
    assert_eq(np.rollaxis(x, start=1, axis=2), da.rollaxis(d, start=1, axis=2))


def test_moveaxis_rollaxis_numpy_api():
    a = da.random.default_rng().random((4, 4, 4), chunks=2)
    result = np.moveaxis(a, 2, 0)
    assert isinstance(result, da.Array)
    assert_eq(result, np.moveaxis(a.compute(), 2, 0))

    result = np.rollaxis(a, 2, 0)
    assert isinstance(result, da.Array)
    assert_eq(result, np.rollaxis(a.compute(), 2, 0))


@pytest.mark.parametrize(
    "funcname, kwargs",
    [
        ("flipud", {}),
        ("fliplr", {}),
        ("flip", {}),
        ("flip", {"axis": 0}),
        ("flip", {"axis": 1}),
        ("flip", {"axis": 2}),
        ("flip", {"axis": -1}),
        ("flip", {"axis": (0, 2)}),
    ],
)
@pytest.mark.parametrize("shape", [tuple(), (4,), (4, 6), (4, 6, 8), (4, 6, 8, 10)])
def test_flip(funcname, kwargs, shape):
    axis = kwargs.get("axis")
    if axis is None:
        if funcname == "flipud":
            axis = (0,)
        elif funcname == "fliplr":
            axis = (1,)
        elif funcname == "flip":
            axis = range(len(shape))
    elif not isinstance(axis, tuple):
        axis = (axis,)

    np_a = np.random.default_rng().random(shape)
    da_a = da.from_array(np_a, chunks=1)

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    try:
        for ax in axis:
            range(np_a.ndim)[ax]
    except IndexError:
        with pytest.raises(ValueError):
            da_func(da_a, **kwargs)
    else:
        np_r = np_func(np_a, **kwargs)
        da_r = da_func(da_a, **kwargs)

        assert_eq(np_r, da_r)


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"axes": (1, 0)},
        {"axes": (2, 3)},
        {"axes": (0, 1, 2)},
    ],
)
@pytest.mark.parametrize(
    "shape",
    [
        tuple(),
        (4,),
        (4, 6),
        (4, 6, 8),
    ],
)
def test_rot90(kwargs, shape):
    axes = kwargs.get("axes", (0, 1))
    np_a = np.random.default_rng().random(shape)
    da_a = da.from_array(np_a, chunks=2)

    np_func = np.rot90
    da_func = da.rot90

    try:
        for axis in axes[:2]:
            range(np_a.ndim)[axis]
    except IndexError:
        with pytest.raises(ValueError):
            da_func(da_a, **kwargs)
    else:
        if len(axes) != 2 or axes[0] == axes[1]:
            with pytest.raises(ValueError):
                da_func(da_a, **kwargs)
        else:
            for k in range(-3, 9):
                np_r = np_func(np_a, k=k, **kwargs)
                da_r = da_func(da_a, k=k, **kwargs)
                assert_eq(np_r, da_r)


@pytest.mark.parametrize(
    "x_shape, y_shape, x_chunks, y_chunks",
    [
        [(), (), (), ()],
        [(), (7,), (), ()],
        [(), (7, 11), (), ()],
        [(), (7, 11, 15), (), ()],
        [(), (7, 11, 15, 19), (), ()],
        [(7,), (), (), ()],
        [(7,), (7,), (), ()],
        [(11,), (11, 7), (), ()],
        [(15,), (7, 15, 11), (), ()],
        [(19,), (7, 11, 19, 15), (), ()],
        [(7, 11), (), (), ()],
        [(7, 11), (11,), (), ()],
        [(7, 11), (11, 7), (), ()],
        [(11, 15), (7, 15, 11), (), ()],
        [(15, 19), (7, 11, 19, 15), (), ()],
        [(7, 11, 15), (), (), ()],
        [(7, 11, 15), (15,), (), ()],
        [(7, 11, 15), (15, 7), (), ()],
        [(7, 11, 15), (7, 15, 11), (), ()],
        [(11, 15, 19), (7, 11, 19, 15), (), ()],
        [(7, 11, 15, 19), (), (), ()],
        [(7, 11, 15, 19), (19,), (), ()],
        [(7, 11, 15, 19), (19, 7), (), ()],
        [(7, 11, 15, 19), (11, 19, 13), (), ()],
        [(7, 11, 15, 19), (7, 11, 19, 15), (), ()],
        # These tests use explicitly special/disparate chunk sizes:
        [(), (7,), (), (5,)],
        [(), (7, 11, 15, 19), (), (1, 3, 5, 19)],
        [(7, 11), (11, 7), (1, 1), (1, 1)],
        [(7, 11), (11, 7), (3, 5), (4, 2)],
        [(7, 11), (11, 7), (7, 11), (11, 7)],
        [(11, 15, 19), (7, 11, 19, 15), (7, 7, 7), (3, 9, 9, 9)],
        [(3, 3, 20, 30), (3, 3, 30, 20), (1, 3, 2, 6), (1, 3, 5, 10)],
    ],
)
def test_matmul(x_shape, y_shape, x_chunks, y_chunks):
    rng = np.random.default_rng(3732)

    x = rng.random(x_shape)[()]
    y = rng.random(y_shape)[()]

    a = da.from_array(x, chunks=x_chunks or tuple((i // 2) for i in x.shape))
    b = da.from_array(y, chunks=y_chunks or tuple((i // 2) for i in y.shape))

    expected = None
    try:
        expected = np.matmul(x, y)
    except ValueError:
        pass

    for d1, d2 in itertools.product([a, x], [b, y]):
        if x.ndim == 0 or y.ndim == 0:
            with pytest.raises(ValueError):
                da.matmul(d1, d2)
        else:
            assert_eq(expected, da.matmul(d1, d2))


def test_tensordot():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 4))
    y = np.arange(200).reshape((20, 10))
    b = da.from_array(y, chunks=(4, 5))

    for axes in [1, (1, 0), (-1, 0)]:
        assert_eq(da.tensordot(a, b, axes=axes), np.tensordot(x, y, axes=axes))
        assert_eq(da.tensordot(x, b, axes=axes), np.tensordot(x, y, axes=axes))
        assert_eq(da.tensordot(a, y, axes=axes), np.tensordot(x, y, axes=axes))

    assert same_keys(da.tensordot(a, b, axes=(1, 0)), da.tensordot(a, b, axes=(1, 0)))

    # Increasing number of chunks warning
    with pytest.warns(da.PerformanceWarning):
        assert not same_keys(da.tensordot(a, b, axes=0), da.tensordot(a, b, axes=1))


@pytest.mark.parametrize(
    "axes", [0, 1, (0, 1), (1, 0), ((1, 0), (2, 1)), ((1, 2), (2, 0)), ((2, 0), (1, 2))]
)
def test_tensordot_2(axes):
    x = np.arange(4 * 4 * 4).reshape((4, 4, 4))
    y = da.from_array(x, chunks=2)

    assert_eq(da.tensordot(y, y, axes=axes), np.tensordot(x, x, axes=axes))


@pytest.mark.parametrize("chunks", ["auto", (4, 6), (2, 3), (4, 3), (2, 6)])
def test_tensordot_double_contraction_neq2(chunks):
    # Regression test for https://github.com/dask/dask/issues/5472
    x = np.arange(24).reshape(4, 6)
    y = da.from_array(x, chunks=chunks)
    assert_eq(da.tensordot(y, y, axes=2), np.tensordot(x, x, axes=2))


def test_tensordot_double_contraction_ngt2():
    # Regression test for https://github.com/dask/dask/issues/5472
    x = np.arange(60.0).reshape(3, 4, 5)
    y = np.arange(60.0).reshape(4, 5, 3)
    u = da.from_array(x)
    v = da.from_array(y)

    assert_eq(da.tensordot(u, v, axes=2), np.tensordot(x, y, axes=2))

    x = np.arange(60.0).reshape(3, 4, 5)
    y = np.arange(60.0).reshape(4, 5, 3)
    u = da.from_array(x, chunks=3)
    v = da.from_array(y)

    assert_eq(da.tensordot(u, v, axes=2), np.tensordot(x, y, axes=2))


def test_tensordot_more_than_26_dims():
    ndim = 27
    x = np.broadcast_to(1, [2] * ndim)
    dx = da.from_array(x, chunks=-1)
    assert_eq(da.tensordot(dx, dx, ndim), np.array(2**ndim))


def test_dot_method():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 5))
    y = np.arange(200).reshape((20, 10))
    b = da.from_array(y, chunks=(5, 5))

    assert_eq(a.dot(b), x.dot(y))


def test_dot_persist_equivalence():
    # Regression test for https://github.com/dask/dask/issues/6907
    x = da.random.default_rng().random((4, 4), chunks=(2, 2))
    x[x < 0.65] = 0
    y = x.persist()
    z = x.compute()
    r1 = da.dot(x, x).compute()
    r2 = da.dot(y, y).compute()
    rr = np.dot(z, z)
    assert np.allclose(rr, r1)
    assert np.allclose(rr, r2)


@pytest.mark.parametrize("shape, chunks", [((20,), (6,)), ((4, 5), (2, 3))])
def test_vdot(shape, chunks):
    rng = np.random.default_rng(1337)

    x = 2 * rng.random((2,) + shape) - 1
    x = x[0] + 1j * x[1]

    y = 2 * rng.random((2,) + shape) - 1
    y = y[0] + 1j * y[1]

    a = da.from_array(x, chunks=chunks)
    b = da.from_array(y, chunks=chunks)

    assert_eq(np.vdot(x, y), da.vdot(a, b))
    assert_eq(np.vdot(y, x), da.vdot(b, a))
    assert_eq(da.vdot(a, b), da.vdot(b, a).conj())


@pytest.mark.parametrize("shape1, shape2", [((20,), (6,)), ((4, 5), (2, 3))])
def test_outer(shape1, shape2):
    rng = np.random.default_rng(1337)

    x = 2 * rng.random(shape1) - 1
    y = 2 * rng.random(shape2) - 1

    a = da.from_array(x, chunks=3)
    b = da.from_array(y, chunks=3)

    assert_eq(np.outer(x, y), da.outer(a, b))
    assert_eq(np.outer(y, x), da.outer(b, a))


@pytest.mark.parametrize(
    "func1d_name, func1d, specify_output_props",
    [
        ["ndim", lambda x: x.ndim, False],
        ["sum", lambda x: x.sum(), False],
        ["range", lambda x: [x.min(), x.max()], False],
        ["range2", lambda x: [[x.min(), x.max()], [x.max(), x.min()]], False],
        ["cumsum", lambda x: np.cumsum(x), True],
    ],
)
@pytest.mark.parametrize(
    "input_shape, axis",
    [[(10, 15, 20), 0], [(10, 15, 20), 1], [(10, 15, 20), 2], [(10, 15, 20), -1]],
)
def test_apply_along_axis(func1d_name, func1d, specify_output_props, input_shape, axis):
    a = np.random.default_rng().integers(0, 10, input_shape)
    d = da.from_array(a, chunks=(len(input_shape) * (5,)))

    output_shape = None
    output_dtype = None

    if specify_output_props:
        slices = [0] * a.ndim
        slices[axis] = slice(None)
        slices = tuple(slices)
        sample = np.array(func1d(a[slices]))
        output_shape = sample.shape
        output_dtype = sample.dtype

    assert_eq(
        da.apply_along_axis(func1d, axis, d, dtype=output_dtype, shape=output_shape),
        np.apply_along_axis(func1d, axis, a),
    )


@pytest.mark.parametrize(
    "func_name, func",
    [
        ["sum0", lambda x, axis: x.sum(axis=axis)],
        ["sum1", lambda x, axis: x.sum(axis=axis, keepdims=True)],
        [
            "range",
            lambda x, axis: np.concatenate(
                [x.min(axis=axis, keepdims=True), x.max(axis=axis, keepdims=True)],
                axis=axis,
            ),
        ],
    ],
)
@pytest.mark.parametrize(
    "shape, axes",
    [
        [(10, 15, 20), tuple()],
        [(10, 15, 20), 0],
        [(10, 15, 20), (1,)],
        [(10, 15, 20), (-1, 1)],
        [(10, 15, 20), (2, 0, 1)],
    ],
)
def test_apply_over_axes(func_name, func, shape, axes):
    a = np.random.default_rng().integers(0, 10, shape)
    d = da.from_array(a, chunks=(len(shape) * (5,)))

    assert_eq(da.apply_over_axes(func, d, axes), np.apply_over_axes(func, a, axes))


@pytest.mark.parametrize(
    "shape, axis",
    [
        [(10, 15, 20), None],
        [(10, 15, 20), 0],
        [(10, 15, 20), 1],
        [(10, 15, 20), 2],
        [(10, 15, 20), -1],
    ],
)
def test_ptp(shape, axis):
    a = np.random.default_rng().integers(0, 10, shape)
    d = da.from_array(a, chunks=(len(shape) * (5,)))

    assert_eq(da.ptp(d, axis), np.ptp(a, axis))


@pytest.mark.parametrize(
    "shape, axis",
    [[(10, 15, 20), 0], [(10, 15, 20), 1], [(10, 15, 20), 2], [(10, 15, 20), -1]],
)
@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff(shape, n, axis):
    x = np.random.default_rng().integers(0, 10, shape)
    a = da.from_array(x, chunks=(len(shape) * (5,)))

    assert_eq(da.diff(a, n, axis), np.diff(x, n, axis))


@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff_prepend(n):
    x = np.arange(5) + 1
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, prepend=0), np.diff(x, n, prepend=0))
    assert_eq(da.diff(a, n, prepend=[0]), np.diff(x, n, prepend=[0]))
    assert_eq(da.diff(a, n, prepend=[-1, 0]), np.diff(x, n, prepend=[-1, 0]))

    x = np.arange(16).reshape(4, 4)
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, axis=1, prepend=0), np.diff(x, n, axis=1, prepend=0))
    assert_eq(
        da.diff(a, n, axis=1, prepend=[[0], [0], [0], [0]]),
        np.diff(x, n, axis=1, prepend=[[0], [0], [0], [0]]),
    )
    assert_eq(da.diff(a, n, axis=0, prepend=0), np.diff(x, n, axis=0, prepend=0))
    assert_eq(
        da.diff(a, n, axis=0, prepend=[[0, 0, 0, 0]]),
        np.diff(x, n, axis=0, prepend=[[0, 0, 0, 0]]),
    )

    if n > 0:
        # When order is 0 the result is the input array, it doesn't raise
        # an error
        with pytest.raises(ValueError):
            da.diff(a, n, prepend=np.zeros((3, 3)))


@pytest.mark.parametrize("n", [0, 1, 2])
def test_diff_append(n):
    x = np.arange(5) + 1
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, append=0), np.diff(x, n, append=0))
    assert_eq(da.diff(a, n, append=[0]), np.diff(x, n, append=[0]))
    assert_eq(da.diff(a, n, append=[-1, 0]), np.diff(x, n, append=[-1, 0]))

    x = np.arange(16).reshape(4, 4)
    a = da.from_array(x, chunks=2)
    assert_eq(da.diff(a, n, axis=1, append=0), np.diff(x, n, axis=1, append=0))
    assert_eq(
        da.diff(a, n, axis=1, append=[[0], [0], [0], [0]]),
        np.diff(x, n, axis=1, append=[[0], [0], [0], [0]]),
    )
    assert_eq(da.diff(a, n, axis=0, append=0), np.diff(x, n, axis=0, append=0))
    assert_eq(
        da.diff(a, n, axis=0, append=[[0, 0, 0, 0]]),
        np.diff(x, n, axis=0, append=[[0, 0, 0, 0]]),
    )

    if n > 0:
        with pytest.raises(ValueError):
            # When order is 0 the result is the input array, it doesn't raise
            # an error
            da.diff(a, n, append=np.zeros((3, 3)))


def test_diff_negative_order():
    with pytest.raises(ValueError):
        da.diff(da.arange(10), -1)


@pytest.mark.parametrize("shape", [(10,), (10, 15)])
@pytest.mark.parametrize("to_end, to_begin", [[None, None], [0, 0], [[1, 2], [3, 4]]])
def test_ediff1d(shape, to_end, to_begin):
    x = np.random.default_rng().integers(0, 10, shape)
    a = da.from_array(x, chunks=(len(shape) * (5,)))

    assert_eq(da.ediff1d(a, to_end, to_begin), np.ediff1d(x, to_end, to_begin))


@pytest.mark.parametrize(
    "shape, varargs, axis",
    [
        [(10, 15, 20), (), None],
        [(10, 15, 20), (2,), None],
        [(10, 15, 20), (1.0, 1.5, 2.0), None],
        [(10, 15, 20), (), 0],
        [(10, 15, 20), (), 1],
        [(10, 15, 20), (), 2],
        [(10, 15, 20), (), -1],
        [(10, 15, 20), (), (0, 2)],
        [(10, 15, 20), (np.exp(np.arange(10)), np.exp(np.arange(20))), (0, 2)],
        [(10, 15, 20), (0.5, np.exp(np.arange(20))), (0, 2)],
        [(10, 15, 20), (np.exp(np.arange(20)),), -1],
    ],
)
@pytest.mark.parametrize("edge_order", [1, 2])
def test_gradient(shape, varargs, axis, edge_order):
    a = np.random.default_rng().integers(0, 10, shape)
    d_a = da.from_array(a, chunks=(len(shape) * (5,)))

    r_a = np.gradient(a, *varargs, axis=axis, edge_order=edge_order)
    r_d_a = da.gradient(d_a, *varargs, axis=axis, edge_order=edge_order)

    if isinstance(axis, Number):
        assert_eq(r_d_a, r_a)
    else:
        assert len(r_d_a) == len(r_a)

        for e_r_d_a, e_r_a in zip(r_d_a, r_a):
            assert_eq(e_r_d_a, e_r_a)

        assert_eq(
            da.sqrt(sum(map(da.square, r_d_a))), np.sqrt(sum(map(np.square, r_a)))
        )


def test_bincount():
    x = np.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2)
    e = da.bincount(d, minlength=6)
    assert_eq(e, np.bincount(x, minlength=6))
    assert same_keys(da.bincount(d, minlength=6), e)
    assert e.shape == (6,)  # shape equal to minlength
    assert e.chunks == ((6,),)

    assert da.bincount(d, minlength=6).name != da.bincount(d, minlength=7).name
    assert da.bincount(d, minlength=6).name == da.bincount(d, minlength=6).name

    expected_output = np.array([0, 2, 2, 0, 0, 1], dtype=e.dtype)
    assert_eq(e[0:], expected_output)  # can bincount result be sliced


@pytest.mark.parametrize(
    "weights",
    [
        np.array([1, 2, 1, 0.5, 1], dtype=np.float32),
        np.array([1, 2, 1, 0, 1], dtype=np.int32),
    ],
)
def test_bincount_with_weights(weights):
    x = np.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2)

    dweights = da.from_array(weights, chunks=2)
    e = da.bincount(d, weights=dweights, minlength=6)
    assert_eq(e, np.bincount(x, weights=dweights.compute(), minlength=6))
    assert same_keys(da.bincount(d, weights=dweights, minlength=6), e)


def test_bincount_unspecified_minlength():
    x = np.array([1, 1, 3, 7, 0])
    d = da.from_array(x, chunks=2)
    e = da.bincount(d)
    assert_eq(e, np.bincount(x))
    assert same_keys(da.bincount(d), e)
    assert len(e.compute()) == 8  # shape is (nan,) so must compute for len()


def test_digitize():
    x = np.array([2, 4, 5, 6, 1])
    bins = np.array([1, 2, 3, 4, 5])
    for chunks in [2, 4]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            assert_eq(
                da.digitize(d, bins, right=right), np.digitize(x, bins, right=right)
            )

    x = np.random.default_rng().random(size=(100, 100))
    bins = np.random.default_rng().random(size=13)
    bins.sort()
    for chunks in [(10, 10), (10, 20), (13, 17), (87, 54)]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            assert_eq(
                da.digitize(d, bins, right=right), np.digitize(x, bins, right=right)
            )


@pytest.mark.parametrize(
    "a, a_chunks, v, v_chunks",
    [
        [[], 1, [], 1],
        [[0], 1, [0], 1],
        [[-10, 0, 10, 20, 30], 3, [11, 30], 2],
        [[-10, 0, 10, 20, 30], 3, [11, 30, -20, 1, -10, 10, 37, 11], 5],
        [[-10, 0, 10, 20, 30], 3, [[11, 30, -20, 1, -10, 10, 37, 11]], 5],
        [[-10, 0, 10, 20, 30], 3, [[7, 0], [-10, 10], [11, -1], [15, 15]], (2, 2)],
    ],
)
@pytest.mark.parametrize("side", ["left", "right"])
def test_searchsorted(a, a_chunks, v, v_chunks, side):
    a = np.array(a)
    v = np.array(v)

    ad = da.asarray(a, chunks=a_chunks)
    vd = da.asarray(v, chunks=v_chunks)

    out = da.searchsorted(ad, vd, side)

    assert out.shape == vd.shape
    assert out.chunks == vd.chunks
    assert_eq(out, np.searchsorted(a, v, side))


def test_searchsorted_sorter_not_implemented():
    with pytest.raises(NotImplementedError):
        da.searchsorted(da.asarray([1, 0]), da.asarray([1]), sorter=da.asarray([1, 0]))


def test_histogram():
    # Test for normal, flattened input
    n = 100
    v = da.random.default_rng().random(n, chunks=10)
    bins = np.arange(0, 1.01, 0.01)
    (a1, b1) = da.histogram(v, bins=bins)
    (a2, b2) = np.histogram(v, bins=bins)

    # Check if the sum of the bins equals the number of samples
    assert a2.sum(axis=0) == n
    assert a1.sum(axis=0) == n
    assert_eq(a1, a2)
    assert same_keys(da.histogram(v, bins=bins)[0], a1)


def test_histogram_alternative_bins_range():
    v = da.random.default_rng().random(100, chunks=10)
    (a1, b1) = da.histogram(v, bins=10, range=(0, 1))
    (a2, b2) = np.histogram(v, bins=10, range=(0, 1))
    assert_eq(a1, a2)
    assert_eq(b1, b2)


def test_histogram_bins_range_with_nan_array():
    # Regression test for issue #3977
    v = da.from_array(np.array([-2, np.nan, 2]), chunks=1)
    (a1, b1) = da.histogram(v, bins=10, range=(-3, 3))
    (a2, b2) = np.histogram(v, bins=10, range=(-3, 3))
    assert_eq(a1, a2)
    assert_eq(b1, b2)


def test_histogram_return_type():
    v = da.random.default_rng().random(100, chunks=10)
    bins = np.arange(0, 1.01, 0.01)
    # Check if return type is same as hist
    bins = np.arange(0, 11, 1, dtype="i4")
    assert_eq(da.histogram(v * 10, bins=bins)[0], np.histogram(v * 10, bins=bins)[0])


def test_histogram_extra_args_and_shapes():
    # Check for extra args and shapes
    bins = np.arange(0, 1.01, 0.01)
    v = da.random.default_rng().random(100, chunks=10)
    data = [
        (v, bins, da.ones(100, chunks=v.chunks) * 5),
        (
            da.random.default_rng().random((50, 50), chunks=10),
            bins,
            da.ones((50, 50), chunks=10) * 5,
        ),
    ]

    for v, bins, w in data:
        # density
        assert_eq(
            da.histogram(v, bins=bins, density=True)[0],
            np.histogram(v, bins=bins, density=True)[0],
        )

        # weights
        assert_eq(
            da.histogram(v, bins=bins, weights=w)[0],
            np.histogram(v, bins=bins, weights=w)[0],
        )

        assert_eq(
            da.histogram(v, bins=bins, weights=w, density=True)[0],
            da.histogram(v, bins=bins, weights=w, density=True)[0],
        )


def test_histogram_normed_deprecation():
    x = da.arange(10)
    with pytest.raises(ValueError) as info:
        da.histogram(x, bins=[1, 2, 3], normed=True)

    assert "density" in str(info.value)
    assert "deprecated" in str(info.value).lower()


@pytest.mark.parametrize(
    "bins, hist_range",
    [
        (None, None),
        (10, None),
        (10, 1),
        (None, (1, 10)),
        (10, [0, 1, 2]),
        (10, [0]),
        (10, np.array([[0, 1]])),
        (10, da.array([[0, 1]])),
        ([[0, 1, 2]], None),
        (np.array([[0, 1, 2]]), None),
        (da.array([[0, 1, 2]]), None),
    ],
)
def test_histogram_bin_range_raises(bins, hist_range):
    data = da.random.default_rng().random(10, chunks=2)
    with pytest.raises((ValueError, TypeError)) as info:
        da.histogram(data, bins=bins, range=hist_range)
    err_msg = str(info.value)
    assert "bins" in err_msg or "range" in err_msg


@pytest.mark.parametrize("density", [True, False])
@pytest.mark.parametrize("weighted", [True, False])
@pytest.mark.parametrize("non_delayed_i", [None, 0, 1])
@pytest.mark.parametrize("delay_n_bins", [False, True])
def test_histogram_delayed_range(density, weighted, non_delayed_i, delay_n_bins):
    n = 100
    v = np.random.default_rng().random(n)
    vd = da.from_array(v, chunks=10)

    if weighted:
        weights = np.random.default_rng().random(n)
        weights_d = da.from_array(weights, chunks=vd.chunks)

    d_range = [vd.min(), vd.max()]
    if non_delayed_i is not None:
        d_range[non_delayed_i] = d_range[non_delayed_i].compute()
    hist_d, bins_d = da.histogram(
        vd,
        bins=da.array(n) if delay_n_bins and not density else n,
        range=d_range,
        density=density,
        weights=weights_d if weighted else None,
    )

    hist, bins = np.histogram(
        v,
        bins=n,
        range=[v.min(), v.max()],
        density=density,
        weights=weights if weighted else None,
    )

    assert_eq(hist_d, hist)
    assert_eq(bins_d, bins)


@pytest.mark.parametrize("density", [True, False])
@pytest.mark.parametrize("weighted", [True, False])
def test_histogram_delayed_bins(density, weighted):
    n = 100
    v = np.random.default_rng().random(n)
    bins = np.array([0, 0.2, 0.5, 0.8, 1])

    vd = da.from_array(v, chunks=10)
    bins_d = da.from_array(bins, chunks=2)

    if weighted:
        weights = np.random.default_rng().random(n)
        weights_d = da.from_array(weights, chunks=vd.chunks)

    hist_d, bins_d2 = da.histogram(
        vd,
        bins=bins_d,
        range=[bins_d[0], bins_d[-1]],
        density=density,
        weights=weights_d if weighted else None,
    )

    hist, bins = np.histogram(
        v,
        bins=bins,
        range=[bins[0], bins[-1]],
        density=density,
        weights=weights if weighted else None,
    )

    assert bins_d is bins_d2
    assert_eq(hist_d, hist)
    assert_eq(bins_d2, bins)


def test_histogram_delayed_n_bins_raises_with_density():
    data = da.random.default_rng().random(10, chunks=2)
    with pytest.raises(
        NotImplementedError, match="`bins` cannot be a scalar Dask object"
    ):
        da.histogram(data, bins=da.array(10), range=[0, 1], density=True)


@pytest.mark.parametrize("weights", [True, False])
@pytest.mark.parametrize("density", [True, False])
@pytest.mark.parametrize("bins", [(5, 6), 5])
def test_histogram2d(weights, density, bins):
    rng = da.random.default_rng()
    n = 800
    b = bins
    r = ((0, 1), (0, 1))
    x = rng.uniform(0, 1, size=(n,), chunks=(200,))
    y = rng.uniform(0, 1, size=(n,), chunks=(200,))
    w = rng.uniform(0.2, 1.1, size=(n,), chunks=(200,)) if weights else None
    a1, b1x, b1y = da.histogram2d(x, y, bins=b, range=r, density=density, weights=w)
    a2, b2x, b2y = np.histogram2d(x, y, bins=b, range=r, density=density, weights=w)
    a3, b3x, b3y = np.histogram2d(
        x.compute(),
        y.compute(),
        bins=b,
        range=r,
        density=density,
        weights=w.compute() if weights else None,
    )
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    if not (weights or density):
        assert a1.sum() == n
        assert a2.sum() == n
    assert same_keys(
        da.histogram2d(x, y, bins=b, range=r, density=density, weights=w)[0],
        a1,
    )
    assert a1.compute().shape == a3.shape


@pytest.mark.parametrize("weights", [True, False])
@pytest.mark.parametrize("density", [True, False])
def test_histogram2d_array_bins(weights, density):
    rng = da.random.default_rng()
    n = 800
    xbins = [0.0, 0.2, 0.6, 0.9, 1.0]
    ybins = [0.0, 0.1, 0.4, 0.5, 1.0]
    b = [xbins, ybins]
    x = rng.uniform(0, 1, size=(n,), chunks=(200,))
    y = rng.uniform(0, 1, size=(n,), chunks=(200,))
    w = rng.uniform(0.2, 1.1, size=(n,), chunks=(200,)) if weights else None
    a1, b1x, b1y = da.histogram2d(x, y, bins=b, density=density, weights=w)
    a2, b2x, b2y = np.histogram2d(x, y, bins=b, density=density, weights=w)
    a3, b3x, b3y = np.histogram2d(
        x.compute(),
        y.compute(),
        bins=b,
        density=density,
        weights=w.compute() if weights else None,
    )
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    if not (weights or density):
        assert a1.sum() == n
        assert a2.sum() == n
    assert same_keys(
        da.histogram2d(x, y, bins=b, density=density, weights=w)[0],
        a1,
    )
    assert a1.compute().shape == a3.shape


def test_histogramdd():
    n1, n2 = 800, 3
    x = da.random.default_rng().uniform(0, 1, size=(n1, n2), chunks=(200, 3))
    bins = [[0, 0.5, 1], [0, 0.25, 0.85, 1], [0, 0.5, 0.8, 1]]
    (a1, b1) = da.histogramdd(x, bins=bins)
    (a2, b2) = np.histogramdd(x, bins=bins)
    (a3, b3) = np.histogramdd(x.compute(), bins=bins)
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    assert a1.sum() == n1
    assert a2.sum() == n1
    assert same_keys(da.histogramdd(x, bins=bins)[0], a1)
    assert a1.compute().shape == a3.shape


def test_histogramdd_seq_of_arrays():
    rng = da.random.default_rng()
    n1 = 800
    x = rng.uniform(size=(n1,), chunks=200)
    y = rng.uniform(size=(n1,), chunks=200)
    bx = [0.0, 0.25, 0.75, 1.0]
    by = [0.0, 0.30, 0.70, 0.8, 1.0]
    (a1, b1) = da.histogramdd([x, y], bins=[bx, by])
    (a2, b2) = np.histogramdd([x, y], bins=[bx, by])
    (a3, b3) = np.histogramdd((x.compute(), y.compute()), bins=[bx, by])
    assert_eq(a1, a2)
    assert_eq(a1, a3)


def test_histogramdd_alternative_bins_range():
    # test for normal input
    n1, n2 = 600, 3
    x = da.random.default_rng().uniform(
        0, 1, size=(n1, n2), chunks=((200, 200, 200), (3,))
    )
    bins = (3, 5, 4)
    ranges = ((0, 1),) * len(bins)
    (a1, b1) = da.histogramdd(x, bins=bins, range=ranges)
    (a2, b2) = np.histogramdd(x, bins=bins, range=ranges)
    (a3, b3) = np.histogramdd(x.compute(), bins=bins, range=ranges)
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    bins = 4
    (a1, b1) = da.histogramdd(x, bins=bins, range=ranges)
    (a2, b2) = np.histogramdd(x, bins=bins, range=ranges)
    assert_eq(a1, a2)

    assert a1.sum() == n1
    assert a2.sum() == n1
    assert same_keys(da.histogramdd(x, bins=bins, range=ranges)[0], a1)


def test_histogramdd_weighted():
    rng = da.random.default_rng()
    # test for normal input
    n1, n2 = 600, 3
    x = rng.uniform(0, 1, size=(n1, n2), chunks=((200, 200, 200), (3,)))
    w = rng.uniform(0.5, 0.8, size=(n1,), chunks=200)
    bins = (3, 5, 4)
    ranges = ((0, 1),) * len(bins)
    (a1, b1) = da.histogramdd(x, bins=bins, range=ranges, weights=w)
    (a2, b2) = np.histogramdd(x, bins=bins, range=ranges, weights=w)
    (a3, b3) = np.histogramdd(x.compute(), bins=bins, range=ranges, weights=w.compute())
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    bins = 4
    (a1, b1) = da.histogramdd(x, bins=bins, range=ranges, weights=w)
    (a2, b2) = np.histogramdd(x, bins=bins, range=ranges, weights=w)
    (a3, b3) = np.histogramdd(x.compute(), bins=bins, range=ranges, weights=w.compute())
    assert_eq(a1, a2)
    assert_eq(a1, a3)


def test_histogramdd_density():
    n1, n2 = 800, 3
    x = da.random.default_rng().uniform(0, 1, size=(n1, n2), chunks=(200, 3))
    bins = [[0, 0.5, 1], [0, 0.25, 0.85, 1], [0, 0.5, 0.8, 1]]
    (a1, b1) = da.histogramdd(x, bins=bins, density=True)
    (a2, b2) = np.histogramdd(x, bins=bins, density=True)
    (a3, b3) = da.histogramdd(x, bins=bins, normed=True)
    (a4, b4) = np.histogramdd(x.compute(), bins=bins, density=True)
    assert_eq(a1, a2)
    assert_eq(a1, a3)
    assert_eq(a1, a4)
    assert same_keys(da.histogramdd(x, bins=bins, density=True)[0], a1)


def test_histogramdd_weighted_density():
    rng = da.random.default_rng()
    n1, n2 = 1200, 4
    x = rng.standard_normal(size=(n1, n2), chunks=(200, 4))
    w = rng.uniform(0.5, 1.2, size=(n1,), chunks=200)
    bins = (5, 6, 7, 8)
    ranges = ((-4, 4),) * len(bins)
    (a1, b1) = da.histogramdd(x, bins=bins, range=ranges, weights=w, density=True)
    (a2, b2) = np.histogramdd(x, bins=bins, range=ranges, weights=w, density=True)
    (a3, b3) = da.histogramdd(x, bins=bins, range=ranges, weights=w, normed=True)
    assert_eq(a1, a2)
    assert_eq(a1, a3)


def test_histogramdd_raises_incompat_sample_chunks():
    data = da.random.default_rng().random(size=(10, 3), chunks=(5, 1))
    with pytest.raises(
        ValueError, match="Input array can only be chunked along the 0th axis"
    ):
        da.histogramdd(data, bins=10, range=((0, 1),) * 3)


def test_histogramdd_raises_incompat_multiarg_chunks():
    rng = da.random.default_rng()
    x = rng.random(size=(10,), chunks=2)
    y = rng.random(size=(10,), chunks=2)
    z = rng.random(size=(10,), chunks=5)
    with pytest.raises(
        ValueError, match="All coordinate arrays must be chunked identically."
    ):
        da.histogramdd((x, y, z), bins=(3,) * 3, range=((0, 1),) * 3)


def test_histogramdd_raises_incompat_weight_chunks():
    rng = da.random.default_rng()
    x = rng.random(size=(10,), chunks=2)
    y = rng.random(size=(10,), chunks=2)
    z = da.atleast_2d((x, y)).T.rechunk((2, 2))
    w = rng.random(size=(10,), chunks=5)
    with pytest.raises(
        ValueError,
        match="Input arrays and weights must have the same shape and chunk structure.",
    ):
        da.histogramdd((x, y), bins=(3,) * 2, range=((0, 1),) * 2, weights=w)
    with pytest.raises(
        ValueError,
        match="Input array and weights must have the same shape and chunk structure along the first dimension.",
    ):
        da.histogramdd(z, bins=(3,) * 2, range=((0, 1),) * 2, weights=w)


def test_histogramdd_raises_incompat_bins_or_range():
    data = da.random.default_rng().random(size=(10, 4), chunks=(5, 4))
    bins = (2, 3, 4, 5)
    ranges = ((0, 1),) * len(bins)

    # bad number of bins defined (should be data.shape[1])
    bins = (2, 3, 4)
    with pytest.raises(
        ValueError,
        match="The dimension of bins must be equal to the dimension of the sample.",
    ):
        da.histogramdd(data, bins=bins, range=ranges)

    # one range per dimension is required.
    bins = (2, 3, 4, 5)
    ranges = ((0, 1),) * 3
    with pytest.raises(
        ValueError,
        match="range argument requires one entry, a min max pair, per dimension.",
    ):
        da.histogramdd(data, bins=bins, range=ranges)

    # has range elements that are not pairs
    with pytest.raises(
        ValueError, match="range argument should be a sequence of pairs"
    ):
        da.histogramdd(data, bins=bins, range=((0, 1), (0, 1, 2), 3, 5))


def test_histogramdd_raise_normed_and_density():
    data = da.random.default_rng().random(size=(10, 3), chunks=(5, 3))
    bins = (4, 5, 6)
    ranges = ((0, 1),) * 3
    with pytest.raises(TypeError, match="Cannot specify both 'normed' and 'density'"):
        da.histogramdd(data, bins=bins, range=ranges, normed=True, density=True)


def test_histogramdd_raise_incompat_shape():
    # 1D
    data = da.random.default_rng().random(size=(10,), chunks=(2,))
    with pytest.raises(
        ValueError, match="Single array input to histogramdd should be columnar"
    ):
        da.histogramdd(data, bins=4, range=((-3, 3),))
    # 3D (not columnar)
    data = da.random.default_rng().random(size=(4, 4, 4), chunks=(2, 2, 2))
    with pytest.raises(
        ValueError, match="Single array input to histogramdd should be columnar"
    ):
        da.histogramdd(data, bins=4, range=((-3, 3),))


def test_histogramdd_edges():
    data = da.random.default_rng().random(size=(10, 3), chunks=(5, 3))
    edges = [
        np.array([0.1, 0.3, 0.8, 1.0]),
        np.array([0.2, 0.3, 0.8, 0.9]),
        np.array([0.1, 0.5, 0.7]),
    ]
    # passing bins as an array of bin edges.
    a1, b1 = da.histogramdd(data, bins=edges)
    a2, b2 = np.histogramdd(data.compute(), bins=edges)
    for ib1, ib2 in zip(b1, b2):
        assert_eq(ib1, ib2)
    # passing bins as an int with range definitions
    a1, b1 = da.histogramdd(data, bins=5, range=((0, 1),) * 3)
    a2, b2 = np.histogramdd(data.compute(), bins=5, range=((0, 1),) * 3)
    for ib1, ib2 in zip(b1, b2):
        assert_eq(ib1, ib2)


def test_cov():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(4, 4))

    assert_eq(da.cov(d), np.cov(x))
    assert_eq(da.cov(d, rowvar=0), np.cov(x, rowvar=0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # dof <= 0 for slice
        assert_eq(da.cov(d, ddof=10), np.cov(x, ddof=10))
    assert_eq(da.cov(d, bias=1), np.cov(x, bias=1))
    assert_eq(da.cov(d, d), np.cov(x, x))

    y = np.arange(8)
    e = da.from_array(y, chunks=(4,))

    assert_eq(da.cov(d, e), np.cov(x, y))
    assert_eq(da.cov(e, d), np.cov(y, x))

    with pytest.raises(ValueError):
        da.cov(d, ddof=1.5)


def test_corrcoef():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(4, 4))

    assert_eq(da.corrcoef(d), np.corrcoef(x))
    assert_eq(da.corrcoef(d, rowvar=0), np.corrcoef(x, rowvar=0))
    assert_eq(da.corrcoef(d, d), np.corrcoef(x, x))

    y = np.arange(8)
    e = da.from_array(y, chunks=(4,))

    assert_eq(da.corrcoef(d, e), np.corrcoef(x, y))
    assert_eq(da.corrcoef(e, d), np.corrcoef(y, x))


def test_round():
    x = np.random.default_rng().random(10)
    d = da.from_array(x, chunks=4)

    for i in (0, 1, 4, 5):
        assert_eq(x.round(i), d.round(i))

    assert_eq(d.round(2), da.round(d, 2))


@pytest.mark.parametrize("return_index", [False, True])
@pytest.mark.parametrize("return_inverse", [False, True])
@pytest.mark.parametrize("return_counts", [False, True])
def test_unique_kwargs(return_index, return_inverse, return_counts):
    kwargs = dict(
        return_index=return_index,
        return_inverse=return_inverse,
        return_counts=return_counts,
    )

    a = np.array([1, 2, 4, 4, 5, 2])
    d = da.from_array(a, chunks=(3,))

    r_a = np.unique(a, **kwargs)
    r_d = da.unique(d, **kwargs)

    if not any([return_index, return_inverse, return_counts]):
        assert isinstance(r_a, np.ndarray)
        assert isinstance(r_d, da.Array)

        r_a = (r_a,)
        r_d = (r_d,)

    assert len(r_a) == len(r_d)

    if return_inverse:
        i = 1 + int(return_index)
        assert (d.size,) == r_d[i].shape

    for e_r_a, e_r_d in zip(r_a, r_d):
        assert_eq(e_r_d, e_r_a)


@pytest.mark.parametrize("seed", [23, 796])
@pytest.mark.parametrize(
    "shape, chunks",
    [[(10,), (5,)], [(10,), (3,)], [(4, 5), (3, 2)], [(20, 20), (4, 5)]],
)
def test_unique_rand(seed, shape, chunks):
    rng = np.random.default_rng(seed)

    a = rng.integers(0, 10, size=shape)
    d = da.from_array(a, chunks=chunks)

    r_a = np.unique(a, return_index=True, return_inverse=True, return_counts=True)
    r_d = da.unique(d, return_index=True, return_inverse=True, return_counts=True)

    assert_eq(r_d[0], r_a[0])
    assert_eq(r_d[1], r_a[1])
    assert_eq(r_d[2], r_a[2])
    assert_eq(r_d[3], r_a[3])


@pytest.mark.parametrize("seed", [23, 796])
@pytest.mark.parametrize("low, high", [[0, 10]])
@pytest.mark.parametrize(
    "elements_shape, elements_chunks",
    [[(10,), (5,)], [(10,), (3,)], [(4, 5), (3, 2)], [(20, 20), (4, 5)]],
)
@pytest.mark.parametrize(
    "test_shape, test_chunks",
    [[(10,), (5,)], [(10,), (3,)], [(4, 5), (3, 2)], [(20, 20), (4, 5)]],
)
@pytest.mark.parametrize("invert", [True, False])
def test_isin_rand(
    seed, low, high, elements_shape, elements_chunks, test_shape, test_chunks, invert
):
    rng = np.random.default_rng(seed)

    a1 = rng.integers(low, high, size=elements_shape)
    d1 = da.from_array(a1, chunks=elements_chunks)

    a2 = rng.integers(low, high, size=test_shape) - 5
    d2 = da.from_array(a2, chunks=test_chunks)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=da.PerformanceWarning)
        r_a = np.isin(a1, a2, invert=invert)
        r_d = da.isin(d1, d2, invert=invert)
    assert_eq(r_a, r_d)


@pytest.mark.parametrize("assume_unique", [True, False])
def test_isin_assume_unique(assume_unique):
    a1 = np.arange(10)
    d1 = da.from_array(a1, chunks=(5,))

    test_elements = np.arange(0, 10, 2)
    r_a = np.isin(a1, test_elements, assume_unique=assume_unique)
    r_d = da.isin(d1, test_elements, assume_unique=assume_unique)
    assert_eq(r_a, r_d)


def _maybe_len(l):
    try:
        return len(l)
    except TypeError:
        return 0


@pytest.mark.parametrize("chunks", [(4, 6), (2, 6)])
@pytest.mark.parametrize("shift", [3, 7, 9, (3, 9), (7, 2)])
@pytest.mark.parametrize("axis", [None, 0, 1, -1, (0, 1), (1, 0)])
def test_roll(chunks, shift, axis):
    x = np.random.default_rng().integers(10, size=(4, 6))
    a = da.from_array(x, chunks=chunks)

    if _maybe_len(shift) != _maybe_len(axis):
        with pytest.raises(TypeError if axis is None else ValueError):
            da.roll(a, shift, axis)
    else:
        assert_eq(np.roll(x, shift, axis), da.roll(a, shift, axis))


def test_roll_always_results_in_a_new_array():
    x = da.arange(2, 3)
    y = da.roll(x, 1)
    assert y is not x


def test_roll_works_even_if_shape_is_0():
    expected = np.roll(np.zeros(0), 0)
    actual = da.roll(da.zeros(0), 0)
    assert_eq(expected, actual)


@pytest.mark.parametrize("shape", [(10,), (5, 10), (5, 10, 10)])
def test_shape_and_ndim(shape):
    x = da.random.default_rng().random(shape)
    assert np.shape(x) == shape

    x = da.random.default_rng().random(shape)
    assert np.ndim(x) == len(shape)


@pytest.mark.parametrize(
    "shape", [((12,), (12,)), ((4, 3), (3, 4)), ((12,), (1, 6, 2))]
)
@pytest.mark.parametrize("reverse", [True, False])
def test_union1d(shape, reverse):
    s1, s2 = shape
    x1 = np.arange(12).reshape(s1)
    x2 = np.arange(6, 18).reshape(s2)

    if reverse:
        x1 = x1[::-1]

    dx1 = da.from_array(x1)
    dx2 = da.from_array(x2)

    result = np.union1d(dx1, dx2)
    expected = np.union1d(x1, x2)

    assert isinstance(result, da.Array)

    assert_eq(result, expected)


def test_ravel():
    x = np.random.default_rng().integers(10, size=(4, 6))

    # 2d
    for chunks in [(4, 6), (2, 6)]:
        a = da.from_array(x, chunks=chunks)
        assert_eq(x.ravel(), a.ravel())
        assert len(a.ravel().dask) == len(a.dask) + len(a.chunks[0])

    # 0d
    assert_eq(x[0, 0].ravel(), a[0, 0].ravel())

    # 1d
    a_flat = a.ravel()
    assert_eq(a_flat.ravel(), a_flat)

    # 3d
    x = np.random.default_rng().integers(10, size=(2, 3, 4))
    for chunks in [4, (1, 3, 4)]:
        a = da.from_array(x, chunks=chunks)
        assert_eq(x.ravel(), a.ravel())

    assert_eq(x.flatten(), a.flatten())
    assert_eq(np.ravel(x), da.ravel(a))


def test_ravel_1D_no_op():
    x = np.random.default_rng().integers(10, size=100)
    dx = da.from_array(x, chunks=10)
    # known dims
    assert_eq(dx.ravel(), x.ravel())
    # Unknown dims
    assert_eq(dx[dx > 2].ravel(), x[x > 2].ravel())


def test_ravel_with_array_like():
    # int
    assert_eq(np.ravel(0), da.ravel(0))
    assert isinstance(da.ravel(0), da.core.Array)

    # list
    assert_eq(np.ravel([0, 0]), da.ravel([0, 0]))
    assert isinstance(da.ravel([0, 0]), da.core.Array)

    # tuple
    assert_eq(np.ravel((0, 0)), da.ravel((0, 0)))
    assert isinstance(da.ravel((0, 0)), da.core.Array)

    # nested i.e. tuples in list
    assert_eq(np.ravel([(0,), (0,)]), da.ravel([(0,), (0,)]))
    assert isinstance(da.ravel([(0,), (0,)]), da.core.Array)


@pytest.mark.parametrize("axis", [None, 0, 1, -1, (0, 1), (0, 2), (1, 2), 2])
def test_expand_dims(axis):
    a = np.arange(10)
    d = da.from_array(a, chunks=(3,))

    if axis is None:
        with pytest.raises(TypeError):
            da.expand_dims(d, axis=axis)
    elif axis == 2:
        with pytest.raises(AxisError):
            da.expand_dims(d, axis=axis)
    else:
        a_e = np.expand_dims(a, axis=axis)
        d_e = da.expand_dims(d, axis=axis)

        assert_eq(d_e, a_e)
        assert same_keys(d_e, da.expand_dims(d, axis=axis))


@pytest.mark.parametrize("is_func", [True, False])
@pytest.mark.parametrize("axis", [None, 0, -1, (0, -1)])
def test_squeeze(is_func, axis):
    a = np.arange(10)[None, :, None, None]
    d = da.from_array(a, chunks=(1, 3, 1, 1))

    if is_func:
        a_s = np.squeeze(a, axis=axis)
        d_s = da.squeeze(d, axis=axis)
    else:
        a_s = a.squeeze(axis=axis)
        d_s = d.squeeze(axis=axis)

    assert_eq(d_s, a_s)
    assert same_keys(d_s, da.squeeze(d, axis=axis))

    if axis is None:
        axis = tuple(range(a.ndim))
    else:
        axis = axis if isinstance(axis, tuple) else (axis,)
        axis = tuple(i % a.ndim for i in axis)
    axis = tuple(i for i, c in enumerate(d.chunks) if i in axis and len(c) == 1)

    exp_d_s_chunks = tuple(c for i, c in enumerate(d.chunks) if i not in axis)
    assert d_s.chunks == exp_d_s_chunks


@pytest.mark.parametrize("shape", [(1,), (1, 1)])
def test_squeeze_1d_array(shape):
    a = np.full(shape=shape, fill_value=2)
    a_s = np.squeeze(a)
    d = da.from_array(a, chunks=(1))
    d_s = da.squeeze(d)
    assert isinstance(d_s, da.Array)
    assert isinstance(d_s.compute(), np.ndarray)
    assert_eq(d_s, a_s)


def test_vstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(np.vstack((x, y)), da.vstack((a, b)))
    assert_eq(np.vstack((x, y[None, :])), da.vstack((a, b[None, :])))


def test_hstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(np.hstack((x[None, :], y[None, :])), da.hstack((a[None, :], b[None, :])))
    assert_eq(np.hstack((x, y)), da.hstack((a, b)))


def test_dstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(
        np.dstack((x[None, None, :], y[None, None, :])),
        da.dstack((a[None, None, :], b[None, None, :])),
    )
    assert_eq(np.dstack((x[None, :], y[None, :])), da.dstack((a[None, :], b[None, :])))
    assert_eq(np.dstack((x, y)), da.dstack((a, b)))


@pytest.mark.parametrize(
    "np_func,dsk_func,nan_chunk",
    [(np.hstack, da.hstack, 0), (np.dstack, da.dstack, 1), (np.vstack, da.vstack, 2)],
)
def test_stack_unknown_chunk_sizes(np_func, dsk_func, nan_chunk):
    shape = (100, 100, 100)
    x = da.ones(shape, chunks=(50, 50, 50))
    y = np.ones(shape)

    tmp = list(x._chunks)
    tmp[nan_chunk] = (np.nan,) * 2
    x._chunks = tuple(tmp)

    with pytest.raises(ValueError):
        dsk_func((x, x))

    np_stacked = np_func((y, y))
    dsk_stacked = dsk_func((x, x), allow_unknown_chunksizes=True)
    assert_eq(np_stacked, dsk_stacked)


def test_take():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 5))

    assert_eq(np.take(x, 3, axis=0), da.take(a, 3, axis=0))
    assert_eq(np.take(x, [3, 4, 5], axis=-1), da.take(a, [3, 4, 5], axis=-1))

    with pytest.raises(ValueError):
        da.take(a, 3, axis=2)

    assert same_keys(da.take(a, [3, 4, 5], axis=-1), da.take(a, [3, 4, 5], axis=-1))


def test_take_dask_from_numpy():
    x = np.arange(5).astype("f8")
    y = da.from_array(np.array([1, 2, 3, 3, 2, 1]), chunks=3)

    z = da.take(x * 2, y)

    assert z.chunks == y.chunks
    assert_eq(z, np.array([2.0, 4.0, 6.0, 6.0, 4.0, 2.0]))


def test_compress():
    x = np.arange(25).reshape((5, 5))
    a = da.from_array(x, chunks=(2, 2))

    c1 = np.array([True, False, True, False, True])
    c2 = np.array([True, False])
    c3 = [True, False]
    dc1 = da.from_array(c1, chunks=3)
    dc2 = da.from_array(c2, chunks=2)

    for c, dc in [(c1, c1), (c2, c2), (c3, c3), (c1, dc1), (c2, dc2), (c3, dc2)]:
        for axis in [None, 0, 1]:
            res = da.compress(dc, a, axis=axis)
            assert_eq(np.compress(c, x, axis=axis), res)
            if isinstance(dc, da.Array):
                # If condition is a dask array then we expect the shape of the
                # compressed array to be nan, because we won't know that until
                # the result is computed.
                axis = axis or 0
                assert np.isnan(res.shape[axis]).all()
                assert np.isnan(res.chunks[axis]).all()
            else:
                # If condition is a not a dask array then we expect the shape of the
                # compressed axis to be known, i.e., not nan.
                axis = axis or 0
                assert np.count_nonzero(dc) == res.shape[axis]
                assert not np.isnan(res.chunks[axis]).any()

    with pytest.raises(ValueError):
        da.compress([True, False], a, axis=100)

    with pytest.raises(ValueError):
        da.compress([[True], [False]], a, axis=100)


def test_extract():
    x = np.arange(25).reshape((5, 5))
    a = da.from_array(x, chunks=(2, 2))

    c1 = np.array([True, False, True, False, True])
    c2 = np.array([[True, False], [True, False]])
    c3 = np.array([True, False])
    dc1 = da.from_array(c1, chunks=3)
    dc2 = da.from_array(c2, chunks=(2, 1))
    dc3 = da.from_array(c3, chunks=2)

    for c, dc in [(c1, c1), (c2, c2), (c3, c3), (c1, dc1), (c2, dc2), (c3, dc3)]:
        res = da.extract(dc, a)
        assert_eq(np.extract(c, x), res)
        if isinstance(dc, da.Array):
            assert np.isnan(res.chunks[0]).all()


def test_isnull():
    x = np.array([1, np.nan])
    a = da.from_array(x, chunks=(2,))
    with contextlib.suppress(ImportError):
        assert_eq(da.isnull(a), np.isnan(x))
        assert_eq(da.notnull(a), ~(np.isnan(x)))


def test_isnull_result_is_an_array():
    # regression test for https://github.com/dask/dask/issues/3822
    arr = da.from_array(np.arange(3, dtype=np.int64), chunks=-1)
    with contextlib.suppress(ImportError):
        result = da.isnull(arr[0]).compute()
        assert type(result) is np.ndarray


def test_isclose():
    x = np.array([0, np.nan, 1, 1.5])
    y = np.array([1e-9, np.nan, 1, 2])
    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))
    assert_eq(da.isclose(a, b, equal_nan=True), np.isclose(x, y, equal_nan=True))


def test_allclose():
    n_a = np.array([0, np.nan, 1, 1.5])
    n_b = np.array([1e-9, np.nan, 1, 2])

    d_a = da.from_array(n_a, chunks=(2,))
    d_b = da.from_array(n_b, chunks=(2,))

    n_r = np.allclose(n_a, n_b, equal_nan=True)
    d_r = da.allclose(d_a, d_b, equal_nan=True)

    assert_eq(np.array(n_r)[()], d_r)


def test_choose():
    # test choose function
    x = np.random.default_rng().integers(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(da.choose(d > 5, [0, d]), np.choose(x > 5, [0, x]))
    assert_eq(da.choose(d > 5, [-d, d]), np.choose(x > 5, [-x, x]))

    # test choose method
    index_dask = d > 5
    index_numpy = x > 5
    assert_eq(index_dask.choose([0, d]), index_numpy.choose([0, x]))
    assert_eq(index_dask.choose([-d, d]), index_numpy.choose([-x, x]))


def test_piecewise():
    rng = np.random.default_rng(1337)

    x = rng.integers(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(
        np.piecewise(x, [x < 5, x >= 5], [lambda e, v, k: e + 1, 5], 1, k=2),
        da.piecewise(d, [d < 5, d >= 5], [lambda e, v, k: e + 1, 5], 1, k=2),
    )


def test_piecewise_otherwise():
    rng = np.random.default_rng(1337)

    x = rng.integers(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(
        np.piecewise(
            x,
            [x > 5, x <= 2],
            [lambda e, v, k: e + 1, lambda e, v, k: v * e, lambda e, v, k: 0],
            1,
            k=2,
        ),
        da.piecewise(
            d,
            [d > 5, d <= 2],
            [lambda e, v, k: e + 1, lambda e, v, k: v * e, lambda e, v, k: 0],
            1,
            k=2,
        ),
    )


def test_select():
    conditions = [
        np.array([False, False, False, False]),
        np.array([False, True, False, True]),
        np.array([False, False, True, True]),
    ]
    choices = [
        np.array([1, 2, 3, 4]),
        np.array([5, 6, 7, 8]),
        np.array([9, 10, 11, 12]),
    ]
    d_conditions = da.from_array(conditions, chunks=(3, 2))
    d_choices = da.from_array(choices)
    assert_eq(np.select(conditions, choices), da.select(d_conditions, d_choices))


def test_select_multidimension():
    x = np.random.default_rng().random((100, 50, 2))
    y = da.from_array(x, chunks=(50, 50, 1))
    res_x = np.select([x < 0, x > 2, x > 1], [x, x * 2, x * 3], default=1)
    res_y = da.select([y < 0, y > 2, y > 1], [y, y * 2, y * 3], default=1)
    assert isinstance(res_y, da.Array)
    assert_eq(res_y, res_x)


def test_select_return_dtype():
    d = np.array([1, 2, 3, np.nan, 5, 7])
    m = np.isnan(d)
    d_d = da.from_array(d)
    d_m = da.isnan(d_d)
    assert_eq(np.select([m], [d]), da.select([d_m], [d_d]), equal_nan=True)


@pytest.mark.xfail(reason="broadcasting in da.select() not implemented yet")
def test_select_broadcasting():
    conditions = [np.array(True), np.array([False, True, False])]
    choices = [1, np.arange(12).reshape(4, 3)]
    d_conditions = da.from_array(conditions)
    d_choices = da.from_array(choices)
    assert_eq(np.select(conditions, choices), da.select(d_conditions, d_choices))
    # default can broadcast too:
    assert_eq(np.select([True], [0], default=[0]), da.select([True], [0], default=[0]))


def test_argwhere():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = np.argwhere(x)
        d_nz = da.argwhere(d)

        assert_eq(d_nz, x_nz)


def test_argwhere_obj():
    x = np.random.default_rng().integers(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_nz = np.argwhere(x)
    d_nz = da.argwhere(d)

    assert_eq(d_nz, x_nz)


def test_argwhere_str():
    # We may have behavior differences with NumPy for strings
    # with just spaces, depending on the version of NumPy.
    # https://github.com/numpy/numpy/issues/9875
    x = np.array(list("Hello world"))
    d = da.from_array(x, chunks=(4,))

    x_nz = np.argwhere(x)
    d_nz = da.argwhere(d)

    assert_eq(d_nz, x_nz)


def test_where():
    rng = np.random.default_rng()
    x = rng.integers(10, size=(15, 14))
    x[5, 5] = x[4, 4] = 0  # Ensure some false elements
    d = da.from_array(x, chunks=(4, 5))
    y = rng.integers(10, size=15).astype(np.uint8)
    e = da.from_array(y, chunks=(4,))

    for c1, c2 in [
        (d > 5, x > 5),
        (d, x),
        (1, 1),
        (0, 0),
        (5, 5),
        (True, True),
        (np.True_, np.True_),
        (False, False),
        (np.False_, np.False_),
    ]:
        for b1, b2 in [(0, 0), (-e[:, None], -y[:, None]), (e[:14], y[:14])]:
            w1 = da.where(c1, d, b1)
            w2 = np.where(c2, x, b2)
            assert_eq(w1, w2)


def test_where_scalar_dtype():
    x = np.int32(3)
    y1 = np.array([4, 5, 6], dtype=np.int16)
    c1 = np.array([1, 0, 1])
    y2 = da.from_array(y1, chunks=2)
    c2 = da.from_array(c1, chunks=2)
    w1 = np.where(c1, x, y1)
    w2 = da.where(c2, x, y2)
    assert_eq(w1, w2)
    # Test again for the bool optimization
    w3 = np.where(True, x, y1)
    w4 = da.where(True, x, y1)
    assert_eq(w3, w4)


def test_where_bool_optimization():
    rng = np.random.default_rng()
    x = rng.integers(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))
    y = rng.integers(10, size=(15, 16))
    e = da.from_array(y, chunks=(4, 5))

    for c in [True, False, np.True_, np.False_, 1, 0]:
        w1 = da.where(c, d, e)
        w2 = np.where(c, x, y)

        assert_eq(w1, w2)

        ex_w1 = d if c else e

        assert w1 is ex_w1


def test_where_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_w = np.where(x)
        d_w = da.where(d)

        assert isinstance(d_w, type(x_w))
        assert len(d_w) == len(x_w)

        for i in range(len(x_w)):
            assert_eq(d_w[i], x_w[i])


def test_where_incorrect_args():
    a = da.ones(5, chunks=3)

    for kwd in ["x", "y"]:
        kwargs = {kwd: a}
        try:
            da.where(a > 0, **kwargs)
        except ValueError as e:
            assert "either both or neither of x and y should be given" in str(e)


def test_count_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_c = np.count_nonzero(x)
        d_c = da.count_nonzero(d)

        if d_c.shape == tuple():
            assert x_c == d_c.compute()
        else:
            assert_eq(x_c, d_c)


@pytest.mark.parametrize("axis", [None, 0, (1,), (0, 1)])
def test_count_nonzero_axis(axis):
    for shape, chunks in [((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_c = np.count_nonzero(x, axis)
        d_c = da.count_nonzero(d, axis)

        if d_c.shape == tuple():
            assert x_c == d_c.compute()
        else:
            assert_eq(x_c, d_c)


def test_count_nonzero_obj():
    x = np.random.default_rng().integers(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_c = np.count_nonzero(x)
    d_c = da.count_nonzero(d)

    if d_c.shape == tuple():
        assert x_c == d_c.compute()
    else:
        assert_eq(x_c, d_c)


@pytest.mark.parametrize("axis", [None, 0, (1,), (0, 1)])
def test_count_nonzero_obj_axis(axis):
    x = np.random.default_rng().integers(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_c = np.count_nonzero(x, axis)
    d_c = da.count_nonzero(d, axis)

    if d_c.shape == tuple():
        assert x_c == d_c.compute()
    else:
        #######################################################
        # Workaround oddness with Windows and object arrays.  #
        #                                                     #
        # xref: https://github.com/numpy/numpy/issues/9468    #
        #######################################################
        assert_eq(x_c.astype(np.intp), d_c)


def test_count_nonzero_str():
    # We may have behavior differences with NumPy for strings
    # with just spaces, depending on the version of NumPy.
    # https://github.com/numpy/numpy/issues/9875
    x = np.array(list("Hello world"))
    d = da.from_array(x, chunks=(4,))

    x_c = np.count_nonzero(x)
    d_c = da.count_nonzero(d)

    assert x_c == d_c.compute()


def test_flatnonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_fnz = np.flatnonzero(x)
        d_fnz = da.flatnonzero(d)

        assert_eq(d_fnz, x_fnz)


def test_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = np.nonzero(x)
        d_nz = da.nonzero(d)

        assert isinstance(d_nz, type(x_nz))
        assert len(d_nz) == len(x_nz)

        for i in range(len(x_nz)):
            assert_eq(d_nz[i], x_nz[i])


def test_nonzero_method():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.default_rng().integers(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = x.nonzero()
        d_nz = d.nonzero()

        assert isinstance(d_nz, type(x_nz))
        assert len(d_nz) == len(x_nz)

        for i in range(len(x_nz)):
            assert_eq(d_nz[i], x_nz[i])


def test_unravel_index_empty():
    shape = tuple()
    findices = np.array(0, dtype=int)
    d_findices = da.from_array(findices, chunks=1)

    indices = np.unravel_index(findices, shape)
    d_indices = da.unravel_index(d_findices, shape)

    assert isinstance(d_indices, type(indices))
    assert len(d_indices) == len(indices) == 0


def test_unravel_index():
    rng = np.random.default_rng()
    for nindices, shape, order in [
        (0, (15,), "C"),
        (1, (15,), "C"),
        (3, (15,), "C"),
        (3, (15,), "F"),
        (2, (15, 16), "C"),
        (2, (15, 16), "F"),
    ]:
        arr = rng.random(shape)
        darr = da.from_array(arr, chunks=1)

        findices = rng.integers(np.prod(shape, dtype=int), size=nindices)
        d_findices = da.from_array(findices, chunks=1)

        indices = np.unravel_index(findices, shape, order)
        d_indices = da.unravel_index(d_findices, shape, order)

        assert isinstance(d_indices, type(indices))
        assert len(d_indices) == len(indices)

        for i in range(len(indices)):
            assert_eq(d_indices[i], indices[i])

        assert_eq(darr.vindex[dask.compute(*d_indices)], arr[indices])


@pytest.mark.parametrize(
    "asarray",
    [
        lambda x: x,
        lambda x: [np.asarray(a) for a in x],
        lambda x: [da.asarray(a) for a in x],
        np.asarray,
        da.from_array,
    ],
)
@pytest.mark.parametrize(
    "arr, chunks, kwargs",
    [
        # Numpy doctests:
        ([[3, 6, 6], [4, 5, 1]], (2, 3), dict(dims=(7, 6), order="C")),
        ([[3, 6, 6], [4, 5, 1]], (2, 1), dict(dims=(7, 6), order="F")),
        ([[3, 6, 6], [4, 5, 1]], 1, dict(dims=(4, 6), mode="clip")),
        ([[3, 6, 6], [4, 5, 1]], (2, 3), dict(dims=(4, 4), mode=("clip", "wrap"))),
        # Shape tests:
        ([[3, 6, 6]], (1, 1), dict(dims=(7), order="C")),
        ([[3, 6, 6], [4, 5, 1], [8, 6, 2]], (3, 1), dict(dims=(7, 6, 9), order="C")),
        # Multi-dimensional index arrays
        (
            np.arange(6).reshape(3, 2, 1).tolist(),
            (1, 2, 1),
            dict(dims=(7, 6, 9), order="C"),
        ),
        # Broadcasting index arrays
        ([1, [2, 3]], None, dict(dims=(8, 9))),
        ([1, [2, 3], [[1, 2], [3, 4], [5, 6], [7, 8]]], None, dict(dims=(8, 9, 10))),
    ],
)
def test_ravel_multi_index(asarray, arr, chunks, kwargs):
    if any(np.isscalar(x) for x in arr) and asarray in (np.asarray, da.from_array):
        pytest.skip()

    if asarray is da.from_array:
        arr = np.asarray(arr)
        input = da.from_array(arr, chunks=chunks)
    else:
        arr = input = asarray(arr)

    assert_eq(
        np.ravel_multi_index(arr, **kwargs),
        da.ravel_multi_index(input, **kwargs),
    )


def test_ravel_multi_index_unknown_shape():
    multi_index = da.from_array([[3, 6, 6], [4, 5, 1], [-1, -1, -1]])
    multi_index = multi_index[(multi_index > 0).all(axis=1)]

    multi_index_np = multi_index.compute()

    assert np.isnan(multi_index.shape).any()
    assert_eq(
        np.ravel_multi_index(multi_index_np, dims=(7, 6)),
        da.ravel_multi_index(multi_index, dims=(7, 6)),
    )


def test_ravel_multi_index_unknown_shape_fails():
    multi_index1 = da.from_array([2, -1, 3, -1], chunks=2)
    multi_index1 = multi_index1[multi_index1 > 0]

    multi_index2 = da.from_array(
        [[1, 2], [-1, -1], [3, 4], [5, 6], [7, 8], [-1, -1]], chunks=(2, 1)
    )
    multi_index2 = multi_index2[(multi_index2 > 0).all(axis=1)]

    multi_index = [1, multi_index1, multi_index2]

    assert np.isnan(multi_index1.shape).any()
    assert np.isnan(multi_index2.shape).any()
    with pytest.raises(ValueError, match="Arrays' chunk sizes"):
        da.ravel_multi_index(multi_index, dims=(8, 9, 10))


@pytest.mark.parametrize("dims", [da.from_array([5, 10]), delayed([5, 10], nout=2)])
@pytest.mark.parametrize("wrap_in_list", [False, True])
def test_ravel_multi_index_delayed_dims(dims, wrap_in_list):
    with pytest.raises(NotImplementedError, match="Dask types are not supported"):
        da.ravel_multi_index((2, 1), [dims[0], dims[1]] if wrap_in_list else dims)


def test_ravel_multi_index_non_int_dtype():
    with pytest.raises(TypeError, match="only int indices permitted"):
        da.ravel_multi_index(
            (1.0, 2),
            (5, 10),
        )


def test_coarsen():
    x = np.random.default_rng().integers(10, size=(24, 24))
    d = da.from_array(x, chunks=(4, 8))

    assert_eq(
        da.chunk.coarsen(np.sum, x, {0: 2, 1: 4}), da.coarsen(np.sum, d, {0: 2, 1: 4})
    )
    assert_eq(
        da.chunk.coarsen(np.sum, x, {0: 2, 1: 4}), da.coarsen(da.sum, d, {0: 2, 1: 4})
    )
    assert_eq(
        da.chunk.coarsen(np.mean, x, {0: 2, 1: 4}, dtype="float32"),
        da.coarsen(da.mean, d, {0: 2, 1: 4}, dtype="float32"),
    )


def test_coarsen_with_excess():
    x = da.arange(10, chunks=5)
    assert_eq(da.coarsen(np.min, x, {0: 5}, trim_excess=True), np.array([0, 5]))
    assert_eq(
        da.coarsen(np.sum, x, {0: 3}, trim_excess=True),
        np.array([0 + 1 + 2, 3 + 4 + 5, 6 + 7 + 8]),
    )


@pytest.mark.parametrize("chunks", [(x,) * 3 for x in range(16, 32)])
def test_coarsen_bad_chunks(chunks):
    x1 = da.arange(np.sum(chunks), chunks=5)
    x2 = x1.rechunk(tuple(chunks))
    assert_eq(
        da.coarsen(np.sum, x1, {0: 10}, trim_excess=True),
        da.coarsen(np.sum, x2, {0: 10}, trim_excess=True),
    )


@pytest.mark.parametrize(
    "chunks, divisor",
    [
        ((1, 1), 1),
        ((1, 1), 2),
        ((1, 1, 1), 2),
        ((10, 1), 10),
        ((20, 10, 15, 23, 24), 10),
        ((20, 10, 15, 23, 24), 8),
        ((10, 20, 30, 40, 2), 10),
        ((20, 10, 15, 42, 23, 24), 16),
        ((20, 10, 15, 47, 23, 24), 10),
        ((2, 10, 15, 47, 23, 24), 4),
    ],
)
def test_aligned_coarsen_chunks(chunks, divisor):
    from dask.array.routines import aligned_coarsen_chunks as acc

    aligned_chunks = acc(chunks, divisor)
    any_remainders = (np.array(aligned_chunks) % divisor) != 0
    valid_chunks = np.where((np.array(chunks) % divisor) == 0)[0]

    # check that total number of elements is conserved
    assert sum(aligned_chunks) == sum(chunks)
    # check that valid chunks are not modified
    assert [chunks[idx] for idx in valid_chunks] == [
        aligned_chunks[idx] for idx in valid_chunks
    ]
    # check that no chunks are 0
    assert (np.array(aligned_chunks) > 0).all()
    # check that at most one chunk was added
    assert len(aligned_chunks) <= len(chunks) + 1
    # check that either 0 or 1 chunks are not divisible by divisor
    assert any_remainders.sum() in (0, 1)
    # check that the only indivisible chunk is the last
    if any_remainders.sum() == 1:
        assert any_remainders[-1] == 1


def test_insert():
    rng = np.random.default_rng()
    x = rng.integers(10, size=(10, 10))
    a = da.from_array(x, chunks=(5, 5))
    y = rng.integers(10, size=(5, 10))
    b = da.from_array(y, chunks=(4, 4))

    assert_eq(np.insert(x, 0, -1, axis=0), da.insert(a, 0, -1, axis=0))
    assert_eq(np.insert(x, 3, -1, axis=-1), da.insert(a, 3, -1, axis=-1))
    assert_eq(np.insert(x, 5, -1, axis=1), da.insert(a, 5, -1, axis=1))
    assert_eq(np.insert(x, -1, -1, axis=-2), da.insert(a, -1, -1, axis=-2))
    assert_eq(np.insert(x, [2, 3, 3], -1, axis=1), da.insert(a, [2, 3, 3], -1, axis=1))
    assert_eq(
        np.insert(x, [2, 3, 8, 8, -2, -2], -1, axis=0),
        da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0),
    )
    assert_eq(
        np.insert(x, slice(1, 4), -1, axis=1), da.insert(a, slice(1, 4), -1, axis=1)
    )
    assert_eq(
        np.insert(x, [2] * 3 + [5] * 2, y, axis=0),
        da.insert(a, [2] * 3 + [5] * 2, b, axis=0),
    )
    assert_eq(np.insert(x, 0, y[0], axis=1), da.insert(a, 0, b[0], axis=1))

    assert same_keys(
        da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0),
        da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0),
    )

    with pytest.raises(NotImplementedError):
        da.insert(a, [4, 2], -1, axis=0)

    with pytest.raises(AxisError):
        da.insert(a, [3], -1, axis=2)

    with pytest.raises(AxisError):
        da.insert(a, [3], -1, axis=-3)


def test_append():
    rng = np.random.default_rng()
    x = rng.integers(10, size=(10, 10))
    a = da.from_array(x, chunks=(5, 5))

    # appendage for axis 1 / -1
    y1 = rng.integers(10, size=(10, 5))
    b1 = da.from_array(y1, chunks=(4, 4))

    # appendage for axis 0 / -2
    y0 = rng.integers(10, size=(5, 10))
    b0 = da.from_array(y0, chunks=(4, 4))

    # test axis None
    assert_eq(np.append(x, x, axis=None), da.append(a, a, axis=None))
    assert_eq(np.append(x, y0, axis=None), da.append(a, b0, axis=None))
    assert_eq(np.append(x, y1, axis=None), da.append(a, b1, axis=None))

    # test axis 0 / -2
    assert_eq(np.append(x, y0, axis=0), da.append(a, b0, axis=0))
    assert_eq(np.append(x, y0, axis=-2), da.append(a, b0, axis=-2))

    # test axis 1 / -1
    assert_eq(np.append(x, y1, axis=1), da.append(a, b1, axis=1))
    assert_eq(np.append(x, y1, axis=-1), da.append(a, b1, axis=-1))

    # test --> treat values as array_likes
    assert_eq(
        np.append(x, ((0,) * 10,) * 10, axis=None),
        da.append(a, ((0,) * 10,) * 10, axis=None),
    )
    assert_eq(
        np.append(x, ((0,) * 10,) * 10, axis=0), da.append(a, ((0,) * 10,) * 10, axis=0)
    )
    assert_eq(
        np.append(x, ((0,) * 10,) * 10, axis=1), da.append(a, ((0,) * 10,) * 10, axis=1)
    )

    # check AxisError
    with pytest.raises(AxisError):
        da.append(a, ((0,) * 10,) * 10, axis=2)
    with pytest.raises(AxisError):
        da.append(a, ((0,) * 10,) * 10, axis=-3)

    # check ValueError if dimensions don't align
    with pytest.raises(ValueError):
        da.append(a, (0,) * 10, axis=0)


def test_multi_insert():
    z = np.random.default_rng().integers(10, size=(1, 2))
    c = da.from_array(z, chunks=(1, 2))
    assert_eq(
        np.insert(np.insert(z, [0, 1], -1, axis=0), [1], -1, axis=1),
        da.insert(da.insert(c, [0, 1], -1, axis=0), [1], -1, axis=1),
    )


def test_delete():
    x = np.random.default_rng().integers(10, size=(10, 10))
    a = da.from_array(x, chunks=(5, 5))

    assert_eq(np.delete(x, 0, axis=0), da.delete(a, 0, axis=0))
    assert_eq(np.delete(x, 3, axis=-1), da.delete(a, 3, axis=-1))
    assert_eq(np.delete(x, 5, axis=1), da.delete(a, 5, axis=1))
    assert_eq(np.delete(x, -1, axis=-2), da.delete(a, -1, axis=-2))
    assert_eq(np.delete(x, [2, 3, 3], axis=1), da.delete(a, [2, 3, 3], axis=1))
    assert_eq(
        np.delete(x, [2, 3, 8, 8], axis=0),
        da.delete(a, [2, 3, 8, 8], axis=0),
    )
    assert_eq(np.delete(x, slice(1, 4), axis=1), da.delete(a, slice(1, 4), axis=1))
    assert_eq(
        np.delete(x, slice(1, 10, -1), axis=1), da.delete(a, slice(1, 10, -1), axis=1)
    )

    assert_eq(np.delete(a, [4, 2], axis=0), da.delete(a, [4, 2], axis=0))

    with pytest.raises(AxisError):
        da.delete(a, [3], axis=2)

    with pytest.raises(AxisError):
        da.delete(a, [3], axis=-3)


def test_result_type():
    a = da.from_array(np.ones(5, np.float32), chunks=(3,))
    b = da.from_array(np.ones(5, np.int16), chunks=(3,))
    c = da.from_array(np.ones(5, np.int64), chunks=(3,))
    x = np.ones(5, np.float32)
    assert da.result_type(b, c) == np.int64
    assert da.result_type(a, b, c) == np.float64
    assert da.result_type(b, np.float32) == np.float32
    assert da.result_type(b, np.dtype(np.float32)) == np.float32
    assert da.result_type(b, x) == np.float32
    # Effect of scalars depends on their value
    assert da.result_type(1, b) == np.int16
    assert da.result_type(1.0, a) == np.float32
    if NUMPY_GE_200:
        assert da.result_type(np.int64(1), b) == np.int64
        assert da.result_type(np.ones((), np.int64), b) == np.int64
        assert da.result_type(1e200, a) == np.float32
    else:
        assert da.result_type(np.int64(1), b) == np.int16
        assert da.result_type(np.ones((), np.int64), b) == np.int16  # 0d array
        assert da.result_type(1e200, a) == np.float64  # 1e200 is too big for float32

    # dask 0d-arrays are NOT treated like scalars
    c = da.from_array(np.ones((), np.float64), chunks=())
    assert da.result_type(a, c) == np.float64


def _numpy_and_dask_inputs(input_sigs):
    # einsum label dimensions
    _dimensions = {
        "a": 5,
        "b": 6,
        "c": 7,
        "d": 5,
        "e": 6,
        "f": 10,
        "g": 1,
        "h": 2,
        "*": 11,
    }

    # dimension chunks sizes
    _chunks = {
        "a": (2, 3),
        "b": (2, 3, 1),
        "c": (2, 3, 2),
        "d": (4, 1),
        "e": (2, 4),
        "f": (1, 2, 3, 4),
        "g": 1,
        "h": (1, 1),
        "*": 11,
    }

    def _shape_from_string(s):
        return tuple(_dimensions[c] for c in s)

    def _chunks_from_string(s):
        return tuple(_chunks[c] for c in s)

    shapes = [_shape_from_string(s) for s in input_sigs]
    chunks = [_chunks_from_string(s) for s in input_sigs]

    np_inputs = [np.random.default_rng().random(s) for s in shapes]
    da_inputs = [da.from_array(i, chunks=c) for i, c in zip(np_inputs, chunks)]

    return np_inputs, da_inputs


@pytest.mark.parametrize(
    "einsum_signature",
    [
        "abc,bad->abcd",
        "abcdef,bcdfg->abcdeg",
        "ea,fb,abcd,gc,hd->efgh",
        "ab,b",
        "aa",
        "a,a->",
        "a,a->a",
        "a,a",
        "a,b",
        "a,b,c",
        "a",
        "ba,b",
        "ba,b->",
        "defab,fedbc->defac",
        "ab...,bc...->ac...",
        "a...a",
        "abc...->cba...",
        "...ab->...a",
        "a...a->a...",
        # Following 2 from # https://stackoverflow.com/a/19203475/1611416
        "...abc,...abcd->...d",
        "ab...,b->ab...",
        # https://github.com/dask/dask/pull/3412#discussion_r182413444
        "aa->a",
        "ab,ab,c->c",
        "aab,bc->ac",
        "aab,bcc->ac",
        "fdf,cdd,ccd,afe->ae",
        "fff,fae,bef,def->abd",
    ],
)
def test_einsum(einsum_signature):
    input_sigs = einsum_signature.split("->")[0].replace("...", "*").split(",")

    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=da.PerformanceWarning)
        assert_eq(
            np.einsum(einsum_signature, *np_inputs),
            da.einsum(einsum_signature, *da_inputs),
        )


def test_einsum_chunksizes():
    arr1 = da.random.random((1024, 8, 8, 8, 8), chunks=(256, 8, 8, 8, 8))
    arr2 = da.random.random((1024, 8, 8, 8, 8), chunks=(256, 8, 8, 8, 8))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=da.PerformanceWarning)
        result = da.einsum("aijkl,amnop->ijklmnop", arr1, arr2)
    assert result.chunks == ((4,) * 2,) * 8

    arr1 = da.random.random((64, 8, 8, 8, 8), chunks=(32, 8, 1, 8, 8))
    arr2 = da.random.random((64, 8, 8, 8, 8), chunks=(32, 8, 8, 1, 8))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=da.PerformanceWarning)
        result = da.einsum("aijkl,amnop->ijklmnop", arr1, arr2)
    assert result.chunks == (
        (4,) * 2,
        (1,) * 8,
        (4,) * 2,
        (4,) * 2,
        (4,) * 2,
        (4,) * 2,
        (1,) * 8,
        (4,) * 2,
    )

    np_arr1 = np.random.random((2, 4, 4))
    np_arr2 = np.random.random((2, 4, 4))

    arr1 = da.from_array(np_arr1, chunks=(1, 2, 2))
    arr2 = da.from_array(np_arr2, chunks=(1, 2, 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=da.PerformanceWarning)
        result = da.einsum("aij,amn->ijmn", arr1, arr2)
    assert result.chunks == ((1,) * 4,) * 4
    assert_eq(np.einsum("aij,amn->ijmn", np_arr1, np_arr2), result)

    # regression test for GH11627
    z = da.ones(
        shape=(40000, 2, 10, 2, 10), dtype=np.float64, chunksize=(40000, 1, 5, 2, 10)
    )
    x = da.ones(shape=(2, 10, 10), dtype=np.float64, chunksize=(2, 10, 10))
    y = da.ones(shape=(2, 10, 10), dtype=np.float64, chunksize=(2, 10, 10))
    res = da.einsum("abcde,bfc,dfe->acef", z, x, y)
    assert res.numblocks == (1, 1, 1, 1)


@pytest.mark.parametrize(
    "optimize_opts", [(True, False), ("greedy", False), ("optimal", False)]
)
def test_einsum_optimize(optimize_opts):
    sig = "ea,fb,abcd,gc,hd->efgh"
    input_sigs = sig.split("->")[0].split(",")
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    opt1, opt2 = optimize_opts

    assert_eq(
        np.einsum(sig, *np_inputs, optimize=opt1),
        da.einsum(sig, *np_inputs, optimize=opt2),
    )

    assert_eq(
        np.einsum(sig, *np_inputs, optimize=opt2),
        da.einsum(sig, *np_inputs, optimize=opt1),
    )


@pytest.mark.parametrize("order", ["C", "F", "A", "K"])
def test_einsum_order(order):
    sig = "ea,fb,abcd,gc,hd->efgh"
    input_sigs = sig.split("->")[0].split(",")
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    assert_eq(
        np.einsum(sig, *np_inputs, order=order), da.einsum(sig, *np_inputs, order=order)
    )


@pytest.mark.parametrize("casting", ["no", "equiv", "safe", "same_kind", "unsafe"])
def test_einsum_casting(casting):
    sig = "ea,fb,abcd,gc,hd->efgh"
    input_sigs = sig.split("->")[0].split(",")
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    assert_eq(
        np.einsum(sig, *np_inputs, casting=casting),
        da.einsum(sig, *np_inputs, casting=casting),
    )


@pytest.mark.parametrize("split_every", [None, 2])
def test_einsum_split_every(split_every):
    np_inputs, da_inputs = _numpy_and_dask_inputs("a")
    assert_eq(
        np.einsum("a", *np_inputs), da.einsum("a", *da_inputs, split_every=split_every)
    )


def test_einsum_invalid_args():
    _, da_inputs = _numpy_and_dask_inputs("a")
    with pytest.raises(TypeError):
        da.einsum("a", *da_inputs, foo=1, bar=2)


def test_einsum_broadcasting_contraction():
    rng = np.random.default_rng()
    a = rng.random((1, 5, 4))
    b = rng.random((4, 6))
    c = rng.random((5, 6))
    d = rng.random(10)

    d_a = da.from_array(a, chunks=(1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum("ijk,kl,jl", a, b, c)
    da_res = da.einsum("ijk,kl,jl", d_a, d_b, d_c)
    assert_eq(np_res, da_res)

    mul_res = da_res * d

    np_res = np.einsum("ijk,kl,jl,i->i", a, b, c, d)
    da_res = da.einsum("ijk,kl,jl,i->i", d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)
    assert_eq(np_res, mul_res)


def test_einsum_broadcasting_contraction2():
    rng = np.random.default_rng()
    a = rng.random((1, 1, 5, 4))
    b = rng.random((4, 6))
    c = rng.random((5, 6))
    d = rng.random((7, 7))

    d_a = da.from_array(a, chunks=(1, 1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum("abjk,kl,jl", a, b, c)
    da_res = da.einsum("abjk,kl,jl", d_a, d_b, d_c)
    assert_eq(np_res, da_res)

    mul_res = da_res * d

    np_res = np.einsum("abjk,kl,jl,ab->ab", a, b, c, d)
    da_res = da.einsum("abjk,kl,jl,ab->ab", d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)
    assert_eq(np_res, mul_res)


def test_einsum_broadcasting_contraction3():
    rng = np.random.default_rng()
    a = rng.random((1, 5, 4))
    b = rng.random((4, 1, 6))
    c = rng.random((5, 6))
    d = rng.random((7, 7))

    d_a = da.from_array(a, chunks=(1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), 1, (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum("ajk,kbl,jl,ab->ab", a, b, c, d)
    da_res = da.einsum("ajk,kbl,jl,ab->ab", d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)


def test_einsum_empty_dimension():
    arr = np.random.random((10, 10))
    darr = da.from_array(arr, chunks=(5, 5))
    darr = darr[:0]
    result = da.einsum("ca,ca->c", darr, darr)
    assert_eq(result, np.einsum("ca,ca->c", arr[:0], arr[:0]))


@pytest.mark.parametrize("a", [np.arange(11), np.arange(6).reshape((3, 2))])
@pytest.mark.parametrize("returned", [True, False])
def test_average(a, returned):
    d_a = da.from_array(a, chunks=2)

    np_avg = np.average(a, returned=returned)
    da_avg = da.average(d_a, returned=returned)

    assert_eq(np_avg, da_avg)


@pytest.mark.parametrize("a", [np.arange(11), np.arange(6).reshape((3, 2))])
def test_average_keepdims(a):
    d_a = da.from_array(a, chunks=2)

    da_avg = da.average(d_a, keepdims=True)

    np_avg = np.average(a, keepdims=True)
    assert_eq(np_avg, da_avg)


@pytest.mark.parametrize("keepdims", [False, True])
def test_average_weights(keepdims):
    a = np.arange(6).reshape((3, 2))
    d_a = da.from_array(a, chunks=2)

    weights = np.array([0.25, 0.75])
    d_weights = da.from_array(weights, chunks=2)

    da_avg = da.average(d_a, weights=d_weights, axis=1, keepdims=keepdims)

    assert_eq(da_avg, np.average(a, weights=weights, axis=1, keepdims=keepdims))


def test_average_raises():
    d_a = da.arange(11, chunks=2)

    with pytest.raises(TypeError):
        da.average(d_a, weights=[1, 2, 3])

    with pytest.warns(RuntimeWarning):
        da.average(d_a, weights=da.zeros_like(d_a)).compute()


def test_iscomplexobj():
    a = da.from_array(np.array([1, 2]), 2)
    assert np.iscomplexobj(a) is False

    a = da.from_array(np.array([1, 2 + 0j]), 2)
    assert np.iscomplexobj(a) is True


def test_tril_triu():
    A = np.random.default_rng().standard_normal((20, 20))
    for chk in [5, 4]:
        dA = da.from_array(A, (chk, chk))

        assert np.allclose(da.triu(dA).compute(), np.triu(A))
        assert np.allclose(da.tril(dA).compute(), np.tril(A))

        for k in [
            -25,
            -20,
            -19,
            -15,
            -14,
            -9,
            -8,
            -6,
            -5,
            -1,
            1,
            4,
            5,
            6,
            8,
            10,
            11,
            15,
            16,
            19,
            20,
            21,
        ]:
            assert np.allclose(da.triu(dA, k).compute(), np.triu(A, k))
            assert np.allclose(da.tril(dA, k).compute(), np.tril(A, k))


def test_tril_ndims():
    A = np.random.default_rng().integers(0, 11, (10, 10, 10))
    dA = da.from_array(A, chunks=(5, 5, 5))
    assert_eq(da.triu(dA), np.triu(A))


def test_tril_triu_non_square_arrays():
    A = np.random.default_rng().integers(0, 11, (30, 35))
    dA = da.from_array(A, chunks=(5, 5))
    assert_eq(da.triu(dA), np.triu(A))
    assert_eq(da.tril(dA), np.tril(A))


@pytest.mark.parametrize(
    "n, k, m, chunks",
    [(3, 0, 3, "auto"), (3, 1, 3, "auto"), (3, -1, 3, "auto"), (5, 0, 5, 1)],
)
def test_tril_triu_indices(n, k, m, chunks):
    actual = da.tril_indices(n=n, k=k, m=m, chunks=chunks)[0]
    expected = np.tril_indices(n=n, k=k, m=m)[0]

    if sys.platform == "win32":
        assert_eq(
            actual.astype(expected.dtype),
            expected,
        )
    else:
        assert_eq(actual, expected)

    actual = da.triu_indices(n=n, k=k, m=m, chunks=chunks)[0]
    expected = np.triu_indices(n=n, k=k, m=m)[0]

    if sys.platform == "win32":
        assert_eq(
            actual.astype(expected.dtype),
            expected,
        )
    else:
        assert_eq(actual, expected)


def test_pickle_vectorized_routines():
    """Test that graphs that internally use np.vectorize can be pickled"""
    a = da.from_array(["foo", "bar", ""])

    b = da.count_nonzero(a)
    assert_eq(b, 2, check_dtype=False)
    b2 = pickle.loads(pickle.dumps(b))
    assert_eq(b2, 2, check_dtype=False)

    c = da.argwhere(a)
    assert_eq(c, [[0], [1]], check_dtype=False)
    c2 = pickle.loads(pickle.dumps(c))
    assert_eq(c2, [[0], [1]], check_dtype=False)
