from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_equal

import dask.array as da
from dask.array.core import Array
from dask.array.gufunc import (
    _parse_gufunc_signature,
    _validate_normalize_axes,
    apply_gufunc,
    as_gufunc,
    gufunc,
)
from dask.array.utils import assert_eq


# Copied from `numpy.lib.test_test_function_base.py`:
def test__parse_gufunc_signature():
    assert_equal(_parse_gufunc_signature("(x)->()"), ([("x",)], ()))
    assert_equal(_parse_gufunc_signature("(x,y)->()"), ([("x", "y")], ()))
    # whitespace
    assert_equal(_parse_gufunc_signature("  (x, y) ->()"), ([("x", "y")], ()))
    assert_equal(_parse_gufunc_signature("(x),(y)->()"), ([("x",), ("y",)], ()))
    assert_equal(_parse_gufunc_signature("(x)->(y)"), ([("x",)], ("y",)))
    assert_equal(_parse_gufunc_signature("(x)->(y),()"), ([("x",)], [("y",), ()]))
    assert_equal(
        _parse_gufunc_signature("(),(a,b,c),(d)->(d,e)"),
        ([(), ("a", "b", "c"), ("d",)], ("d", "e")),
    )
    with pytest.raises(ValueError):
        _parse_gufunc_signature("(x)(y)->()")
    with pytest.raises(ValueError):
        _parse_gufunc_signature("(x),(y)->")
    with pytest.raises(ValueError):
        _parse_gufunc_signature("((x))->(x)")
    with pytest.raises(ValueError):
        _parse_gufunc_signature("(x)->(x),")


def test_apply_gufunc_axes_input_validation_01():
    def foo(x):
        return np.mean(x, axis=-1)

    a = da.random.default_rng().normal(size=(20, 30), chunks=30)

    with pytest.raises(ValueError):
        apply_gufunc(foo, "(i)->()", a, axes=0)

    apply_gufunc(foo, "(i)->()", a, axes=[0])
    apply_gufunc(foo, "(i)->()", a, axes=[(0,)])
    apply_gufunc(foo, "(i)->()", a, axes=[0, tuple()])
    apply_gufunc(foo, "(i)->()", a, axes=[(0,), tuple()])

    with pytest.raises(ValueError):
        apply_gufunc(foo, "(i)->()", a, axes=[(0, 1)])

    with pytest.raises(ValueError):
        apply_gufunc(foo, "(i)->()", a, axes=[0, 0])


def test_apply_gufunc_axes_args_validation():
    def add(x, y):
        return x + y

    a = da.from_array(np.array([1, 2, 3]), chunks=2, name="a")
    b = da.from_array(np.array([1, 2, 3]), chunks=2, name="b")
    with pytest.raises(ValueError):
        apply_gufunc(add, "(),()->()", a, b, 0, output_dtypes=a.dtype)


def test__validate_normalize_axes_01():
    with pytest.raises(ValueError):
        _validate_normalize_axes([(1, 0)], None, False, [("i", "j")], ("j",))

    with pytest.raises(ValueError):
        _validate_normalize_axes([0, 0], None, False, [("i", "j")], ("j",))

    with pytest.raises(ValueError):
        _validate_normalize_axes([(0,), 0], None, False, [("i", "j")], ("j",))

    i, o = _validate_normalize_axes([(1, 0), 0], None, False, [("i", "j")], ("j",))
    assert i == [(1, 0)]
    assert o == [(0,)]


def test__validate_normalize_axes_02():
    i, o = _validate_normalize_axes(None, 0, False, [("i",), ("i",)], ())
    assert i == [(0,), (0,)]
    assert o == [()]

    i, o = _validate_normalize_axes(None, 0, False, [("i",)], ("i",))
    assert i == [(0,)]
    assert o == [(0,)]

    i, o = _validate_normalize_axes(None, 0, True, [("i",), ("i",)], ())
    assert i == [(0,), (0,)]
    assert o == [(0,)]

    with pytest.raises(ValueError):
        _validate_normalize_axes(None, (0,), False, [("i",), ("i",)], ())

    with pytest.raises(ValueError):
        _validate_normalize_axes(None, 0, False, [("i",), ("j",)], ())

    with pytest.raises(ValueError):
        _validate_normalize_axes(None, 0, False, [("i",), ("j",)], ("j",))


def test__validate_normalize_axes_03():
    i, o = _validate_normalize_axes(None, 0, True, [("i",)], ())
    assert i == [(0,)]
    assert o == [(0,)]

    with pytest.raises(ValueError):
        _validate_normalize_axes(None, 0, True, [("i",)], ("i",))

    with pytest.raises(ValueError):
        _validate_normalize_axes([(0, 1), (0, 1)], None, True, [("i", "j")], ("i", "j"))

    with pytest.raises(ValueError):
        _validate_normalize_axes([(0,), (0,)], None, True, [("i",), ("j",)], ())


def test_apply_gufunc_01():
    def stats(x):
        return np.mean(x, axis=-1), np.std(x, axis=-1)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=(5, 5, 30))
    result = apply_gufunc(stats, "(i)->(),()", a, output_dtypes=2 * (a.dtype,))
    mean, std = result
    assert isinstance(result, tuple)
    assert mean.compute().shape == (10, 20)
    assert std.compute().shape == (10, 20)


def test_apply_gufunc_01b():
    def stats(x):
        return np.mean(x, axis=-1), np.std(x, axis=-1)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=5)
    mean, std = apply_gufunc(
        stats, "(i)->(),()", a, output_dtypes=2 * (a.dtype,), allow_rechunk=True
    )
    assert mean.compute().shape == (10, 20)
    assert std.compute().shape == (10, 20)


@pytest.mark.parametrize("vectorize", [False, True])
def test_apply_gufunc_output_dtypes_string(vectorize):
    def stats(x):
        return np.mean(x, axis=-1)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=(5, 5, 30))
    mean = apply_gufunc(stats, "(i)->()", a, output_dtypes="f", vectorize=vectorize)
    assert mean.compute().shape == (10, 20)


@pytest.mark.parametrize("vectorize", [False, True])
def test_apply_gufunc_output_dtypes_string_many_outputs(vectorize):
    def stats(x):
        return np.mean(x, axis=-1), np.std(x, axis=-1), np.min(x, axis=-1)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=(5, 5, 30))
    mean, std, min = apply_gufunc(
        stats, "(i)->(),(),()", a, output_dtypes=("f", "f", "f"), vectorize=vectorize
    )
    assert mean.compute().shape == (10, 20)
    assert std.compute().shape == (10, 20)
    assert min.compute().shape == (10, 20)


def test_apply_gufunc_pass_additional_kwargs():
    def foo(x, bar):
        assert bar == 2
        return x

    ret = apply_gufunc(foo, "()->()", 1.0, output_dtypes=float, bar=2)
    assert_eq(ret, np.array(1.0, dtype=float))


def test_apply_gufunc_02():
    def outer_product(x, y):
        return np.einsum("...i,...j->...ij", x, y)

    rng = da.random.default_rng()
    a = rng.normal(size=(20, 30), chunks=(5, 30))
    b = rng.normal(size=(10, 1, 40), chunks=(10, 1, 40))
    c = apply_gufunc(outer_product, "(i),(j)->(i,j)", a, b, output_dtypes=a.dtype)

    assert c.compute().shape == (10, 20, 30, 40)


def test_apply_gufunc_scalar_output():
    def foo():
        return 1

    x = apply_gufunc(foo, "->()", output_dtypes=int)
    assert x.compute() == 1


def test_apply_gufunc_elemwise_01():
    def add(x, y):
        return x + y

    a = da.from_array(np.array([1, 2, 3]), chunks=2, name="a")
    b = da.from_array(np.array([1, 2, 3]), chunks=2, name="b")
    z = apply_gufunc(add, "(),()->()", a, b, output_dtypes=a.dtype)
    assert_eq(z, np.array([2, 4, 6]))


def test_apply_gufunc_elemwise_01b():
    def add(x, y):
        return x + y

    a = da.from_array(np.array([1, 2, 3]), chunks=2, name="a")
    b = da.from_array(np.array([1, 2, 3]), chunks=1, name="b")
    with pytest.raises(ValueError):
        apply_gufunc(add, "(),()->()", a, b, output_dtypes=a.dtype)


def test_apply_gufunc_elemwise_02():
    def addmul(x, y):
        assert x.shape in ((2,), (1,))
        return x + y, x * y

    a = da.from_array(np.array([1, 2, 3]), chunks=2, name="a")
    b = da.from_array(np.array([1, 2, 3]), chunks=2, name="b")
    z1, z2 = apply_gufunc(addmul, "(),()->(),()", a, b, output_dtypes=2 * (a.dtype,))
    assert_eq(z1, np.array([2, 4, 6]))
    assert_eq(z2, np.array([1, 4, 9]))


def test_gufunc_vector_output():
    def foo():
        return np.array([1, 2, 3], dtype=int)

    x = apply_gufunc(foo, "->(i_0)", output_dtypes=int, output_sizes={"i_0": 3})
    assert x.chunks == ((3,),)
    assert_eq(x, np.array([1, 2, 3]))


def test_apply_gufunc_elemwise_loop():
    def foo(x):
        assert x.shape in ((2,), (1,))
        return 2 * x

    a = da.from_array(np.array([1, 2, 3]), chunks=2, name="a")
    z = apply_gufunc(foo, "()->()", a, output_dtypes=int)
    assert z.chunks == ((2, 1),)
    assert_eq(z, np.array([2, 4, 6]))


def test_apply_gufunc_elemwise_core():
    def foo(x):
        assert x.shape == (3,)
        return 2 * x

    a = da.from_array(np.array([1, 2, 3]), chunks=3, name="a")
    z = apply_gufunc(foo, "(i)->(i)", a, output_dtypes=int)
    assert z.chunks == ((3,),)
    assert_eq(z, np.array([2, 4, 6]))


# TODO: In case single tuple output will get enabled:
# def test_apply_gufunc_one_scalar_output():
#     def foo():
#         return 1,
#     x, = apply_gufunc(foo, "->(),", output_dtypes=(int,))
#     assert x.compute() == 1


def test_apply_gufunc_two_scalar_output():
    def foo():
        return 1, 2

    x, y = apply_gufunc(foo, "->(),()", output_dtypes=(int, int))
    assert x.compute() == 1
    assert y.compute() == 2


def test_apply_gufunc_two_mixed_outputs():
    def foo():
        return 1, np.ones((2, 3), dtype=float)

    x, y = apply_gufunc(
        foo, "->(),(i,j)", output_dtypes=(int, float), output_sizes={"i": 2, "j": 3}
    )
    assert x.compute() == 1
    assert y.chunks == ((2,), (3,))
    assert_eq(y, np.ones((2, 3), dtype=float))


@pytest.mark.parametrize("output_dtypes", [int, (int,)])
def test_apply_gufunc_output_dtypes(output_dtypes):
    def foo(x):
        return y

    x = np.random.default_rng().standard_normal(10)
    y = x.astype(int)
    dy = apply_gufunc(foo, "()->()", x, output_dtypes=output_dtypes)
    # print(x, x.compute())
    assert_eq(y, dy)


def test_gufunc_two_inputs():
    def foo(x, y):
        return np.einsum("...ij,...jk->ik", x, y)

    a = da.ones((2, 3), chunks=100, dtype=int)
    b = da.ones((3, 4), chunks=100, dtype=int)
    x = apply_gufunc(foo, "(i,j),(j,k)->(i,k)", a, b, output_dtypes=int)
    assert_eq(x, 3 * np.ones((2, 4), dtype=int))


def test_gufunc_mixed_inputs():
    def foo(x, y):
        return x + y

    a = np.ones((2, 1), dtype=int)
    b = da.ones((1, 8), chunks=(2, 3), dtype=int)
    x = apply_gufunc(foo, "(),()->()", a, b, output_dtypes=int)
    assert_eq(x, 2 * np.ones((2, 8), dtype=int))


def test_gufunc_mixed_inputs_vectorize():
    def foo(x, y):
        return (x + y).sum(axis=1)

    a = da.ones((8, 3, 5), chunks=(2, 3, 5), dtype=int)
    b = np.ones(5, dtype=int)
    x = apply_gufunc(foo, "(m,n),(n)->(m)", a, b, vectorize=True)

    assert_eq(x, np.full((8, 3), 10, dtype=int))


def test_gufunc_vectorize_whitespace():
    # Regression test for https://github.com/dask/dask/issues/7972.
    # NumPy versions before https://github.com/numpy/numpy/pull/19627
    # would not ignore whitespace characters in `signature` like they
    # are supposed to. We remove the whitespace in Dask as a workaround.

    def foo(x, y):
        return (x + y).sum(axis=1)

    a = da.ones((8, 3, 5), chunks=(2, 3, 5), dtype=int)
    b = np.ones(5, dtype=int)
    x = apply_gufunc(foo, "(m, n),(n)->(m)", a, b, vectorize=True)

    assert_eq(x, np.full((8, 3), 10, dtype=int))

    a = da.random.default_rng().random((6, 5, 5))

    @da.as_gufunc(signature="(n, n)->(n, n)", output_dtypes=float, vectorize=True)
    def gufoo(x):
        return np.linalg.inv(x)

    # Previously calling `gufoo` would raise an error due to the whitespace
    # in its `signature`. Let's make sure it doesn't raise here.
    gufoo(a)


def test_gufunc():
    x = da.random.default_rng().normal(size=(10, 5), chunks=(2, 5))

    def foo(x):
        return np.mean(x, axis=-1)

    gufoo = gufunc(
        foo,
        signature="(i)->()",
        axis=-1,
        keepdims=False,
        output_dtypes=float,
        vectorize=True,
    )

    y = gufoo(x)
    valy = y.compute()
    assert isinstance(y, Array)
    assert valy.shape == (10,)


def test_as_gufunc():
    x = da.random.default_rng().normal(size=(10, 5), chunks=(2, 5))

    @as_gufunc("(i)->()", axis=-1, keepdims=False, output_dtypes=float, vectorize=True)
    def foo(x):
        return np.mean(x, axis=-1)

    y = foo(x)
    valy = y.compute()
    assert isinstance(y, Array)
    assert valy.shape == (10,)


def test_apply_gufunc_broadcasting_loopdims():
    def foo(x, y):
        assert len(x.shape) == 2
        assert len(y.shape) == 3
        x, y = np.broadcast_arrays(x, y)
        return x, y, x * y

    rng = da.random.default_rng()
    a = rng.normal(size=(10, 30), chunks=(8, 30))
    b = rng.normal(size=(20, 1, 30), chunks=(3, 1, 30))

    x, y, z = apply_gufunc(
        foo, "(i),(i)->(i),(i),(i)", a, b, output_dtypes=3 * (float,), vectorize=False
    )

    assert x.compute().shape == (20, 10, 30)
    assert y.compute().shape == (20, 10, 30)
    assert z.compute().shape == (20, 10, 30)


def test_apply_gufunc_check_same_dimsizes():
    def foo(x, y):
        return x + y

    rng = da.random.default_rng()
    a = rng.normal(size=(3,), chunks=(2,))
    b = rng.normal(size=(4,), chunks=(2,))

    with pytest.raises(ValueError) as excinfo:
        apply_gufunc(foo, "(),()->()", a, b, output_dtypes=float, allow_rechunk=True)
    assert "different lengths in arrays" in str(excinfo.value)


def test_apply_gufunc_check_coredim_chunksize():
    def foo(x):
        return np.sum(x, axis=-1)

    a = da.random.default_rng().normal(size=(8,), chunks=3)
    with pytest.raises(ValueError) as excinfo:
        da.apply_gufunc(foo, "(i)->()", a, output_dtypes=float, allow_rechunk=False)
    assert "consists of multiple chunks" in str(excinfo.value)


def test_apply_gufunc_check_inhomogeneous_chunksize():
    def foo(x, y):
        return x + y

    rng = da.random.default_rng()
    a = rng.normal(size=(8,), chunks=((2, 2, 2, 2),))
    b = rng.normal(size=(8,), chunks=((2, 3, 3),))

    with pytest.raises(ValueError) as excinfo:
        da.apply_gufunc(
            foo, "(),()->()", a, b, output_dtypes=float, allow_rechunk=False
        )
    assert "with different chunksize present" in str(excinfo.value)


def test_apply_gufunc_infer_dtype():
    x = np.arange(50).reshape((5, 10))
    y = np.arange(10)
    dx = da.from_array(x, chunks=5)
    dy = da.from_array(y, chunks=5)

    def foo(x, *args, **kwargs):
        cast = kwargs.pop("cast", "i8")
        return (x + sum(args)).astype(cast)

    dz = apply_gufunc(foo, "(),(),()->()", dx, dy, 1)
    z = foo(dx, dy, 1)
    assert_eq(dz, z)

    dz = apply_gufunc(foo, "(),(),()->()", dx, dy, 1, cast="f8")
    z = foo(dx, dy, 1, cast="f8")
    assert_eq(dz, z)

    dz = apply_gufunc(foo, "(),(),()->()", dx, dy, 1, cast="f8", output_dtypes="f8")
    z = foo(dx, dy, 1, cast="f8")
    assert_eq(dz, z)

    def foo(x):
        raise RuntimeError("Woops")

    with pytest.raises(ValueError) as e:
        apply_gufunc(foo, "()->()", dx)
    msg = str(e.value)
    assert msg.startswith("`dtype` inference failed")
    assert "Please specify the dtype explicitly" in msg
    assert "RuntimeError" in msg

    # Multiple outputs
    def foo(x, y):
        return x + y, x - y

    z0, z1 = apply_gufunc(foo, "(),()->(),()", dx, dy)

    assert_eq(z0, dx + dy)
    assert_eq(z1, dx - dy)


@pytest.mark.parametrize("keepdims", [False, True])
def test_apply_gufunc_axis_01(keepdims):
    def mymedian(x):
        return np.median(x, axis=-1)

    a = np.random.default_rng().standard_normal((10, 5))
    da_ = da.from_array(a, chunks=2)

    m = np.median(a, axis=0, keepdims=keepdims)
    dm = apply_gufunc(
        mymedian, "(i)->()", da_, axis=0, keepdims=keepdims, allow_rechunk=True
    )
    assert_eq(m, dm)


def test_apply_gufunc_axis_02():
    def myfft(x):
        return np.fft.fft(x, axis=-1)

    a = np.random.default_rng().standard_normal((10, 5))
    da_ = da.from_array(a, chunks=2)

    m = np.fft.fft(a, axis=0)
    dm = apply_gufunc(myfft, "(i)->(i)", da_, axis=0, allow_rechunk=True)
    assert_eq(m, dm)


def test_apply_gufunc_axis_02b():
    def myfilter(x, cn=10, axis=-1):
        y = np.fft.fft(x, axis=axis)
        y[cn:-cn] = 0
        nx = np.fft.ifft(y, axis=axis)
        return np.real(nx)

    a = np.random.default_rng().standard_normal((3, 6, 4))
    da_ = da.from_array(a, chunks=2)

    m = myfilter(a, axis=1)
    dm = apply_gufunc(myfilter, "(i)->(i)", da_, axis=1, allow_rechunk=True)
    assert_eq(m, dm)


def test_apply_gufunc_axis_03():
    def mydiff(x):
        return np.diff(x, axis=-1)

    a = np.random.default_rng().standard_normal((3, 6, 4))
    da_ = da.from_array(a, chunks=2)

    m = np.diff(a, axis=1)
    dm = apply_gufunc(
        mydiff, "(i)->(i)", da_, axis=1, output_sizes={"i": 5}, allow_rechunk=True
    )
    assert_eq(m, dm)


@pytest.mark.parametrize("axis", [-2, -1, None])
def test_apply_gufunc_axis_keepdims(axis):
    def mymedian(x):
        return np.median(x, axis=-1)

    a = np.random.default_rng().standard_normal((10, 5))
    da_ = da.from_array(a, chunks=2)

    m = np.median(a, axis=-1 if not axis else axis, keepdims=True)
    dm = apply_gufunc(
        mymedian, "(i)->()", da_, axis=axis, keepdims=True, allow_rechunk=True
    )
    assert_eq(m, dm)


@pytest.mark.parametrize("axes", [[0, 1], [(0,), (1,)]])
def test_apply_gufunc_axes_01(axes):
    def mystats(x, y):
        return np.std(x, axis=-1) * np.mean(y, axis=-1)

    rng = np.random.default_rng()
    a = rng.standard_normal((10, 5))
    b = rng.standard_normal((5, 6))
    da_ = da.from_array(a, chunks=2)
    db_ = da.from_array(b, chunks=2)

    m = np.std(a, axis=0) * np.mean(b, axis=1)
    dm = apply_gufunc(mystats, "(i),(j)->()", da_, db_, axes=axes, allow_rechunk=True)
    assert_eq(m, dm)


def test_apply_gufunc_axes_02():
    def matmul(x, y):
        return np.einsum("...ij,...jk->...ik", x, y)

    rng = np.random.default_rng()
    a = rng.standard_normal((3, 2, 1))
    b = rng.standard_normal((3, 7, 5))

    da_ = da.from_array(a, chunks=2)
    db = da.from_array(b, chunks=3)

    m = np.einsum("jiu,juk->uik", a, b)
    dm = apply_gufunc(
        matmul,
        "(i,j),(j,k)->(i,k)",
        da_,
        db,
        axes=[(1, 0), (0, -1), (-2, -1)],
        allow_rechunk=True,
    )
    assert_eq(m, dm)


def test_apply_gufunc_axes_two_kept_coredims():
    rng = da.random.default_rng()
    a = rng.normal(size=(20, 30), chunks=(10, 30))
    b = rng.normal(size=(10, 1, 40), chunks=(5, 1, 40))

    def outer_product(x, y):
        return np.einsum("i,j->ij", x, y)

    c = apply_gufunc(outer_product, "(i),(j)->(i,j)", a, b, vectorize=True)
    assert c.compute().shape == (10, 20, 30, 40)


def test_apply_gufunc_via_numba_01():
    numba = pytest.importorskip("numba")

    @numba.guvectorize(
        [(numba.float64[:], numba.float64[:], numba.float64[:])], "(n),(n)->(n)"
    )
    def g(x, y, res):
        for i in range(x.shape[0]):
            res[i] = x[i] + y[i]

    rng = da.random.default_rng()
    a = rng.normal(size=(20, 30), chunks=30)
    b = rng.normal(size=(20, 30), chunks=30)

    x = a + b
    y = g(a, b, axis=0)

    assert_eq(x, y)


def test_apply_gufunc_via_numba_02():
    numba = pytest.importorskip("numba")

    @numba.guvectorize([(numba.float64[:], numba.float64[:])], "(n)->()")
    def mysum(x, res):
        res[0] = 0.0
        for i in range(x.shape[0]):
            res[0] += x[i]

    a = da.random.default_rng().normal(size=(20, 30), chunks=30)

    x = a.sum(axis=0, keepdims=True)
    y = mysum(a, axis=0, keepdims=True)
    assert_eq(x, y)


def test_preserve_meta_type():
    sparse = pytest.importorskip("sparse")

    def stats(x):
        return np.sum(x, axis=-1), np.mean(x, axis=-1)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=(5, 5, 30))
    a = a.map_blocks(sparse.COO.from_numpy)
    sum, mean = apply_gufunc(stats, "(i)->(),()", a, output_dtypes=2 * (a.dtype,))

    assert isinstance(a._meta, sparse.COO)
    assert isinstance(sum._meta, sparse.COO)
    assert isinstance(mean._meta, sparse.COO)

    assert_eq(sum, sum)
    assert_eq(mean, mean)


def test_apply_gufunc_with_meta():
    def stats(x):
        return np.mean(x, axis=-1), np.std(x, axis=-1, dtype=np.float32)

    a = da.random.default_rng().normal(size=(10, 20, 30), chunks=(5, 5, 30))
    meta = (np.ones(0, dtype=np.float64), np.ones(0, dtype=np.float32))
    result = apply_gufunc(stats, "(i)->(),()", a, meta=meta)
    expected = stats(a.compute())
    assert_eq(expected[0], result[0])
    assert_eq(expected[1], result[1])


def test_as_gufunc_with_meta():
    stack = da.ones((1, 50, 60), chunks=(1, -1, -1))
    expected = (stack, stack.max())

    meta = (np.array((), dtype=np.float64), np.array((), dtype=np.float64))

    @da.as_gufunc(signature="(i,j) ->(i,j), ()", meta=meta)
    def array_and_max(arr):
        return arr, np.atleast_1d(arr.max())

    result = array_and_max(stack)
    assert_eq(expected[0], result[0])

    # Because `np.max` returns a scalar instead of an `np.ndarray`, we cast
    # the expected output to a `np.ndarray`, as `meta` defines the output
    # should be.
    assert_eq(np.array([expected[1].compute()]), result[1].compute())
