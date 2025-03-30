from __future__ import annotations

import pickle
import warnings
from functools import partial
from operator import add

import pytest

np = pytest.importorskip("numpy")

import dask.array as da
from dask.array.ufunc import da_frompyfunc
from dask.array.utils import assert_eq
from dask.base import tokenize

DISCLAIMER = """
This docstring was copied from numpy.{name}.

Some inconsistencies with the Dask version may exist.
"""


@pytest.mark.parametrize("name", ["log", "modf", "frexp"])
def test_ufunc_meta(name):
    disclaimer = DISCLAIMER.format(name=name)
    skip_test = "  # doctest: +SKIP"
    ufunc = getattr(da, name)
    assert ufunc.__name__ == name
    assert disclaimer in ufunc.__doc__

    assert (
        ufunc.__doc__.replace(disclaimer, "").replace(skip_test, "")
        == getattr(np, name).__doc__
    )


def test_ufunc():
    for attr in ["nin", "nargs", "nout", "ntypes", "identity", "signature", "types"]:
        assert getattr(da.log, attr) == getattr(np.log, attr)

    with pytest.raises(AttributeError):
        da.log.not_an_attribute

    assert repr(da.log) == repr(np.log)
    assert "nin" in dir(da.log)
    assert "outer" in dir(da.log)


binary_ufuncs = [
    "add",
    "arctan2",
    "copysign",
    "divide",
    "equal",
    "bitwise_and",
    "bitwise_or",
    "bitwise_xor",
    "floor_divide",
    "fmax",
    "fmin",
    "fmod",
    "greater",
    "greater_equal",
    "hypot",
    "ldexp",
    "left_shift",
    "less",
    "less_equal",
    "logaddexp",
    "logaddexp2",
    "logical_and",
    "logical_or",
    "logical_xor",
    "maximum",
    "minimum",
    "mod",
    "multiply",
    "nextafter",
    "not_equal",
    "power",
    "remainder",
    "right_shift",
    "subtract",
    "true_divide",
    "float_power",
]

unary_ufuncs = [
    "abs",
    "absolute",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctanh",
    "bitwise_not",
    "cbrt",
    "ceil",
    "conj",
    "cos",
    "cosh",
    "deg2rad",
    "degrees",
    "exp",
    "exp2",
    "expm1",
    "fabs",
    "fix",
    "floor",
    "invert",
    "isfinite",
    "isinf",
    "isnan",
    "log",
    "log10",
    "log1p",
    "log2",
    "logical_not",
    "negative",
    "positive",
    "rad2deg",
    "radians",
    "reciprocal",
    "rint",
    "sign",
    "signbit",
    "sin",
    "sinh",
    "spacing",
    "sqrt",
    "square",
    "tan",
    "tanh",
    "trunc",
]


@pytest.mark.parametrize("ufunc", unary_ufuncs)
def test_unary_ufunc(ufunc):
    if ufunc == "fix":
        pytest.skip("fix calls floor in a way that we do not yet support")
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr = np.random.randint(1, 100, size=(20, 20))
    darr = da.from_array(arr, 3)

    with warnings.catch_warnings():  # some invalid values (arccos, arcsin, etc.)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # applying Dask ufunc doesn't trigger computation
        assert isinstance(dafunc(darr), da.Array)
        assert_eq(dafunc(darr), npfunc(arr), equal_nan=True)

    with warnings.catch_warnings():  # some invalid values (arccos, arcsin, etc.)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # applying NumPy ufunc is lazy
        if isinstance(npfunc, np.ufunc):
            assert isinstance(npfunc(darr), da.Array)
        else:
            assert isinstance(npfunc(darr), np.ndarray)
        assert_eq(npfunc(darr), npfunc(arr), equal_nan=True)

    with warnings.catch_warnings():  # some invalid values (arccos, arcsin, etc.)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # applying Dask ufunc to normal ndarray triggers computation
        assert isinstance(dafunc(arr), np.ndarray)
        assert_eq(dafunc(arr), npfunc(arr), equal_nan=True)


@pytest.mark.parametrize("ufunc", binary_ufuncs)
def test_binary_ufunc(ufunc):
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr1 = np.random.randint(1, 100, size=(20, 20))
    darr1 = da.from_array(arr1, 3)

    arr2 = np.random.randint(1, 100, size=(20, 20))
    darr2 = da.from_array(arr2, 3)

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(darr1, darr2), da.Array)
    assert_eq(dafunc(darr1, darr2), npfunc(arr1, arr2))

    # applying NumPy ufunc triggers computation or is lazy
    assert isinstance(npfunc(darr1, darr2), da.Array)
    assert_eq(npfunc(darr1, darr2), npfunc(arr1, arr2))

    # applying Dask ufunc to normal ndarray triggers computation
    assert isinstance(dafunc(arr1, arr2), np.ndarray)
    assert_eq(dafunc(arr1, arr2), npfunc(arr1, arr2))

    # with scalar
    assert isinstance(dafunc(darr1, 10), da.Array)
    assert_eq(dafunc(darr1, 10), npfunc(arr1, 10))

    with warnings.catch_warnings():  # overflow in ldexp
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        assert isinstance(dafunc(10, darr1), da.Array)
        assert_eq(dafunc(10, darr1), npfunc(10, arr1))

    assert isinstance(dafunc(arr1, 10), np.ndarray)
    assert_eq(dafunc(arr1, 10), npfunc(arr1, 10))

    with warnings.catch_warnings():  # overflow in ldexp
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        assert isinstance(dafunc(10, arr1), np.ndarray)
        assert_eq(dafunc(10, arr1), npfunc(10, arr1))


def test_ufunc_outer():
    arr1 = np.random.randint(1, 100, size=20)
    darr1 = da.from_array(arr1, 3)

    arr2 = np.random.randint(1, 100, size=(10, 3))
    darr2 = da.from_array(arr2, 3)

    # Check output types
    assert isinstance(da.add.outer(darr1, darr2), da.Array)
    assert isinstance(da.add.outer(arr1, darr2), da.Array)
    assert isinstance(da.add.outer(darr1, arr2), da.Array)
    assert isinstance(da.add.outer(arr1, arr2), np.ndarray)

    # Check mix of dimensions, dtypes, and numpy/dask/object
    cases = [
        ((darr1, darr2), (arr1, arr2)),
        ((darr2, darr1), (arr2, arr1)),
        ((darr2, darr1.astype("f8")), (arr2, arr1.astype("f8"))),
        ((darr1, arr2), (arr1, arr2)),
        ((darr1, 1), (arr1, 1)),
        ((1, darr2), (1, arr2)),
        ((1.5, darr2), (1.5, arr2)),
        (([1, 2, 3], darr2), ([1, 2, 3], arr2)),
        ((darr1.sum(), darr2), (arr1.sum(), arr2)),
        ((np.array(1), darr2), (np.array(1), arr2)),
    ]

    for (dA, dB), (A, B) in cases:
        assert_eq(da.add.outer(dA, dB), np.add.outer(A, B))

    # Check dtype kwarg works
    assert_eq(
        da.add.outer(darr1, darr2, dtype="f8"), np.add.outer(arr1, arr2, dtype="f8")
    )

    with pytest.raises(ValueError):
        da.add.outer(darr1, darr2, out=arr1)

    with pytest.raises(ValueError):
        da.sin.outer(darr1, darr2)


@pytest.mark.parametrize("ufunc", ["isreal", "iscomplex", "real", "imag"])
def test_complex(ufunc):
    dafunc = getattr(da, ufunc)
    # Note that these functions are not NumPy ufuncs
    npfunc = getattr(np, ufunc)

    real = np.random.randint(1, 100, size=(20, 20))
    imag = np.random.randint(1, 100, size=(20, 20)) * 1j
    comp = real + imag

    dareal = da.from_array(real, 3)
    daimag = da.from_array(imag, 3)
    dacomp = da.from_array(comp, 3)

    assert_eq(dacomp.real, comp.real)
    assert_eq(dacomp.imag, comp.imag)
    assert_eq(dacomp.conj(), comp.conj())

    for darr, arr in [(dacomp, comp), (dareal, real), (daimag, imag)]:
        # applying Dask ufunc doesn't trigger computation
        assert isinstance(dafunc(darr), da.Array)
        assert_eq(dafunc(darr), npfunc(arr))
        assert_eq(npfunc(darr), npfunc(arr))

        # applying Dask ufunc to normal ndarray triggers computation
        assert isinstance(dafunc(arr), np.ndarray)
        assert_eq(dafunc(arr), npfunc(arr))


@pytest.mark.parametrize("ufunc", ["frexp", "modf"])
def test_ufunc_2results(ufunc):
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr = np.random.randint(1, 100, size=(20, 20))
    darr = da.from_array(arr, 3)

    # applying Dask ufunc doesn't trigger computation
    res1, res2 = dafunc(darr)
    assert isinstance(res1, da.Array)
    assert isinstance(res2, da.Array)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)

    # applying NumPy ufunc is now lazy
    res1, res2 = npfunc(darr)
    assert isinstance(res1, da.Array)
    assert isinstance(res2, da.Array)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)

    # applying Dask ufunc to normal ndarray triggers computation
    res1, res2 = dafunc(arr)
    assert isinstance(res1, da.Array)
    assert isinstance(res2, da.Array)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)


def test_clip():
    x = np.random.normal(0, 10, size=(10, 10))
    d = da.from_array(x, chunks=(3, 4))

    assert_eq(x.clip(5), d.clip(5))
    assert_eq(x.clip(1, 5), d.clip(1, 5))
    assert_eq(x.clip(min=5), d.clip(min=5))
    assert_eq(x.clip(max=5), d.clip(max=5))
    assert_eq(x.clip(max=1, min=5), d.clip(max=1, min=5))
    assert_eq(x.clip(min=1, max=5), d.clip(min=1, max=5))


def test_angle():
    real = np.random.randint(1, 100, size=(20, 20))
    imag = np.random.randint(1, 100, size=(20, 20)) * 1j
    comp = real + imag
    dacomp = da.from_array(comp, 3)

    assert_eq(da.angle(dacomp), np.angle(comp))
    assert_eq(da.angle(dacomp, deg=True), np.angle(comp, deg=True))
    assert isinstance(da.angle(comp), np.ndarray)
    assert_eq(da.angle(comp), np.angle(comp))


def test_issignedinf():
    with np.errstate(invalid="ignore", divide="ignore"):
        arr = np.random.randint(-1, 2, size=(20, 20)).astype(float) / 0
    darr = da.from_array(arr, 3)

    assert_eq(np.isneginf(arr), da.isneginf(darr))
    assert_eq(np.isposinf(arr), da.isposinf(darr))


@pytest.mark.parametrize("func", ["i0", "sinc", "nan_to_num"])
def test_non_ufunc_others(func):
    arr = np.random.randint(1, 100, size=(20, 20))
    darr = da.from_array(arr, 3)

    dafunc = getattr(da, func)
    npfunc = getattr(np, func)

    assert_eq(dafunc(darr), npfunc(arr), equal_nan=True)


def test_frompyfunc():
    myadd = da.frompyfunc(add, 2, 1)
    np_myadd = np.frompyfunc(add, 2, 1)

    x = np.random.normal(0, 10, size=(10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    y = np.random.normal(0, 10, size=10)
    dy = da.from_array(y, chunks=2)

    assert_eq(myadd(dx, dy), np_myadd(x, y))
    assert_eq(myadd.outer(dx, dy), np_myadd.outer(x, y))

    with pytest.raises(NotImplementedError):
        da.frompyfunc(lambda x, y: (x + y, x - y), 2, 2)


def test_frompyfunc_wrapper():
    f = da_frompyfunc(add, 2, 1)
    np_f = np.frompyfunc(add, 2, 1)
    x = np.array([1, 2, 3])

    # Callable
    np.testing.assert_equal(f(x, 1), np_f(x, 1))

    # picklable
    f2 = pickle.loads(pickle.dumps(f))
    np.testing.assert_equal(f2(x, 1), np_f(x, 1))

    # Attributes
    assert f.ntypes == np_f.ntypes
    with pytest.raises(AttributeError):
        f.not_an_attribute

    # Tab completion
    assert "ntypes" in dir(f)

    # Methods
    np.testing.assert_equal(f.outer(x, x), np_f.outer(x, x))

    # funcname
    assert f.__name__ == "frompyfunc-add"

    # repr
    assert repr(f) == "da.frompyfunc<add, 2, 1>"

    # tokenize
    assert tokenize(da_frompyfunc(add, 2, 1)) == tokenize(da_frompyfunc(add, 2, 1))


def test_array_ufunc():
    x = np.arange(24).reshape((4, 6))
    d = da.from_array(x, chunks=(2, 3))

    for func in [np.sin, np.sum, np.negative, partial(np.prod, axis=0)]:
        assert isinstance(func(d), da.Array)
        assert_eq(func(d), func(x))


def test_array_ufunc_binop():
    x = np.arange(25).reshape((5, 5))
    d = da.from_array(x, chunks=(2, 2))

    for func in [np.add, np.multiply]:
        assert isinstance(func(d, d), da.Array)
        assert_eq(func(d, d), func(x, x))

        assert isinstance(func.outer(d, d), da.Array)
        assert_eq(func.outer(d, d), func.outer(x, x))


def test_array_ufunc_out():
    x = da.arange(10, chunks=(5,))
    np.sin(x, out=x)
    np.add(x, 10, out=x)
    assert_eq(x, np.sin(np.arange(10)) + 10)


def test_unsupported_ufunc_methods():
    x = da.arange(10, chunks=(5,))
    with pytest.raises(TypeError):
        assert np.add.reduce(x)


def test_out_numpy():
    x = da.arange(10, chunks=(5,))
    empty = np.empty(10, dtype=x.dtype)
    with pytest.raises((TypeError, NotImplementedError)) as info:
        np.add(x, 1, out=empty)

    assert "ndarray" in str(info.value)
    assert "Array" in str(info.value)


def test_out_shape_mismatch():
    x = da.arange(10, chunks=(5,))
    y = da.arange(15, chunks=(5,))
    with pytest.raises(ValueError):
        assert np.log(x, out=y)


def test_divmod():
    arr1 = np.random.randint(1, 100, size=(20, 20))
    arr2 = np.random.randint(1, 100, size=(20, 20))

    darr1 = da.from_array(arr1, 3)
    darr2 = da.from_array(arr2, 3)

    result = np.divmod(darr1, 2.0)
    expected = np.divmod(arr1, 2.0)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])

    result = np.divmod(darr1, darr2)
    expected = np.divmod(arr1, arr2)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])

    result = divmod(darr1, 2.0)
    expected = divmod(arr1, 2.0)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])

    result = divmod(darr1, darr2)
    expected = divmod(arr1, arr2)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])


@pytest.mark.parametrize("dt", ["float64", "float32", "int32", "int64"])
def test_dtype_kwarg(dt):
    arr1 = np.array([1, 2, 3])
    arr2 = np.array([4, 5, 6])

    darr1 = da.from_array(arr1)
    darr2 = da.from_array(arr2)

    expected = np.add(arr1, arr2, dtype=dt)
    result = np.add(darr1, darr2, dtype=dt)
    assert_eq(expected, result)

    result = da.add(darr1, darr2, dtype=dt)
    assert_eq(expected, result)


@pytest.mark.parametrize("dtype", [None, "f8"])
@pytest.mark.parametrize("left_is_da", [False, True])
@pytest.mark.parametrize("right_is_da", [False, True])
@pytest.mark.parametrize("where_kind", [True, False, "numpy", "dask"])
def test_ufunc_where(dtype, left_is_da, right_is_da, where_kind):
    left = np.arange(12).reshape((3, 4))
    right = np.arange(4)
    out = np.zeros_like(left, dtype=dtype)
    d_out = da.zeros_like(left, dtype=dtype)

    if where_kind in (True, False):
        d_where = where = where_kind
    else:
        d_where = where = np.array([False, True, True, False])
        if where_kind == "dask":
            d_where = da.from_array(where, chunks=2)

    d_left = da.from_array(left, chunks=2) if left_is_da else left
    d_right = da.from_array(right, chunks=2) if right_is_da else right

    expected = np.add(left, right, where=where, out=out, dtype=dtype)
    result = da.add(d_left, d_right, where=d_where, out=d_out, dtype=dtype)
    assert result is d_out
    assert_eq(expected, result)


@pytest.mark.parametrize("left_is_da", [False, True])
@pytest.mark.parametrize("right_is_da", [False, True])
@pytest.mark.parametrize("where_is_da", [False, True])
def test_ufunc_where_broadcasts(left_is_da, right_is_da, where_is_da):
    left = np.arange(4)
    right = np.arange(4, 8)
    where = np.array([[0, 1, 1, 0], [1, 0, 0, 1], [0, 1, 0, 1]]).astype("bool")
    out = np.zeros(where.shape, dtype=left.dtype)

    d_out = da.zeros(where.shape, dtype=left.dtype)
    d_where = da.from_array(where, chunks=2) if where_is_da else where
    d_left = da.from_array(left, chunks=2) if left_is_da else left
    d_right = da.from_array(right, chunks=2) if right_is_da else right

    expected = np.add(left, right, where=where, out=out)
    result = da.add(d_left, d_right, where=d_where, out=d_out)
    assert result is d_out
    assert_eq(expected, result)


def test_ufunc_where_no_out():
    left = np.arange(4)
    right = np.arange(4, 8)
    where = np.array([[0, 1, 1, 0], [1, 0, 0, 1], [0, 1, 0, 1]]).astype("bool")

    d_where = da.from_array(where, chunks=2)
    d_left = da.from_array(left, chunks=2)
    d_right = da.from_array(right, chunks=2)

    expected = np.add(left, right, where=where)
    result = da.add(d_left, d_right, where=d_where)

    # If no `out` is provided, numpy leaves elements that don't match `where`
    # uninitialized, so they effectively may be any random value.  We test that
    # the set values match, and that the unset values aren't equal to if
    # `where` wasn't provided (to test that `where` was actually passed).

    expected_masked = np.where(where, expected, 0)
    result_masked = np.where(where, expected, 0)
    assert_eq(expected_masked, result_masked)

    expected_no_where = np.add(left, right)
    assert not np.equal(result.compute(), expected_no_where).all()


def test_ufunc_where_doesnt_mutate_out():
    """Dask array's are immutable, ensure that the backing numpy array for
    `out` isn't actually mutated"""
    left = da.from_array(np.arange(4, dtype="i8"), chunks=2)
    right = da.from_array(np.arange(4, 8, dtype="i8"), chunks=2)
    where = da.from_array(np.array([1, 0, 0, 1], dtype="bool"), chunks=2)
    out_np = np.zeros(4, dtype="i8")
    out = da.from_array(out_np, chunks=2)
    result = da.add(left, right, where=where, out=out)
    assert out is result
    assert_eq(out, np.array([4, 0, 0, 10], dtype="i8"))

    # Check that original `out` array isn't mutated
    assert np.equal(out_np, 0).all()
