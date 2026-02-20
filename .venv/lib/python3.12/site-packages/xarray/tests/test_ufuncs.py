from __future__ import annotations

import pickle
from unittest.mock import patch

import numpy as np
import numpy.typing as npt
import pytest

import xarray as xr
import xarray.ufuncs as xu
from xarray.tests import assert_allclose, assert_array_equal, mock, requires_dask
from xarray.tests import assert_identical as assert_identical_


def assert_identical(a, b):
    assert type(a) is type(b) or float(a) == float(b)
    if isinstance(a, xr.DataArray | xr.Dataset | xr.Variable):
        assert_identical_(a, b)
    else:
        assert_array_equal(a, b)


@pytest.mark.parametrize(
    "a",
    [
        xr.Variable(["x"], [0, 0]),
        xr.DataArray([0, 0], dims="x"),
        xr.Dataset({"y": ("x", [0, 0])}),
    ],
)
def test_unary(a):
    assert_allclose(a + 1, np.cos(a))


def test_binary():
    args: list[int | float | npt.NDArray | xr.Variable | xr.DataArray | xr.Dataset] = [
        0,
        np.zeros(2),
        xr.Variable(["x"], [0, 0]),
        xr.DataArray([0, 0], dims="x"),
        xr.Dataset({"y": ("x", [0, 0])}),
    ]
    for n, t1 in enumerate(args):
        for t2 in args[n:]:
            assert_identical(t2 + 1, np.maximum(t1, t2 + 1))
            assert_identical(t2 + 1, np.maximum(t2, t1 + 1))
            assert_identical(t2 + 1, np.maximum(t1 + 1, t2))
            assert_identical(t2 + 1, np.maximum(t2 + 1, t1))


def test_binary_out():
    args: list[int | float | npt.NDArray | xr.Variable | xr.DataArray | xr.Dataset] = [
        1,
        np.ones(2),
        xr.Variable(["x"], [1, 1]),
        xr.DataArray([1, 1], dims="x"),
        xr.Dataset({"y": ("x", [1, 1])}),
    ]
    for arg in args:
        actual_mantissa, actual_exponent = np.frexp(arg)
        assert_identical(actual_mantissa, 0.5 * arg)
        assert_identical(actual_exponent, arg)


def test_binary_coord_attrs():
    t = xr.Variable("t", np.arange(2, 4), attrs={"units": "s"})
    x = xr.DataArray(t.values**2, coords={"t": t}, attrs={"units": "s^2"})
    y = xr.DataArray(t.values**3, coords={"t": t}, attrs={"units": "s^3"})
    z1 = xr.apply_ufunc(np.add, x, y, keep_attrs=True)
    assert z1.coords["t"].attrs == {"units": "s"}
    z2 = xr.apply_ufunc(np.add, x, y, keep_attrs=False)
    assert z2.coords["t"].attrs == {}
    # Check also that input array's coordinate attributes weren't affected
    assert t.attrs == {"units": "s"}
    assert x.coords["t"].attrs == {"units": "s"}


def test_groupby():
    ds = xr.Dataset({"a": ("x", [0, 0, 0])}, {"c": ("x", [0, 0, 1])})
    ds_grouped = ds.groupby("c")
    group_mean = ds_grouped.mean("x")
    arr_grouped = ds["a"].groupby("c")

    assert_identical(ds, np.maximum(ds_grouped, group_mean))  # type: ignore[call-overload]
    assert_identical(ds, np.maximum(group_mean, ds_grouped))  # type: ignore[call-overload]

    assert_identical(ds, np.maximum(arr_grouped, group_mean))  # type: ignore[call-overload]
    assert_identical(ds, np.maximum(group_mean, arr_grouped))  # type: ignore[call-overload]

    assert_identical(ds, np.maximum(ds_grouped, group_mean["a"]))  # type: ignore[call-overload]
    assert_identical(ds, np.maximum(group_mean["a"], ds_grouped))  # type: ignore[call-overload]

    assert_identical(ds.a, np.maximum(arr_grouped, group_mean.a))  # type: ignore[call-overload]
    assert_identical(ds.a, np.maximum(group_mean.a, arr_grouped))  # type: ignore[call-overload]

    with pytest.raises(ValueError, match=r"mismatched lengths for dimension"):
        np.maximum(ds.a.variable, ds_grouped)  # type: ignore[call-overload]


def test_alignment():
    ds1 = xr.Dataset({"a": ("x", [1, 2])}, {"x": [0, 1]})
    ds2 = xr.Dataset({"a": ("x", [2, 3]), "b": 4}, {"x": [1, 2]})

    actual = np.add(ds1, ds2)
    expected = xr.Dataset({"a": ("x", [4])}, {"x": [1]})
    assert_identical_(actual, expected)

    with xr.set_options(arithmetic_join="outer"):
        actual = np.add(ds1, ds2)
        expected = xr.Dataset(
            {"a": ("x", [np.nan, 4, np.nan]), "b": np.nan}, coords={"x": [0, 1, 2]}
        )
        assert_identical_(actual, expected)


def test_kwargs():
    x = xr.DataArray(0)
    result = np.add(x, 1, dtype=np.float64)
    assert result.dtype == np.float64


def test_xarray_defers_to_unrecognized_type():
    class Other:
        def __array_ufunc__(self, *args, **kwargs):
            return "other"

    xarray_obj = xr.DataArray([1, 2, 3])
    other = Other()
    assert np.maximum(xarray_obj, other) == "other"  # type: ignore[call-overload]
    assert np.sin(xarray_obj, out=other) == "other"  # type: ignore[call-overload]


def test_xarray_handles_dask():
    da = pytest.importorskip("dask.array")
    x = xr.DataArray(np.ones((2, 2)), dims=["x", "y"])
    y = da.ones((2, 2), chunks=(2, 2))
    result = np.add(x, y)
    assert result.chunks == ((2,), (2,))
    assert isinstance(result, xr.DataArray)


def test_dask_defers_to_xarray():
    da = pytest.importorskip("dask.array")
    x = xr.DataArray(np.ones((2, 2)), dims=["x", "y"])
    y = da.ones((2, 2), chunks=(2, 2))
    result = np.add(y, x)
    assert result.chunks == ((2,), (2,))
    assert isinstance(result, xr.DataArray)


def test_gufunc_methods():
    xarray_obj = xr.DataArray([1, 2, 3])
    with pytest.raises(NotImplementedError, match=r"reduce method"):
        np.add.reduce(xarray_obj, 1)


def test_out():
    xarray_obj = xr.DataArray([1, 2, 3])

    # xarray out arguments should raise
    with pytest.raises(NotImplementedError, match=r"`out` argument"):
        np.add(xarray_obj, 1, out=xarray_obj)  # type: ignore[call-overload]

    # but non-xarray should be OK
    other = np.zeros((3,))
    np.add(other, xarray_obj, out=other)
    assert_identical(other, np.array([1, 2, 3]))


def test_gufuncs():
    xarray_obj = xr.DataArray([1, 2, 3])
    fake_gufunc = mock.Mock(signature="(n)->()", autospec=np.sin)
    with pytest.raises(NotImplementedError, match=r"generalized ufuncs"):
        xarray_obj.__array_ufunc__(fake_gufunc, "__call__", xarray_obj)


class DuckArray(np.ndarray):
    # Minimal subclassed duck array with its own self-contained namespace,
    # which implements a few ufuncs
    def __new__(cls, array):
        obj = np.asarray(array).view(cls)
        return obj

    def __array_namespace__(self, *, api_version=None):
        return DuckArray

    @staticmethod
    def sin(x):
        return np.sin(x)

    @staticmethod
    def add(x, y):
        return x + y


class DuckArray2(DuckArray):
    def __array_namespace__(self, *, api_version=None):
        return DuckArray2


class TestXarrayUfuncs:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.x = xr.DataArray([1, 2, 3])
        self.xd = xr.DataArray(DuckArray([1, 2, 3]))
        self.xd2 = xr.DataArray(DuckArray2([1, 2, 3]))
        self.xt = xr.DataArray(np.datetime64("2021-01-01", "ns"))

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    @pytest.mark.parametrize("name", xu.__all__)
    def test_ufuncs(self, name, request):
        xu_func = getattr(xu, name)
        np_func = getattr(np, name, None)
        if np_func is None and np.lib.NumpyVersion(np.__version__) < "2.0.0":
            pytest.skip(f"Ufunc {name} is not available in numpy {np.__version__}.")

        if name == "isnat":
            args = (self.xt,)
        elif hasattr(np_func, "nin") and np_func.nin == 2:  # type: ignore[union-attr]
            args = (self.x, self.x)  # type: ignore[assignment]
        else:
            args = (self.x,)

        expected = np_func(*args)  # type: ignore[misc]
        actual = xu_func(*args)

        if name in ["angle", "iscomplex"]:
            np.testing.assert_equal(expected, actual.values)
        else:
            assert_identical(actual, expected)

    def test_ufunc_pickle(self):
        a = 1.0
        cos_pickled = pickle.loads(pickle.dumps(xu.cos))
        assert_identical(cos_pickled(a), xu.cos(a))

    def test_ufunc_scalar(self):
        actual = xu.sin(1)
        assert isinstance(actual, float)

    def test_ufunc_duck_array_dataarray(self):
        actual = xu.sin(self.xd)
        assert isinstance(actual.data, DuckArray)

    def test_ufunc_duck_array_variable(self):
        actual = xu.sin(self.xd.variable)
        assert isinstance(actual.data, DuckArray)

    def test_ufunc_duck_array_dataset(self):
        ds = xr.Dataset({"a": self.xd})
        actual = xu.sin(ds)
        assert isinstance(actual.a.data, DuckArray)

    @requires_dask
    def test_ufunc_duck_dask(self):
        import dask.array as da

        x = xr.DataArray(da.from_array(DuckArray(np.array([1, 2, 3]))))
        actual = xu.sin(x)
        assert isinstance(actual.data._meta, DuckArray)

    @requires_dask
    @pytest.mark.xfail(reason="dask ufuncs currently dispatch to numpy")
    def test_ufunc_duck_dask_no_array_ufunc(self):
        import dask.array as da

        # dask ufuncs currently only preserve duck arrays that implement __array_ufunc__
        with patch.object(DuckArray, "__array_ufunc__", new=None, create=True):
            x = xr.DataArray(da.from_array(DuckArray(np.array([1, 2, 3]))))
            actual = xu.sin(x)
            assert isinstance(actual.data._meta, DuckArray)

    def test_ufunc_mixed_arrays_compatible(self):
        actual = xu.add(self.xd, self.x)
        assert isinstance(actual.data, DuckArray)

    def test_ufunc_mixed_arrays_incompatible(self):
        with pytest.raises(ValueError, match=r"Mixed array types"):
            xu.add(self.xd, self.xd2)
