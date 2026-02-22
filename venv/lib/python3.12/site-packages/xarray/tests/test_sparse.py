from __future__ import annotations

import math
import pickle
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import xarray as xr
import xarray.ufuncs as xu
from xarray import DataArray, Variable
from xarray.namedarray.pycompat import array_type
from xarray.tests import assert_equal, assert_identical, requires_dask

filterwarnings = pytest.mark.filterwarnings
param = pytest.param
xfail = pytest.mark.xfail

sparse = pytest.importorskip("sparse")
sparse_array_type = array_type("sparse")


def assert_sparse_equal(a, b):
    assert isinstance(a, sparse_array_type)
    assert isinstance(b, sparse_array_type)
    np.testing.assert_equal(a.todense(), b.todense())


def make_ndarray(shape):
    return np.arange(math.prod(shape)).reshape(shape)


def make_sparray(shape):
    return sparse.random(shape, density=0.1, random_state=0)


def make_xrvar(dim_lengths):
    return xr.Variable(
        tuple(dim_lengths.keys()), make_sparray(shape=tuple(dim_lengths.values()))
    )


def make_xrarray(dim_lengths, coords=None, name="test"):
    if coords is None:
        coords = {d: np.arange(n) for d, n in dim_lengths.items()}
    return xr.DataArray(
        make_sparray(shape=tuple(dim_lengths.values())),
        dims=tuple(coords.keys()),
        coords=coords,
        name=name,
    )


class do:
    def __init__(self, meth, *args, **kwargs):
        self.meth = meth
        self.args = args
        self.kwargs = kwargs

    def __call__(self, obj):
        # cannot pass np.sum when using pytest-xdist
        kwargs = self.kwargs.copy()
        if "func" in self.kwargs:
            kwargs["func"] = getattr(np, kwargs["func"])

        return getattr(obj, self.meth)(*self.args, **kwargs)

    def __repr__(self):
        return f"obj.{self.meth}(*{self.args}, **{self.kwargs})"


@pytest.mark.parametrize(
    "prop",
    [
        "chunks",
        "data",
        "dims",
        "dtype",
        "encoding",
        "imag",
        "nbytes",
        "ndim",
        param("values", marks=xfail(reason="Coercion to dense")),
    ],
)
def test_variable_property(prop):
    var = make_xrvar({"x": 10, "y": 5})
    getattr(var, prop)


@pytest.mark.parametrize(
    "func,sparse_output",
    [
        (do("all"), False),
        (do("any"), False),
        (do("astype", dtype=int), True),
        (do("clip", min=0, max=1), True),
        (do("coarsen", windows={"x": 2}, func="sum"), True),
        (do("compute"), True),
        (do("conj"), True),
        (do("copy"), True),
        (do("count"), False),
        (do("get_axis_num", dim="x"), False),
        (do("isel", x=slice(2, 4)), True),
        (do("isnull"), True),
        (do("load"), True),
        (do("mean"), False),
        (do("notnull"), True),
        (do("roll"), True),
        (do("round"), True),
        (do("set_dims", dim=("x", "y", "z")), True),
        (do("stack", dim={"flat": ("x", "y")}), True),
        (do("to_base_variable"), True),
        (do("transpose"), True),
        (do("unstack", dim={"x": {"x1": 5, "x2": 2}}), True),
        (do("broadcast_equals", make_xrvar({"x": 10, "y": 5})), False),
        (do("equals", make_xrvar({"x": 10, "y": 5})), False),
        (do("identical", make_xrvar({"x": 10, "y": 5})), False),
        param(
            do("argmax"),
            True,
            marks=[
                xfail(reason="Missing implementation for np.argmin"),
                filterwarnings("ignore:Behaviour of argmin/argmax"),
            ],
        ),
        param(
            do("argmin"),
            True,
            marks=[
                xfail(reason="Missing implementation for np.argmax"),
                filterwarnings("ignore:Behaviour of argmin/argmax"),
            ],
        ),
        param(
            do("argsort"),
            True,
            marks=xfail(reason="'COO' object has no attribute 'argsort'"),
        ),
        param(
            do(
                "concat",
                variables=[
                    make_xrvar({"x": 10, "y": 5}),
                    make_xrvar({"x": 10, "y": 5}),
                ],
            ),
            True,
        ),
        param(
            do("conjugate"),
            True,
            marks=xfail(reason="'COO' object has no attribute 'conjugate'"),
        ),
        param(
            do("cumprod"),
            True,
            marks=xfail(reason="Missing implementation for np.nancumprod"),
        ),
        param(
            do("cumsum"),
            True,
            marks=xfail(reason="Missing implementation for np.nancumsum"),
        ),
        (do("fillna", 0), True),
        param(
            do("item", (1, 1)),
            False,
            marks=xfail(reason="'COO' object has no attribute 'item'"),
        ),
        param(
            do("median"),
            False,
            marks=xfail(reason="Missing implementation for np.nanmedian"),
        ),
        param(do("max"), False),
        param(do("min"), False),
        param(
            do("no_conflicts", other=make_xrvar({"x": 10, "y": 5})),
            True,
            marks=xfail(reason="mixed sparse-dense operation"),
        ),
        param(
            do("pad", mode="constant", pad_widths={"x": (1, 1)}, fill_value=5),
            True,
            marks=xfail(reason="Missing implementation for np.pad"),
        ),
        (do("prod"), False),
        param(
            do("quantile", q=0.5),
            True,
            marks=xfail(reason="Missing implementation for np.nanpercentile"),
        ),
        param(
            do("rank", dim="x"),
            False,
            marks=xfail(reason="Only implemented for NumPy arrays (via bottleneck)"),
        ),
        param(
            do("reduce", func="sum", dim="x"),
            True,
        ),
        param(
            do("rolling_window", dim="x", window=2, window_dim="x_win"),
            True,
            marks=xfail(reason="Missing implementation for np.pad"),
        ),
        param(
            do("shift", x=2), True, marks=xfail(reason="mixed sparse-dense operation")
        ),
        param(
            do("std"), False, marks=xfail(reason="Missing implementation for np.nanstd")
        ),
        (do("sum"), False),
        param(
            do("var"), False, marks=xfail(reason="Missing implementation for np.nanvar")
        ),
        param(do("to_dict"), False),
        (do("where", cond=make_xrvar({"x": 10, "y": 5}) > 0.5), True),
    ],
    ids=repr,
)
def test_variable_method(func, sparse_output):
    var_s = make_xrvar({"x": 10, "y": 5})
    var_d = xr.Variable(var_s.dims, var_s.data.todense())
    ret_s = func(var_s)
    ret_d = func(var_d)

    # TODO: figure out how to verify the results of each method
    if isinstance(ret_d, xr.Variable) and isinstance(ret_d.data, sparse.SparseArray):
        ret_d = ret_d.copy(data=ret_d.data.todense())

    if sparse_output:
        assert isinstance(ret_s.data, sparse.SparseArray)
        assert np.allclose(ret_s.data.todense(), ret_d.data, equal_nan=True)
    elif func.meth != "to_dict":
        assert np.allclose(ret_s, ret_d)
    else:
        # pop the arrays from the dict
        arr_s, arr_d = ret_s.pop("data"), ret_d.pop("data")

        assert np.allclose(arr_s, arr_d)
        assert ret_s == ret_d


@pytest.mark.parametrize(
    "func,sparse_output",
    [
        (do("squeeze"), True),
        param(do("to_index"), False, marks=xfail(reason="Coercion to dense")),
        param(do("to_index_variable"), False, marks=xfail(reason="Coercion to dense")),
        param(
            do("searchsorted", 0.5),
            True,
            marks=xfail(reason="'COO' object has no attribute 'searchsorted'"),
        ),
    ],
)
def test_1d_variable_method(func, sparse_output):
    var_s = make_xrvar({"x": 10})
    var_d = xr.Variable(var_s.dims, var_s.data.todense())
    ret_s = func(var_s)
    ret_d = func(var_d)

    if sparse_output:
        assert isinstance(ret_s.data, sparse.SparseArray)
        assert np.allclose(ret_s.data.todense(), ret_d.data)
    else:
        assert np.allclose(ret_s, ret_d)


class TestSparseVariable:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.data = sparse.random((4, 6), random_state=0, density=0.5)
        self.var = xr.Variable(("x", "y"), self.data)

    def test_nbytes(self):
        assert self.var.nbytes == self.data.nbytes

    def test_unary_op(self):
        assert_sparse_equal(-self.var.data, -self.data)
        assert_sparse_equal(abs(self.var).data, abs(self.data))
        assert_sparse_equal(self.var.round().data, self.data.round())

    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_univariate_ufunc(self):
        assert_sparse_equal(np.sin(self.data), np.sin(self.var).data)

    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_bivariate_ufunc(self):
        assert_sparse_equal(np.maximum(self.data, 0), np.maximum(self.var, 0).data)
        assert_sparse_equal(np.maximum(self.data, 0), np.maximum(0, self.var).data)

    def test_univariate_xufunc(self):
        assert_sparse_equal(xu.sin(self.var).data, np.sin(self.data))

    def test_bivariate_xufunc(self):
        assert_sparse_equal(xu.multiply(self.var, 0).data, np.multiply(self.data, 0))
        assert_sparse_equal(xu.multiply(0, self.var).data, np.multiply(0, self.data))

    def test_repr(self):
        expected = dedent(
            """\
            <xarray.Variable (x: 4, y: 6)> Size: 288B
            <COO: shape=(4, 6), dtype=float64, nnz=12, fill_value=0.0>"""
        )
        assert expected == repr(self.var)

    def test_pickle(self):
        v1 = self.var
        v2 = pickle.loads(pickle.dumps(v1))
        assert_sparse_equal(v1.data, v2.data)

    def test_missing_values(self):
        a = np.array([0, 1, np.nan, 3])
        s = sparse.COO.from_numpy(a)
        var_s = Variable("x", s)
        assert np.all(var_s.fillna(2).data.todense() == np.arange(4))
        assert np.all(var_s.count() == 3)


@pytest.mark.parametrize(
    "prop",
    [
        "attrs",
        "chunks",
        "coords",
        "data",
        "dims",
        "dtype",
        "encoding",
        "imag",
        "indexes",
        "loc",
        "name",
        "nbytes",
        "ndim",
        "plot",
        "real",
        "shape",
        "size",
        "sizes",
        "str",
        "variable",
    ],
)
def test_dataarray_property(prop):
    arr = make_xrarray({"x": 10, "y": 5})
    getattr(arr, prop)


@pytest.mark.parametrize(
    "func,sparse_output",
    [
        (do("all"), False),
        (do("any"), False),
        (do("assign_attrs", {"foo": "bar"}), True),
        (do("assign_coords", x=make_xrarray({"x": 10}).x + 1), True),
        (do("astype", int), True),
        (do("clip", min=0, max=1), True),
        (do("compute"), True),
        (do("conj"), True),
        (do("copy"), True),
        (do("count"), False),
        (do("diff", "x"), True),
        (do("drop_vars", "x"), True),
        (do("expand_dims", {"z": 2}, axis=2), True),
        (do("get_axis_num", "x"), False),
        (do("get_index", "x"), False),
        (do("identical", make_xrarray({"x": 5, "y": 5})), False),
        (do("integrate", "x"), True),
        (do("isel", {"x": slice(0, 3), "y": slice(2, 4)}), True),
        (do("isnull"), True),
        (do("load"), True),
        (do("mean"), False),
        (do("persist"), True),
        (do("reindex", {"x": [1, 2, 3]}), True),
        (do("rename", "foo"), True),
        (do("reorder_levels"), True),
        (do("reset_coords", drop=True), True),
        (do("reset_index", "x"), True),
        (do("round"), True),
        (do("sel", x=[0, 1, 2]), True),
        (do("shift"), True),
        (do("sortby", "x", ascending=False), True),
        (do("stack", z=["x", "y"]), True),
        (do("transpose"), True),
        # TODO
        # set_index
        # swap_dims
        (do("broadcast_equals", make_xrvar({"x": 10, "y": 5})), False),
        (do("equals", make_xrvar({"x": 10, "y": 5})), False),
        param(
            do("argmax"),
            True,
            marks=[
                xfail(reason="Missing implementation for np.argmax"),
                filterwarnings("ignore:Behaviour of argmin/argmax"),
            ],
        ),
        param(
            do("argmin"),
            True,
            marks=[
                xfail(reason="Missing implementation for np.argmin"),
                filterwarnings("ignore:Behaviour of argmin/argmax"),
            ],
        ),
        param(
            do("argsort"),
            True,
            marks=xfail(reason="'COO' object has no attribute 'argsort'"),
        ),
        param(
            do("bfill", dim="x"),
            False,
            marks=xfail(reason="Missing implementation for np.flip"),
        ),
        (do("combine_first", make_xrarray({"x": 10, "y": 5})), True),
        param(
            do("conjugate"),
            False,
            marks=xfail(reason="'COO' object has no attribute 'conjugate'"),
        ),
        param(
            do("cumprod"),
            True,
            marks=xfail(reason="Missing implementation for np.nancumprod"),
        ),
        param(
            do("cumsum"),
            True,
            marks=xfail(reason="Missing implementation for np.nancumsum"),
        ),
        param(
            do("differentiate", "x"),
            False,
            marks=xfail(reason="Missing implementation for np.gradient"),
        ),
        param(
            do("dot", make_xrarray({"x": 10, "y": 5})),
            True,
            marks=xfail(reason="Missing implementation for np.einsum"),
        ),
        param(do("dropna", "x"), False, marks=xfail(reason="Coercion to dense")),
        param(do("ffill", "x"), False, marks=xfail(reason="Coercion to dense")),
        (do("fillna", 0), True),
        param(
            do("interp", coords={"x": np.arange(10) + 0.5}),
            True,
            marks=xfail(reason="Coercion to dense"),
        ),
        param(
            do(
                "interp_like",
                make_xrarray(
                    {"x": 10, "y": 5},
                    coords={"x": np.arange(10) + 0.5, "y": np.arange(5) + 0.5},
                ),
            ),
            True,
            marks=xfail(reason="Indexing COO with more than one iterable index"),
        ),
        param(do("interpolate_na", "x"), True, marks=xfail(reason="Coercion to dense")),
        param(
            do("isin", [1, 2, 3]),
            False,
            marks=xfail(reason="Missing implementation for np.isin"),
        ),
        param(
            do("item", (1, 1)),
            False,
            marks=xfail(reason="'COO' object has no attribute 'item'"),
        ),
        param(do("max"), False),
        param(do("min"), False),
        param(
            do("median"),
            False,
            marks=xfail(reason="Missing implementation for np.nanmedian"),
        ),
        (do("notnull"), True),
        (do("pipe", func="sum", axis=1), True),
        (do("prod"), False),
        param(
            do("quantile", q=0.5),
            False,
            marks=xfail(reason="Missing implementation for np.nanpercentile"),
        ),
        param(
            do("rank", "x"),
            False,
            marks=xfail(reason="Only implemented for NumPy arrays (via bottleneck)"),
        ),
        param(
            do("reduce", func="sum", dim="x"),
            False,
            marks=xfail(reason="Coercion to dense"),
        ),
        param(
            do(
                "reindex_like",
                make_xrarray(
                    {"x": 10, "y": 5},
                    coords={"x": np.arange(10) + 0.5, "y": np.arange(5) + 0.5},
                ),
            ),
            True,
            marks=xfail(reason="Indexing COO with more than one iterable index"),
        ),
        (do("roll", x=2, roll_coords=True), True),
        param(
            do("sel", x=[0, 1, 2], y=[2, 3]),
            True,
            marks=xfail(reason="Indexing COO with more than one iterable index"),
        ),
        param(
            do("std"), False, marks=xfail(reason="Missing implementation for np.nanstd")
        ),
        (do("sum"), False),
        param(
            do("var"), False, marks=xfail(reason="Missing implementation for np.nanvar")
        ),
        param(
            do("where", make_xrarray({"x": 10, "y": 5}) > 0.5),
            False,
            marks=xfail(reason="Conversion of dense to sparse when using sparse mask"),
        ),
    ],
    ids=repr,
)
def test_dataarray_method(func, sparse_output):
    arr_s = make_xrarray(
        {"x": 10, "y": 5}, coords={"x": np.arange(10), "y": np.arange(5)}
    )
    arr_d = xr.DataArray(arr_s.data.todense(), coords=arr_s.coords, dims=arr_s.dims)
    ret_s = func(arr_s)
    ret_d = func(arr_d)

    if sparse_output:
        assert isinstance(ret_s.data, sparse.SparseArray)
        assert np.allclose(ret_s.data.todense(), ret_d.data, equal_nan=True)
    else:
        assert np.allclose(ret_s, ret_d, equal_nan=True)


@pytest.mark.parametrize(
    "func,sparse_output",
    [
        (do("squeeze"), True),
        param(
            do("searchsorted", [1, 2, 3]),
            False,
            marks=xfail(reason="'COO' object has no attribute 'searchsorted'"),
        ),
    ],
)
def test_datarray_1d_method(func, sparse_output):
    arr_s = make_xrarray({"x": 10}, coords={"x": np.arange(10)})
    arr_d = xr.DataArray(arr_s.data.todense(), coords=arr_s.coords, dims=arr_s.dims)
    ret_s = func(arr_s)
    ret_d = func(arr_d)

    if sparse_output:
        assert isinstance(ret_s.data, sparse.SparseArray)
        assert np.allclose(ret_s.data.todense(), ret_d.data, equal_nan=True)
    else:
        assert np.allclose(ret_s, ret_d, equal_nan=True)


class TestSparseDataArrayAndDataset:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.sp_ar = sparse.random((4, 6), random_state=0, density=0.5)
        self.sp_xr = xr.DataArray(
            self.sp_ar, coords={"x": range(4)}, dims=("x", "y"), name="foo"
        )
        self.ds_ar = self.sp_ar.todense()
        self.ds_xr = xr.DataArray(
            self.ds_ar, coords={"x": range(4)}, dims=("x", "y"), name="foo"
        )

    def test_to_dataset_roundtrip(self):
        x = self.sp_xr
        assert_equal(x, x.to_dataset("x").to_dataarray("x"))

    def test_align(self):
        a1 = xr.DataArray(
            sparse.COO.from_numpy(np.arange(4)),
            dims=["x"],
            coords={"x": ["a", "b", "c", "d"]},
        )
        b1 = xr.DataArray(
            sparse.COO.from_numpy(np.arange(4)),
            dims=["x"],
            coords={"x": ["a", "b", "d", "e"]},
        )
        a2, b2 = xr.align(a1, b1, join="inner")
        assert isinstance(a2.data, sparse.SparseArray)
        assert isinstance(b2.data, sparse.SparseArray)
        assert np.all(a2.coords["x"].data == ["a", "b", "d"])
        assert np.all(b2.coords["x"].data == ["a", "b", "d"])

    @pytest.mark.xfail(
        reason="COO objects currently do not accept more than one "
        "iterable index at a time"
    )
    def test_align_2d(self):
        A1 = xr.DataArray(
            self.sp_ar,
            dims=["x", "y"],
            coords={
                "x": np.arange(self.sp_ar.shape[0]),
                "y": np.arange(self.sp_ar.shape[1]),
            },
        )

        A2 = xr.DataArray(
            self.sp_ar,
            dims=["x", "y"],
            coords={
                "x": np.arange(1, self.sp_ar.shape[0] + 1),
                "y": np.arange(1, self.sp_ar.shape[1] + 1),
            },
        )

        B1, B2 = xr.align(A1, A2, join="inner")
        assert np.all(B1.coords["x"] == np.arange(1, self.sp_ar.shape[0]))
        assert np.all(B1.coords["y"] == np.arange(1, self.sp_ar.shape[0]))
        assert np.all(B1.coords["x"] == B2.coords["x"])
        assert np.all(B1.coords["y"] == B2.coords["y"])

    def test_align_outer(self):
        a1 = xr.DataArray(
            sparse.COO.from_numpy(np.arange(4)),
            dims=["x"],
            coords={"x": ["a", "b", "c", "d"]},
        )
        b1 = xr.DataArray(
            sparse.COO.from_numpy(np.arange(4)),
            dims=["x"],
            coords={"x": ["a", "b", "d", "e"]},
        )
        a2, b2 = xr.align(a1, b1, join="outer")
        assert isinstance(a2.data, sparse.SparseArray)
        assert isinstance(b2.data, sparse.SparseArray)
        assert np.all(a2.coords["x"].data == ["a", "b", "c", "d", "e"])
        assert np.all(b2.coords["x"].data == ["a", "b", "c", "d", "e"])

    def test_concat(self):
        ds1 = xr.Dataset(data_vars={"d": self.sp_xr})
        ds2 = xr.Dataset(data_vars={"d": self.sp_xr})
        ds3 = xr.Dataset(data_vars={"d": self.sp_xr})
        out = xr.concat([ds1, ds2, ds3], dim="x")
        assert_sparse_equal(
            out["d"].data,
            sparse.concatenate([self.sp_ar, self.sp_ar, self.sp_ar], axis=0),
        )

        out_concat = xr.concat([self.sp_xr, self.sp_xr, self.sp_xr], dim="y")
        assert_sparse_equal(
            out_concat.data,
            sparse.concatenate([self.sp_ar, self.sp_ar, self.sp_ar], axis=1),
        )

    def test_stack(self):
        arr = make_xrarray({"w": 2, "x": 3, "y": 4})
        stacked = arr.stack(z=("x", "y"))

        z = pd.MultiIndex.from_product(
            [list(range(3)), list(range(4))], names=["x", "y"]
        )

        expected = xr.DataArray(
            arr.data.reshape((2, -1)), {"w": [0, 1], "z": z}, dims=["w", "z"]
        )

        assert_equal(expected, stacked)

        roundtripped = stacked.unstack()
        assert_identical(arr, roundtripped)

    def test_dataarray_repr(self):
        a = xr.DataArray(
            sparse.COO.from_numpy(np.ones(4)),
            dims=["x"],
            coords={"y": ("x", sparse.COO.from_numpy(np.arange(4, dtype="i8")))},
        )
        expected = dedent(
            """\
            <xarray.DataArray (x: 4)> Size: 64B
            <COO: shape=(4,), dtype=float64, nnz=4, fill_value=0.0>
            Coordinates:
                y        (x) int64 48B <COO: nnz=3, fill_value=0>
            Dimensions without coordinates: x"""
        )
        assert expected == repr(a)

    def test_dataset_repr(self):
        ds = xr.Dataset(
            data_vars={"a": ("x", sparse.COO.from_numpy(np.ones(4)))},
            coords={"y": ("x", sparse.COO.from_numpy(np.arange(4, dtype="i8")))},
        )
        expected = dedent(
            """\
            <xarray.Dataset> Size: 112B
            Dimensions:  (x: 4)
            Coordinates:
                y        (x) int64 48B <COO: nnz=3, fill_value=0>
            Dimensions without coordinates: x
            Data variables:
                a        (x) float64 64B <COO: nnz=4, fill_value=0.0>"""
        )
        assert expected == repr(ds)

    @requires_dask
    def test_sparse_dask_dataset_repr(self):
        ds = xr.Dataset(
            data_vars={"a": ("x", sparse.COO.from_numpy(np.ones(4)))}
        ).chunk()
        if Version(sparse.__version__) >= Version("0.16.0"):
            meta = "sparse.numba_backend._coo.core.COO"
        else:
            meta = "sparse.COO"
        expected = dedent(
            f"""\
            <xarray.Dataset> Size: 32B
            Dimensions:  (x: 4)
            Dimensions without coordinates: x
            Data variables:
                a        (x) float64 32B dask.array<chunksize=(4,), meta={meta}>"""
        )
        assert expected == repr(ds)

    def test_dataarray_pickle(self):
        a1 = xr.DataArray(
            sparse.COO.from_numpy(np.ones(4)),
            dims=["x"],
            coords={"y": ("x", sparse.COO.from_numpy(np.arange(4)))},
        )
        a2 = pickle.loads(pickle.dumps(a1))
        assert_identical(a1, a2)

    def test_dataset_pickle(self):
        ds1 = xr.Dataset(
            data_vars={"a": ("x", sparse.COO.from_numpy(np.ones(4)))},
            coords={"y": ("x", sparse.COO.from_numpy(np.arange(4)))},
        )
        ds2 = pickle.loads(pickle.dumps(ds1))
        assert_identical(ds1, ds2)

    def test_coarsen(self):
        a1 = self.ds_xr
        a2 = self.sp_xr
        m1 = a1.coarsen(x=2, boundary="trim").mean()  # type: ignore[attr-defined]
        m2 = a2.coarsen(x=2, boundary="trim").mean()  # type: ignore[attr-defined]

        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail(reason="No implementation of np.pad")
    def test_rolling(self):
        a1 = self.ds_xr
        a2 = self.sp_xr
        m1 = a1.rolling(x=2, center=True).mean()
        m2 = a2.rolling(x=2, center=True).mean()

        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail(reason="Coercion to dense")
    def test_rolling_exp(self):
        a1 = self.ds_xr
        a2 = self.sp_xr
        m1 = a1.rolling_exp(x=2, center=True).mean()
        m2 = a2.rolling_exp(x=2, center=True).mean()

        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail(reason="No implementation of np.einsum")
    def test_dot(self):
        a1 = self.sp_xr.dot(self.sp_xr[0])
        a2 = self.sp_ar.dot(self.sp_ar[0])
        assert_equal(a1, a2)

    @pytest.mark.xfail(reason="Groupby reductions produce dense output")
    def test_groupby(self):
        x1 = self.ds_xr
        x2 = self.sp_xr
        m1 = x1.groupby("x").mean(...)
        m2 = x2.groupby("x").mean(...)
        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail(reason="Groupby reductions produce dense output")
    def test_groupby_first(self):
        x = self.sp_xr.copy()
        x.coords["ab"] = ("x", ["a", "a", "b", "b"])
        x.groupby("ab").first()
        x.groupby("ab").first(skipna=False)

    @pytest.mark.xfail(reason="Groupby reductions produce dense output")
    def test_groupby_bins(self):
        x1 = self.ds_xr
        x2 = self.sp_xr
        m1 = x1.groupby_bins("x", bins=[0, 3, 7, 10]).sum(...)
        m2 = x2.groupby_bins("x", bins=[0, 3, 7, 10]).sum(...)
        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail(reason="Resample produces dense output")
    def test_resample(self):
        t1 = xr.DataArray(
            np.linspace(0, 11, num=12),
            coords=[
                pd.date_range("1999-12-15", periods=12, freq=pd.DateOffset(months=1))
            ],
            dims="time",
        )
        t2 = t1.copy()
        t2.data = sparse.COO(t2.data)
        m1 = t1.resample(time="QS-DEC").mean()
        m2 = t2.resample(time="QS-DEC").mean()
        assert isinstance(m2.data, sparse.SparseArray)
        assert np.allclose(m1.data, m2.data.todense())

    @pytest.mark.xfail
    def test_reindex(self):
        x1 = self.ds_xr
        x2 = self.sp_xr
        for kwargs in [
            {"x": [2, 3, 4]},
            {"x": [1, 100, 2, 101, 3]},
            {"x": [2.5, 3, 3.5], "y": [2, 2.5, 3]},
        ]:
            m1 = x1.reindex(**kwargs)  # type: ignore[arg-type]
            m2 = x2.reindex(**kwargs)  # type: ignore[arg-type]
            assert np.allclose(m1, m2, equal_nan=True)

    @pytest.mark.xfail
    def test_merge(self):
        x = self.sp_xr
        y = xr.merge([x, x.rename("bar")]).to_dataarray()
        assert isinstance(y, sparse.SparseArray)

    @pytest.mark.xfail
    def test_where(self):
        a = np.arange(10)
        cond = a > 3
        xr.DataArray(a).where(cond)

        s = sparse.COO.from_numpy(a)
        cond2 = s > 3
        xr.DataArray(s).where(cond2)

        x = xr.DataArray(s)
        cond3: DataArray = x > 3
        x.where(cond3)


class TestSparseCoords:
    @pytest.mark.xfail(reason="Coercion of coords to dense")
    def test_sparse_coords(self):
        xr.DataArray(
            sparse.COO.from_numpy(np.arange(4)),
            dims=["x"],
            coords={"x": sparse.COO.from_numpy([1, 2, 3, 4])},
        )


@requires_dask
def test_chunk():
    s = sparse.COO.from_numpy(np.array([0, 0, 1, 2]))
    a = DataArray(s)
    ac = a.chunk(2)
    assert ac.chunks == ((2, 2),)
    assert isinstance(ac.data._meta, sparse.COO)
    assert_identical(ac, a)

    ds = a.to_dataset(name="a")
    dsc = ds.chunk(2)
    assert dsc.chunks == {"dim_0": (2, 2)}
    assert_identical(dsc, ds)


@requires_dask
def test_dask_token():
    import dask

    s = sparse.COO.from_numpy(np.array([0, 0, 1, 2]))
    a = DataArray(s)
    t1 = dask.base.tokenize(a)
    t2 = dask.base.tokenize(a)
    t3 = dask.base.tokenize(a + 1)
    assert t1 == t2
    assert t3 != t2
    assert isinstance(a.data, sparse.COO)

    ac = a.chunk(2)
    t4 = dask.base.tokenize(ac)
    t5 = dask.base.tokenize(ac + 1)
    assert t4 != t5
    assert isinstance(ac.data._meta, sparse.COO)


@requires_dask
def test_apply_ufunc_check_meta_coherence():
    s = sparse.COO.from_numpy(np.array([0, 0, 1, 2]))
    a = DataArray(s)
    ac = a.chunk(2)
    sparse_meta = ac.data._meta

    result = xr.apply_ufunc(lambda x: x, ac, dask="parallelized").data._meta

    assert_sparse_equal(result, sparse_meta)
