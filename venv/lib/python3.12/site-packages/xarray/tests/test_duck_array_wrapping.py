import numpy as np
import pandas as pd
import pytest

import xarray as xr

# Don't run cupy in CI because it requires a GPU
NAMESPACE_ARRAYS = {
    "cupy": {
        "attrs": {
            "array": "ndarray",
            "constructor": "asarray",
        },
        "xfails": {"quantile": "no nanquantile"},
    },
    "dask.array": {
        "attrs": {
            "array": "Array",
            "constructor": "from_array",
        },
        "xfails": {
            "argsort": "no argsort",
            "conjugate": "conj but no conjugate",
            "searchsorted": "dask.array.searchsorted but no Array.searchsorted",
        },
    },
    "jax.numpy": {
        "attrs": {
            "array": "ndarray",
            "constructor": "asarray",
        },
        "xfails": {
            "rolling_construct": "no sliding_window_view",
            "rolling_reduce": "no sliding_window_view",
            "cumulative_construct": "no sliding_window_view",
            "cumulative_reduce": "no sliding_window_view",
        },
    },
    "pint": {
        "attrs": {
            "array": "Quantity",
            "constructor": "Quantity",
        },
        "xfails": {
            "all": "returns a bool",
            "any": "returns a bool",
            "argmax": "returns an int",
            "argmin": "returns an int",
            "argsort": "returns an int",
            "count": "returns an int",
            "dot": "no tensordot",
            "full_like": "should work, see: https://github.com/hgrecco/pint/pull/1669",
            "idxmax": "returns the coordinate",
            "idxmin": "returns the coordinate",
            "isin": "returns a bool",
            "isnull": "returns a bool",
            "notnull": "returns a bool",
            "rolling_reduce": "no dispatch for numbagg/bottleneck",
            "cumulative_reduce": "no dispatch for numbagg/bottleneck",
            "searchsorted": "returns an int",
            "weighted": "no tensordot",
        },
    },
    "sparse": {
        "attrs": {
            "array": "COO",
            "constructor": "COO",
        },
        "xfails": {
            "cov": "dense output",
            "corr": "no nanstd",
            "cross": "no cross",
            "count": "dense output",
            "dot": "fails on some platforms/versions",
            "isin": "no isin",
            "rolling_construct": "no sliding_window_view",
            "rolling_reduce": "no sliding_window_view",
            "cumulative_construct": "no sliding_window_view",
            "cumulative_reduce": "no sliding_window_view",
            "coarsen_construct": "pad constant_values must be fill_value",
            "coarsen_reduce": "pad constant_values must be fill_value",
            "weighted": "fill_value error",
            "coarsen": "pad constant_values must be fill_value",
            "quantile": "no non skipping version",
            "differentiate": "no gradient",
            "argmax": "no nan skipping version",
            "argmin": "no nan skipping version",
            "idxmax": "no nan skipping version",
            "idxmin": "no nan skipping version",
            "median": "no nan skipping version",
            "std": "no nan skipping version",
            "var": "no nan skipping version",
            "cumsum": "no cumsum",
            "cumprod": "no cumprod",
            "argsort": "no argsort",
            "conjugate": "no conjugate",
            "searchsorted": "no searchsorted",
            "shift": "pad constant_values must be fill_value",
            "pad": "pad constant_values must be fill_value",
        },
    },
}

try:
    import jax  # type: ignore[import-not-found,unused-ignore]

    # enable double-precision
    jax.config.update("jax_enable_x64", True)
except ImportError:
    pass


class _BaseTest:
    def setup_for_test(self, request, namespace):
        self.namespace = namespace
        self.xp = pytest.importorskip(namespace)
        self.Array = getattr(self.xp, NAMESPACE_ARRAYS[namespace]["attrs"]["array"])
        self.constructor = getattr(
            self.xp, NAMESPACE_ARRAYS[namespace]["attrs"]["constructor"]
        )
        xarray_method = request.node.name.split("test_")[1].split("[")[0]
        if xarray_method in NAMESPACE_ARRAYS[namespace]["xfails"]:
            reason = NAMESPACE_ARRAYS[namespace]["xfails"][xarray_method]
            pytest.xfail(f"xfail for {self.namespace}: {reason}")

    def get_test_dataarray(self):
        data = np.asarray([[1, 2, 3, np.nan, 5]])
        x = np.arange(5)
        data = self.constructor(data)
        return xr.DataArray(
            data,
            dims=["y", "x"],
            coords={"y": [1], "x": x},
            name="foo",
        )


@pytest.mark.parametrize("namespace", NAMESPACE_ARRAYS)
class TestTopLevelMethods(_BaseTest):
    @pytest.fixture(autouse=True)
    def setUp(self, request, namespace):
        self.setup_for_test(request, namespace)
        self.x1 = self.get_test_dataarray()
        self.x2 = self.get_test_dataarray().assign_coords(x=np.arange(2, 7))

    def test_apply_ufunc(self):
        func = lambda x: x + 1
        result = xr.apply_ufunc(func, self.x1, dask="parallelized")
        assert isinstance(result.data, self.Array)

    def test_align(self):
        result = xr.align(self.x1, self.x2)
        assert isinstance(result[0].data, self.Array)
        assert isinstance(result[1].data, self.Array)

    def test_broadcast(self):
        result = xr.broadcast(self.x1, self.x2)
        assert isinstance(result[0].data, self.Array)
        assert isinstance(result[1].data, self.Array)

    def test_concat(self):
        result = xr.concat([self.x1, self.x2], dim="x")
        assert isinstance(result.data, self.Array)

    def test_merge(self):
        result = xr.merge([self.x1, self.x2], compat="override", join="outer")
        assert isinstance(result.foo.data, self.Array)

    def test_where(self):
        x1, x2 = xr.align(self.x1, self.x2, join="inner")
        result = xr.where(x1 > 2, x1, x2)
        assert isinstance(result.data, self.Array)

    def test_full_like(self):
        result = xr.full_like(self.x1, 0)
        assert isinstance(result.data, self.Array)

    def test_cov(self):
        result = xr.cov(self.x1, self.x2)
        assert isinstance(result.data, self.Array)

    def test_corr(self):
        result = xr.corr(self.x1, self.x2)
        assert isinstance(result.data, self.Array)

    def test_cross(self):
        x1, x2 = xr.align(self.x1.squeeze(), self.x2.squeeze(), join="inner")
        result = xr.cross(x1, x2, dim="x")
        assert isinstance(result.data, self.Array)

    def test_dot(self):
        result = xr.dot(self.x1, self.x2)
        assert isinstance(result.data, self.Array)

    def test_map_blocks(self):
        result = xr.map_blocks(lambda x: x + 1, self.x1)
        assert isinstance(result.data, self.Array)


@pytest.mark.parametrize("namespace", NAMESPACE_ARRAYS)
class TestDataArrayMethods(_BaseTest):
    @pytest.fixture(autouse=True)
    def setUp(self, request, namespace):
        self.setup_for_test(request, namespace)
        self.x = self.get_test_dataarray()

    def test_loc(self):
        result = self.x.loc[{"x": slice(1, 3)}]
        assert isinstance(result.data, self.Array)

    def test_isel(self):
        result = self.x.isel(x=slice(1, 3))
        assert isinstance(result.data, self.Array)

    def test_sel(self):
        result = self.x.sel(x=slice(1, 3))
        assert isinstance(result.data, self.Array)

    def test_squeeze(self):
        result = self.x.squeeze("y")
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="interp uses numpy and scipy")
    def test_interp(self):
        # TODO: some cases could be made to work
        result = self.x.interp(x=2.5)
        assert isinstance(result.data, self.Array)

    def test_isnull(self):
        result = self.x.isnull()
        assert isinstance(result.data, self.Array)

    def test_notnull(self):
        result = self.x.notnull()
        assert isinstance(result.data, self.Array)

    def test_count(self):
        result = self.x.count()
        assert isinstance(result.data, self.Array)

    def test_dropna(self):
        result = self.x.dropna(dim="x")
        assert isinstance(result.data, self.Array)

    def test_fillna(self):
        result = self.x.fillna(0)
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="ffill uses bottleneck or numbagg")
    def test_ffill(self):
        result = self.x.ffill()
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="bfill uses bottleneck or numbagg")
    def test_bfill(self):
        result = self.x.bfill()
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="interpolate_na uses numpy and scipy")
    def test_interpolate_na(self):
        result = self.x.interpolate_na()
        assert isinstance(result.data, self.Array)

    def test_where(self):
        result = self.x.where(self.x > 2)
        assert isinstance(result.data, self.Array)

    def test_isin(self):
        test_elements = self.constructor(np.asarray([1]))
        result = self.x.isin(test_elements)
        assert isinstance(result.data, self.Array)

    def test_groupby(self):
        result = self.x.groupby("x").mean()
        assert isinstance(result.data, self.Array)

    def test_groupby_bins(self):
        result = self.x.groupby_bins("x", bins=[0, 2, 4, 6]).mean()
        assert isinstance(result.data, self.Array)

    def test_rolling_iter(self):
        result = self.x.rolling(x=3)
        elem = next(iter(result))[1]
        assert isinstance(elem.data, self.Array)

    def test_rolling_construct(self):
        result = self.x.rolling(x=3).construct(x="window")
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_rolling_reduce(self, skipna):
        result = self.x.rolling(x=3).mean(skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="rolling_exp uses numbagg")
    def test_rolling_exp_reduce(self):
        result = self.x.rolling_exp(x=3).mean()
        assert isinstance(result.data, self.Array)

    def test_cumulative_iter(self):
        result = self.x.cumulative("x")
        elem = next(iter(result))[1]
        assert isinstance(elem.data, self.Array)

    def test_cumulative_construct(self):
        result = self.x.cumulative("x").construct(x="window")
        assert isinstance(result.data, self.Array)

    def test_cumulative_reduce(self):
        result = self.x.cumulative("x").sum()
        assert isinstance(result.data, self.Array)

    def test_weighted(self):
        result = self.x.weighted(self.x.fillna(0)).mean()
        assert isinstance(result.data, self.Array)

    def test_coarsen_construct(self):
        result = self.x.coarsen(x=2, boundary="pad").construct(x=["a", "b"])
        assert isinstance(result.data, self.Array)

    def test_coarsen_reduce(self):
        result = self.x.coarsen(x=2, boundary="pad").mean()
        assert isinstance(result.data, self.Array)

    def test_resample(self):
        time_coord = pd.date_range("2000-01-01", periods=5)
        result = self.x.assign_coords(x=time_coord).resample(x="D").mean()
        assert isinstance(result.data, self.Array)

    def test_diff(self):
        result = self.x.diff("x")
        assert isinstance(result.data, self.Array)

    def test_dot(self):
        result = self.x.dot(self.x)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_quantile(self, skipna):
        result = self.x.quantile(0.5, skipna=skipna)
        assert isinstance(result.data, self.Array)

    def test_differentiate(self):
        # edge_order is not implemented in jax, and only supports passing None
        edge_order = None if self.namespace == "jax.numpy" else 1
        result = self.x.differentiate("x", edge_order=edge_order)
        assert isinstance(result.data, self.Array)

    def test_integrate(self):
        result = self.x.integrate("x")
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="polyfit uses numpy linalg")
    def test_polyfit(self):
        # TODO: this could work, there are just a lot of different linalg calls
        result = self.x.polyfit("x", 1)
        assert isinstance(result.polyfit_coefficients.data, self.Array)

    def test_map_blocks(self):
        result = self.x.map_blocks(lambda x: x + 1)
        assert isinstance(result.data, self.Array)

    def test_all(self):
        result = self.x.all(dim="x")
        assert isinstance(result.data, self.Array)

    def test_any(self):
        result = self.x.any(dim="x")
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_argmax(self, skipna):
        result = self.x.argmax(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_argmin(self, skipna):
        result = self.x.argmin(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_idxmax(self, skipna):
        result = self.x.idxmax(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_idxmin(self, skipna):
        result = self.x.idxmin(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_max(self, skipna):
        result = self.x.max(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_min(self, skipna):
        result = self.x.min(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_mean(self, skipna):
        result = self.x.mean(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_median(self, skipna):
        result = self.x.median(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_prod(self, skipna):
        result = self.x.prod(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_sum(self, skipna):
        result = self.x.sum(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_std(self, skipna):
        result = self.x.std(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_var(self, skipna):
        result = self.x.var(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_cumsum(self, skipna):
        result = self.x.cumsum(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_cumprod(self, skipna):
        result = self.x.cumprod(dim="x", skipna=skipna)
        assert isinstance(result.data, self.Array)

    def test_argsort(self):
        result = self.x.argsort()
        assert isinstance(result.data, self.Array)

    def test_astype(self):
        result = self.x.astype(int)
        assert isinstance(result.data, self.Array)

    def test_clip(self):
        result = self.x.clip(min=2.0, max=4.0)
        assert isinstance(result.data, self.Array)

    def test_conj(self):
        result = self.x.conj()
        assert isinstance(result.data, self.Array)

    def test_conjugate(self):
        result = self.x.conjugate()
        assert isinstance(result.data, self.Array)

    def test_imag(self):
        result = self.x.imag
        assert isinstance(result.data, self.Array)

    def test_searchsorted(self):
        v = self.constructor(np.asarray([3]))
        result = self.x.squeeze().searchsorted(v)
        assert isinstance(result, self.Array)

    def test_round(self):
        result = self.x.round()
        assert isinstance(result.data, self.Array)

    def test_real(self):
        result = self.x.real
        assert isinstance(result.data, self.Array)

    def test_T(self):
        result = self.x.T
        assert isinstance(result.data, self.Array)

    @pytest.mark.xfail(reason="rank uses bottleneck")
    def test_rank(self):
        # TODO: scipy has rankdata, as does jax, so this can work
        result = self.x.rank()
        assert isinstance(result.data, self.Array)

    def test_transpose(self):
        result = self.x.transpose()
        assert isinstance(result.data, self.Array)

    def test_stack(self):
        result = self.x.stack(z=("x", "y"))
        assert isinstance(result.data, self.Array)

    def test_unstack(self):
        result = self.x.stack(z=("x", "y")).unstack("z")
        assert isinstance(result.data, self.Array)

    def test_shift(self):
        result = self.x.shift(x=1)
        assert isinstance(result.data, self.Array)

    def test_roll(self):
        result = self.x.roll(x=1)
        assert isinstance(result.data, self.Array)

    def test_pad(self):
        result = self.x.pad(x=1)
        assert isinstance(result.data, self.Array)

    def test_sortby(self):
        result = self.x.sortby("x")
        assert isinstance(result.data, self.Array)

    def test_broadcast_like(self):
        result = self.x.broadcast_like(self.x)
        assert isinstance(result.data, self.Array)
