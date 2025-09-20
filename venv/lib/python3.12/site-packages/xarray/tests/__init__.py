from __future__ import annotations

import importlib
import platform
import string
import warnings
from contextlib import contextmanager, nullcontext
from unittest import mock  # noqa: F401

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_equal  # noqa: F401
from packaging.version import Version
from pandas.testing import assert_frame_equal  # noqa: F401

import xarray.testing
from xarray import Dataset
from xarray.coding.times import _STANDARD_CALENDARS as _STANDARD_CALENDARS_UNSORTED
from xarray.core.duck_array_ops import allclose_or_equiv  # noqa: F401
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.options import set_options
from xarray.core.variable import IndexVariable
from xarray.testing import (  # noqa: F401
    assert_chunks_equal,
    assert_duckarray_allclose,
    assert_duckarray_equal,
)
from xarray.tests.arrays import (  # noqa: F401
    ConcatenatableArray,
    DuckArrayWrapper,
    FirstElementAccessibleArray,
    InaccessibleArray,
    UnexpectedDataAccess,
)

# import mpl and change the backend before other mpl imports
try:
    import matplotlib as mpl

    # Order of imports is important here.
    # Using a different backend makes Travis CI work
    mpl.use("Agg")
except ImportError:
    pass

# https://github.com/pydata/xarray/issues/7322
warnings.filterwarnings("ignore", "'urllib3.contrib.pyopenssl' module is deprecated")
warnings.filterwarnings("ignore", "Deprecated call to `pkg_resources.declare_namespace")
warnings.filterwarnings("ignore", "pkg_resources is deprecated as an API")
warnings.filterwarnings("ignore", message="numpy.ndarray size changed")

arm_xfail = pytest.mark.xfail(
    platform.machine() == "aarch64" or "arm" in platform.machine(),
    reason="expected failure on ARM",
)


def assert_writeable(ds):
    readonly = [
        name
        for name, var in ds.variables.items()
        if not isinstance(var, IndexVariable)
        and not isinstance(
            var.data, PandasExtensionArray | pd.api.extensions.ExtensionArray
        )
        and not var.data.flags.writeable
    ]
    assert not readonly, readonly


def _importorskip(
    modname: str, minversion: str | None = None
) -> tuple[bool, pytest.MarkDecorator]:
    try:
        mod = importlib.import_module(modname)
        has = True
        if minversion is not None:
            v = getattr(mod, "__version__", "999")
            if Version(v) < Version(minversion):
                raise ImportError("Minimum version not satisfied")
    except ImportError:
        has = False

    reason = f"requires {modname}"
    if minversion is not None:
        reason += f">={minversion}"
    func = pytest.mark.skipif(not has, reason=reason)
    return has, func


has_matplotlib, requires_matplotlib = _importorskip("matplotlib")
has_scipy, requires_scipy = _importorskip("scipy")
has_scipy_ge_1_13, requires_scipy_ge_1_13 = _importorskip("scipy", "1.13")
with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="'cgi' is deprecated and slated for removal in Python 3.13",
        category=DeprecationWarning,
    )
    has_pydap, requires_pydap = _importorskip("pydap.client")
has_netCDF4, requires_netCDF4 = _importorskip("netCDF4")
with warnings.catch_warnings():
    # see https://github.com/pydata/xarray/issues/8537
    warnings.filterwarnings(
        "ignore",
        message="h5py is running against HDF5 1.14.3",
        category=UserWarning,
    )

    has_h5netcdf, requires_h5netcdf = _importorskip("h5netcdf")
has_cftime, requires_cftime = _importorskip("cftime")
has_dask, requires_dask = _importorskip("dask")
has_dask_ge_2024_08_1, requires_dask_ge_2024_08_1 = _importorskip(
    "dask", minversion="2024.08.1"
)
has_dask_ge_2024_11_0, requires_dask_ge_2024_11_0 = _importorskip("dask", "2024.11.0")
has_dask_ge_2025_1_0, requires_dask_ge_2025_1_0 = _importorskip("dask", "2025.1.0")
if has_dask_ge_2025_1_0:
    has_dask_expr = True
    requires_dask_expr = pytest.mark.skipif(not has_dask_expr, reason="should not skip")
else:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="The current Dask DataFrame implementation is deprecated.",
            category=DeprecationWarning,
        )
        has_dask_expr, requires_dask_expr = _importorskip("dask_expr")
has_bottleneck, requires_bottleneck = _importorskip("bottleneck")
has_rasterio, requires_rasterio = _importorskip("rasterio")
has_zarr, requires_zarr = _importorskip("zarr")
has_zarr_v3, requires_zarr_v3 = _importorskip("zarr", "3.0.0")
has_zarr_v3_dtypes, requires_zarr_v3_dtypes = _importorskip("zarr", "3.1.0")
has_zarr_v3_async_oindex, requires_zarr_v3_async_oindex = _importorskip("zarr", "3.1.2")
if has_zarr_v3:
    import zarr

    # manual update by checking attrs for now
    # TODO: use version specifier
    # installing from git main is giving me a lower version than the
    # most recently released zarr
    has_zarr_v3_dtypes = hasattr(zarr.core, "dtype")
    has_zarr_v3_async_oindex = hasattr(zarr.AsyncArray, "oindex")

    requires_zarr_v3_dtypes = pytest.mark.skipif(
        not has_zarr_v3_dtypes, reason="requires zarr>3.1.0"
    )
    requires_zarr_v3_async_oindex = pytest.mark.skipif(
        not has_zarr_v3_async_oindex, reason="requires zarr>3.1.1"
    )


has_fsspec, requires_fsspec = _importorskip("fsspec")
has_iris, requires_iris = _importorskip("iris")
has_numbagg, requires_numbagg = _importorskip("numbagg")
has_pyarrow, requires_pyarrow = _importorskip("pyarrow")
with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="is_categorical_dtype is deprecated and will be removed in a future version.",
        category=DeprecationWarning,
    )
    # seaborn uses the deprecated `pandas.is_categorical_dtype`
    has_seaborn, requires_seaborn = _importorskip("seaborn")
has_sparse, requires_sparse = _importorskip("sparse")
has_cupy, requires_cupy = _importorskip("cupy")
has_cartopy, requires_cartopy = _importorskip("cartopy")
has_pint, requires_pint = _importorskip("pint")
has_numexpr, requires_numexpr = _importorskip("numexpr")
has_flox, requires_flox = _importorskip("flox")
has_netcdf, requires_netcdf = _importorskip("netcdf")
has_pandas_ge_2_2, requires_pandas_ge_2_2 = _importorskip("pandas", "2.2")
has_pandas_3, requires_pandas_3 = _importorskip("pandas", "3.0.0.dev0")


# some special cases
has_scipy_or_netCDF4 = has_scipy or has_netCDF4
requires_scipy_or_netCDF4 = pytest.mark.skipif(
    not has_scipy_or_netCDF4, reason="requires scipy or netCDF4"
)
has_numbagg_or_bottleneck = has_numbagg or has_bottleneck
requires_numbagg_or_bottleneck = pytest.mark.skipif(
    not has_numbagg_or_bottleneck, reason="requires numbagg or bottleneck"
)
has_numpy_2, requires_numpy_2 = _importorskip("numpy", "2.0.0")
has_flox_0_9_12, requires_flox_0_9_12 = _importorskip("flox", "0.9.12")

has_array_api_strict, requires_array_api_strict = _importorskip("array_api_strict")

parametrize_zarr_format = pytest.mark.parametrize(
    "zarr_format",
    [
        pytest.param(2, id="zarr_format=2"),
        pytest.param(
            3,
            marks=pytest.mark.skipif(
                not has_zarr_v3,
                reason="zarr-python v2 cannot understand the zarr v3 format",
            ),
            id="zarr_format=3",
        ),
    ],
)


def _importorskip_h5netcdf_ros3(has_h5netcdf: bool):
    if not has_h5netcdf:
        return has_h5netcdf, pytest.mark.skipif(
            not has_h5netcdf, reason="requires h5netcdf"
        )

    import h5py

    h5py_with_ros3 = h5py.get_config().ros3

    return h5py_with_ros3, pytest.mark.skipif(
        not h5py_with_ros3,
        reason="requires h5netcdf>=1.3.0 and h5py with ros3 support",
    )


has_h5netcdf_ros3, requires_h5netcdf_ros3 = _importorskip_h5netcdf_ros3(has_h5netcdf)
has_netCDF4_1_6_2_or_above, requires_netCDF4_1_6_2_or_above = _importorskip(
    "netCDF4", "1.6.2"
)

has_h5netcdf_1_4_0_or_above, requires_h5netcdf_1_4_0_or_above = _importorskip(
    "h5netcdf", "1.4.0.dev"
)

has_netCDF4_1_7_0_or_above, requires_netCDF4_1_7_0_or_above = _importorskip(
    "netCDF4", "1.7.0"
)

# change some global options for tests
set_options(warn_for_unclosed_files=True)

if has_dask:
    import dask


class CountingScheduler:
    """Simple dask scheduler counting the number of computes.

    Reference: https://stackoverflow.com/questions/53289286/"""

    def __init__(self, max_computes=0):
        self.total_computes = 0
        self.max_computes = max_computes

    def __call__(self, dsk, keys, **kwargs):
        self.total_computes += 1
        if self.total_computes > self.max_computes:
            raise RuntimeError(
                f"Too many computes. Total: {self.total_computes} > max: {self.max_computes}."
            )
        return dask.get(dsk, keys, **kwargs)


def raise_if_dask_computes(max_computes=0):
    # return a dummy context manager so that this can be used for non-dask objects
    if not has_dask:
        return nullcontext()
    scheduler = CountingScheduler(max_computes)
    return dask.config.set(scheduler=scheduler)


flaky = pytest.mark.flaky
network = pytest.mark.network


class ReturnItem:
    def __getitem__(self, key):
        return key


class IndexerMaker:
    def __init__(self, indexer_cls):
        self._indexer_cls = indexer_cls

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        return self._indexer_cls(key)


def source_ndarray(array):
    """Given an ndarray, return the base object which holds its memory, or the
    object itself.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "DatetimeIndex.base")
        warnings.filterwarnings("ignore", "TimedeltaIndex.base")
        base = getattr(array, "base", np.asarray(array).base)
    if base is None:
        base = array
    return base


def format_record(record) -> str:
    """Format warning record like `FutureWarning('Function will be deprecated...')`"""
    return f"{str(record.category)[8:-2]}('{record.message}'))"


@contextmanager
def assert_no_warnings():
    with warnings.catch_warnings(record=True) as record:
        yield record
        assert len(record) == 0, (
            f"Got {len(record)} unexpected warning(s): {[format_record(r) for r in record]}"
        )


# Internal versions of xarray's test functions that validate additional
# invariants


def assert_equal(a, b, check_default_indexes=True):
    __tracebackhide__ = True
    xarray.testing.assert_equal(a, b)
    xarray.testing._assert_internal_invariants(a, check_default_indexes)
    xarray.testing._assert_internal_invariants(b, check_default_indexes)


def assert_identical(a, b, check_default_indexes=True):
    __tracebackhide__ = True
    xarray.testing.assert_identical(a, b)
    xarray.testing._assert_internal_invariants(a, check_default_indexes)
    xarray.testing._assert_internal_invariants(b, check_default_indexes)


def assert_allclose(a, b, check_default_indexes=True, **kwargs):
    __tracebackhide__ = True
    xarray.testing.assert_allclose(a, b, **kwargs)
    xarray.testing._assert_internal_invariants(a, check_default_indexes)
    xarray.testing._assert_internal_invariants(b, check_default_indexes)


_DEFAULT_TEST_DIM_SIZES = (8, 9, 10)


def create_test_data(
    seed: int = 12345,
    add_attrs: bool = True,
    dim_sizes: tuple[int, int, int] = _DEFAULT_TEST_DIM_SIZES,
    use_extension_array: bool = False,
) -> Dataset:
    rs = np.random.default_rng(seed)
    _vars = {
        "var1": ["dim1", "dim2"],
        "var2": ["dim1", "dim2"],
        "var3": ["dim3", "dim1"],
    }
    _dims = {"dim1": dim_sizes[0], "dim2": dim_sizes[1], "dim3": dim_sizes[2]}

    obj = Dataset()
    obj["dim2"] = ("dim2", 0.5 * np.arange(_dims["dim2"]))
    if _dims["dim3"] > 26:
        raise RuntimeError(
            f"Not enough letters for filling this dimension size ({_dims['dim3']})"
        )
    obj["dim3"] = ("dim3", list(string.ascii_lowercase[0 : _dims["dim3"]]))
    obj["time"] = (
        "time",
        pd.date_range(
            "2000-01-01",
            periods=20,
            unit="ns",
        ),
    )
    for v, dims in sorted(_vars.items()):
        data = rs.normal(size=tuple(_dims[d] for d in dims))
        obj[v] = (dims, data)
        if add_attrs:
            obj[v].attrs = {"foo": "variable"}
    if use_extension_array:
        obj["var4"] = (
            "dim1",
            pd.Categorical(
                rs.choice(
                    list(string.ascii_lowercase[: rs.integers(1, 5)]),
                    size=dim_sizes[0],
                )
            ),
        )
        if has_pyarrow:
            obj["var5"] = (
                "dim1",
                pd.array(
                    rs.integers(1, 10, size=dim_sizes[0]).tolist(),
                    dtype="int64[pyarrow]",
                ),
            )
    if dim_sizes == _DEFAULT_TEST_DIM_SIZES:
        numbers_values = np.array([0, 1, 2, 0, 0, 1, 1, 2, 2, 3], dtype="int64")
    else:
        numbers_values = rs.integers(0, 3, _dims["dim3"], dtype="int64")
    obj.coords["numbers"] = ("dim3", numbers_values)
    obj.encoding = {"foo": "bar"}
    assert_writeable(obj)
    return obj


_STANDARD_CALENDAR_NAMES = sorted(_STANDARD_CALENDARS_UNSORTED)
_NON_STANDARD_CALENDAR_NAMES = {
    "noleap",
    "365_day",
    "360_day",
    "julian",
    "all_leap",
    "366_day",
}
_NON_STANDARD_CALENDARS = [
    pytest.param(cal, marks=requires_cftime)
    for cal in sorted(_NON_STANDARD_CALENDAR_NAMES)
]
_STANDARD_CALENDARS = [
    pytest.param(cal, marks=requires_cftime if cal != "standard" else ())
    for cal in _STANDARD_CALENDAR_NAMES
]
_ALL_CALENDARS = sorted(_STANDARD_CALENDARS + _NON_STANDARD_CALENDARS)
_CFTIME_CALENDARS = [
    pytest.param(*p.values, marks=requires_cftime) for p in _ALL_CALENDARS
]


def _all_cftime_date_types():
    import cftime

    return {
        "noleap": cftime.DatetimeNoLeap,
        "365_day": cftime.DatetimeNoLeap,
        "360_day": cftime.Datetime360Day,
        "julian": cftime.DatetimeJulian,
        "all_leap": cftime.DatetimeAllLeap,
        "366_day": cftime.DatetimeAllLeap,
        "gregorian": cftime.DatetimeGregorian,
        "proleptic_gregorian": cftime.DatetimeProlepticGregorian,
    }
