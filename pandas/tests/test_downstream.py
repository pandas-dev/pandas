"""
Testing that we work in the downstream packages
"""

import array
from functools import partial
import importlib
import subprocess
import sys

import numpy as np
import pytest

from pandas.errors import IntCastingNaNError

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    TimedeltaIndex,
)
import pandas._testing as tm
from pandas.util.version import Version


@pytest.fixture
def df():
    return DataFrame({"A": [1, 2, 3]})


def test_dask(df):
    # dask sets "compute.use_numexpr" to False, so catch the current value
    # and ensure to reset it afterwards to avoid impacting other tests
    olduse = pd.get_option("compute.use_numexpr")

    try:
        pytest.importorskip("toolz")
        dd = pytest.importorskip("dask.dataframe")

        ddf = dd.from_pandas(df, npartitions=3)
        assert ddf.A is not None
        assert ddf.compute() is not None
    finally:
        pd.set_option("compute.use_numexpr", olduse)


# TODO(CoW) see https://github.com/pandas-dev/pandas/pull/51082
@pytest.mark.skip(reason="not implemented with CoW")
def test_dask_ufunc():
    # dask sets "compute.use_numexpr" to False, so catch the current value
    # and ensure to reset it afterwards to avoid impacting other tests
    olduse = pd.get_option("compute.use_numexpr")

    try:
        da = pytest.importorskip("dask.array")
        dd = pytest.importorskip("dask.dataframe")

        s = Series([1.5, 2.3, 3.7, 4.0])
        ds = dd.from_pandas(s, npartitions=2)

        result = da.log(ds).compute()
        expected = np.log(s)
        tm.assert_series_equal(result, expected)
    finally:
        pd.set_option("compute.use_numexpr", olduse)


def test_construct_dask_float_array_int_dtype_match_ndarray():
    # GH#40110 make sure we treat a float-dtype dask array with the same
    #  rules we would for an ndarray
    dd = pytest.importorskip("dask.dataframe")

    arr = np.array([1, 2.5, 3])
    darr = dd.from_array(arr)

    res = Series(darr)
    expected = Series(arr)
    tm.assert_series_equal(res, expected)

    # GH#49599 in 2.0 we raise instead of silently ignoring the dtype
    msg = "Trying to coerce float values to integers"
    with pytest.raises(ValueError, match=msg):
        Series(darr, dtype="i8")

    msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
    arr[2] = np.nan
    with pytest.raises(IntCastingNaNError, match=msg):
        Series(darr, dtype="i8")
    # which is the same as we get with a numpy input
    with pytest.raises(IntCastingNaNError, match=msg):
        Series(arr, dtype="i8")


def test_xarray(df):
    pytest.importorskip("xarray")

    assert df.to_xarray() is not None


def test_xarray_cftimeindex_nearest():
    # https://github.com/pydata/xarray/issues/3751
    cftime = pytest.importorskip("cftime")
    xarray = pytest.importorskip("xarray")

    times = xarray.date_range("0001", periods=2, use_cftime=True)
    key = cftime.DatetimeGregorian(2000, 1, 1)
    result = times.get_indexer([key], method="nearest")
    expected = 1
    assert result == expected


@pytest.mark.single_cpu
def test_oo_optimizable():
    # GH 21071
    subprocess.check_call([sys.executable, "-OO", "-c", "import pandas"])


@pytest.mark.single_cpu
def test_oo_optimized_datetime_index_unpickle():
    # GH 42866
    subprocess.check_call(
        [
            sys.executable,
            "-OO",
            "-c",
            (
                "import pandas as pd, pickle; "
                "pickle.loads(pickle.dumps(pd.date_range('2021-01-01', periods=1)))"
            ),
        ]
    )


def test_statsmodels():
    smf = pytest.importorskip("statsmodels.formula.api")

    df = DataFrame(
        {"Lottery": range(5), "Literacy": range(5), "Pop1831": range(100, 105)}
    )
    smf.ols("Lottery ~ Literacy + np.log(Pop1831)", data=df).fit()


def test_scikit_learn():
    pytest.importorskip("sklearn")
    from sklearn import (
        datasets,
        svm,
    )

    digits = datasets.load_digits()
    clf = svm.SVC(gamma=0.001, C=100.0)
    clf.fit(digits.data[:-1], digits.target[:-1])
    clf.predict(digits.data[-1:])


def test_seaborn(mpl_cleanup):
    seaborn = pytest.importorskip("seaborn")
    tips = DataFrame(
        {"day": pd.date_range("2023", freq="D", periods=5), "total_bill": range(5)}
    )
    seaborn.stripplot(x="day", y="total_bill", data=tips)


@pytest.mark.xfail(reason="pandas_datareader uses old variant of deprecate_kwarg")
def test_pandas_datareader():
    # https://github.com/pandas-dev/pandas/pull/61468
    # https://github.com/pydata/pandas-datareader/issues/1005
    pytest.importorskip("pandas_datareader")


@pytest.mark.filterwarnings("ignore:Passing a BlockManager:DeprecationWarning")
def test_pyarrow(df):
    pyarrow = pytest.importorskip("pyarrow")
    table = pyarrow.Table.from_pandas(df)
    result = table.to_pandas()
    tm.assert_frame_equal(result, df)


def test_yaml_dump(df):
    # GH#42748
    yaml = pytest.importorskip("yaml")

    dumped = yaml.dump(df)

    loaded = yaml.load(dumped, Loader=yaml.Loader)
    tm.assert_frame_equal(df, loaded)

    loaded2 = yaml.load(dumped, Loader=yaml.UnsafeLoader)
    tm.assert_frame_equal(df, loaded2)


@pytest.mark.parametrize("dependency", ["numpy", "dateutil"])
def test_missing_required_dependency(monkeypatch, dependency):
    # GH#61030
    original_import = __import__
    mock_error = ImportError(f"Mock error for {dependency}")

    def mock_import(name, *args, **kwargs):
        if name == dependency:
            raise mock_error
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", mock_import)

    with pytest.raises(ImportError, match=dependency):
        importlib.reload(importlib.import_module("pandas"))


def test_frame_setitem_dask_array_into_new_col(request):
    # GH#47128

    # dask sets "compute.use_numexpr" to False, so catch the current value
    # and ensure to reset it afterwards to avoid impacting other tests
    olduse = pd.get_option("compute.use_numexpr")

    try:
        dask = pytest.importorskip("dask")
        da = pytest.importorskip("dask.array")
        if Version(dask.__version__) <= Version("2025.1.0") and Version(
            np.__version__
        ) >= Version("2.1"):
            request.applymarker(
                pytest.mark.xfail(reason="loc.__setitem__ incorrectly mutated column c")
            )

        dda = da.array([1, 2])
        df = DataFrame({"a": ["a", "b"]})
        df["b"] = dda
        df["c"] = dda
        df.loc[[False, True], "b"] = 100
        result = df.loc[[1], :]
        expected = DataFrame({"a": ["b"], "b": [100], "c": [2]}, index=[1])
        tm.assert_frame_equal(result, expected)
    finally:
        pd.set_option("compute.use_numexpr", olduse)


def test_pandas_priority():
    # GH#48347

    class MyClass:
        __pandas_priority__ = 5000

        def __radd__(self, other):
            return self

    left = MyClass()
    right = Series(range(3))

    assert right.__add__(left) is NotImplemented
    assert right + left is left


@pytest.mark.parametrize("dtype", ["M8[ns]", "m8[ns]"])
@pytest.mark.parametrize(
    "box", [memoryview, partial(array.array, "i"), "dask", "xarray"]
)
def test_from_obscure_array(dtype, box):
    # GH#24539 recognize e.g xarray, dask, ...
    # Note: we dont do this for PeriodArray bc _from_sequence won't accept
    #  an array of integers
    # TODO: could check with arraylike of Period objects
    # GH#24539 recognize e.g xarray, dask, ...
    arr = np.array([1, 2, 3], dtype=np.int64)
    if box == "dask":
        da = pytest.importorskip("dask.array")
        data = da.array(arr)
    elif box == "xarray":
        xr = pytest.importorskip("xarray")
        data = xr.DataArray(arr)
    else:
        data = box(arr)

    func = {"M8[ns]": pd.to_datetime, "m8[ns]": pd.to_timedelta}[dtype]
    result = func(arr).array
    expected = func(data).array
    tm.assert_equal(result, expected)

    # Let's check the Indexes while we're here
    idx_cls = {"M8[ns]": DatetimeIndex, "m8[ns]": TimedeltaIndex}[dtype]
    result = idx_cls(arr)
    expected = idx_cls(data)
    tm.assert_index_equal(result, expected)


def test_xarray_coerce_unit():
    # GH44053
    xr = pytest.importorskip("xarray")

    arr = xr.DataArray([1, 2, 3])
    result = pd.to_datetime(arr, unit="ns")
    expected = DatetimeIndex(
        [
            "1970-01-01 00:00:00.000000001",
            "1970-01-01 00:00:00.000000002",
            "1970-01-01 00:00:00.000000003",
        ],
        dtype="datetime64[ns]",
        freq=None,
    )
    tm.assert_index_equal(result, expected)
