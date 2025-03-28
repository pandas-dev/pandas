from __future__ import annotations

import contextlib
import string
import warnings

from dask._compatibility import import_optional_dependency

import_optional_dependency("pandas")
import_optional_dependency("numpy")
import_optional_dependency("pyarrow")

import numpy as np
import pandas as pd
import pyarrow as pa
from packaging.version import Version

PANDAS_VERSION = Version(pd.__version__)
PANDAS_GE_201 = PANDAS_VERSION.release >= (2, 0, 1)
PANDAS_GE_202 = PANDAS_VERSION.release >= (2, 0, 2)
PANDAS_GE_210 = PANDAS_VERSION.release >= (2, 1, 0)
PANDAS_GE_211 = PANDAS_VERSION.release >= (2, 1, 1)
PANDAS_GE_220 = PANDAS_VERSION.release >= (2, 2, 0)
PANDAS_GE_230 = PANDAS_VERSION.release >= (2, 3, 0)
PANDAS_GE_300 = PANDAS_VERSION.major >= 3

PYARROW_VERSION = Version(pa.__version__)
PYARROW_GE_1500 = PYARROW_VERSION.release >= (15, 0, 0)

import pandas.testing as tm


def assert_categorical_equal(left, right, *args, **kwargs):
    tm.assert_extension_array_equal(left, right, *args, **kwargs)
    assert isinstance(
        left.dtype, pd.CategoricalDtype
    ), f"{left} is not categorical dtype"
    assert isinstance(
        right.dtype, pd.CategoricalDtype
    ), f"{right} is not categorical dtype"


def assert_numpy_array_equal(left, right):
    left_na = pd.isna(left)
    right_na = pd.isna(right)
    np.testing.assert_array_equal(left_na, right_na)

    left_valid = left[~left_na]
    right_valid = right[~right_na]
    np.testing.assert_array_equal(left_valid, right_valid)


def makeDataFrame():
    data = np.random.randn(30, 4)
    index = list(string.ascii_letters)[:30]
    return pd.DataFrame(data, index=index, columns=list("ABCD"))


def makeTimeDataFrame():
    data = makeDataFrame()
    data.index = makeDateIndex()
    return data


def makeTimeSeries():
    return makeTimeDataFrame()["A"]


def makeDateIndex(k=30, freq="B"):
    return pd.date_range("2000", periods=k, freq=freq)


def makeTimedeltaIndex(k=30, freq="D"):
    return pd.timedelta_range("1 day", periods=k, freq=freq)


def makeMissingDataframe():
    df = makeDataFrame()
    data = df.values
    data = np.where(data > 1, np.nan, data)
    return pd.DataFrame(data, index=df.index, columns=df.columns)


def makeMixedDataFrame():
    df = pd.DataFrame(
        {
            "A": [0.0, 1, 2, 3, 4],
            "B": [0.0, 1, 0, 1, 0],
            "C": [f"foo{i}" for i in range(5)],
            "D": pd.date_range("2009-01-01", periods=5),
        }
    )
    return df


@contextlib.contextmanager
def check_groupby_axis_deprecation():
    if PANDAS_GE_210:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                ".*Call without passing 'axis' instead|.*Operate on the un-grouped DataFrame instead",
                FutureWarning,
            )
            yield
    else:
        yield


@contextlib.contextmanager
def check_observed_deprecation():
    if PANDAS_GE_210:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The default of observed=False",
                category=FutureWarning,
            )
            yield
    else:
        yield


@contextlib.contextmanager
def check_convert_dtype_deprecation():
    if PANDAS_GE_210:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="the convert_dtype parameter",
                category=FutureWarning,
            )
            yield
    else:
        yield


@contextlib.contextmanager
def check_apply_dataframe_deprecation():
    if PANDAS_GE_210:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Returning a DataFrame",
                category=FutureWarning,
            )
            yield
    else:
        yield


@contextlib.contextmanager
def check_reductions_runtime_warning():
    if not PANDAS_GE_201:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="invalid value encountered in double_scalars|Degrees of freedom <= 0 for slice",
                category=RuntimeWarning,
            )
            yield
    else:
        yield


IndexingError = pd.errors.IndexingError


def is_any_real_numeric_dtype(arr_or_dtype) -> bool:
    return pd.api.types.is_any_real_numeric_dtype(arr_or_dtype)


def is_string_dtype(arr_or_dtype) -> bool:
    # is_string_dtype did not recognize pyarrow strings before 2.0
    # Can remove once 2.0 is minimum version for us
    if hasattr(arr_or_dtype, "dtype"):
        dtype = arr_or_dtype.dtype
    else:
        dtype = arr_or_dtype
    return pd.api.types.is_string_dtype(dtype)
