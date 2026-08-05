import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    _testing as tm,
    date_range,
    period_range,
    read_hdf,
)

pytestmark = pytest.mark.single_cpu


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_store_datetime_fractional_secs(temp_hdfstore, unit):
    dt = datetime.datetime(2012, 1, 2, 3, 4, 5, 123456)
    dti = DatetimeIndex([dt], dtype=f"M8[{unit}]")
    series = Series([0], index=dti)
    temp_hdfstore["a"] = series
    assert temp_hdfstore["a"].index[0] == dt


@pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
def test_tseries_indices_series(temp_hdfstore):
    idx = date_range("2020-01-01", periods=10)
    ser = Series(np.random.default_rng(2).standard_normal(len(idx)), idx)
    temp_hdfstore["a"] = ser
    result = temp_hdfstore["a"]

    tm.assert_series_equal(result, ser)
    assert result.index.freq == ser.index.freq
    tm.assert_class_equal(result.index, ser.index, obj="series index")

    idx = period_range("2020-01-01", periods=10, freq="D")
    ser = Series(np.random.default_rng(2).standard_normal(len(idx)), idx)
    temp_hdfstore["a"] = ser
    result = temp_hdfstore["a"]

    tm.assert_series_equal(result, ser)
    assert result.index.freq == ser.index.freq
    tm.assert_class_equal(result.index, ser.index, obj="series index")


@pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
def test_tseries_indices_frame(temp_hdfstore):
    idx = date_range("2020-01-01", periods=10)
    df = DataFrame(np.random.default_rng(2).standard_normal((len(idx), 3)), index=idx)
    temp_hdfstore["a"] = df
    result = temp_hdfstore["a"]

    tm.assert_frame_equal(result, df)
    assert result.index.freq == df.index.freq
    tm.assert_class_equal(result.index, df.index, obj="dataframe index")

    idx = period_range("2020-01-01", periods=10, freq="D")
    df = DataFrame(np.random.default_rng(2).standard_normal((len(idx), 3)), idx)
    temp_hdfstore["a"] = df
    result = temp_hdfstore["a"]

    tm.assert_frame_equal(result, df)
    assert result.index.freq == df.index.freq
    tm.assert_class_equal(result.index, df.index, obj="dataframe index")


def test_to_hdf_dateoffset_freq(temp_file):
    # GH#45790: pytables stores attrs at pickle protocol 0; DateOffset
    #  (i.e. RelativeDeltaOffset) previously failed at that protocol
    df = DataFrame(
        {"a": [0]},
        index=DatetimeIndex(["2020-01-01"], freq=pd.DateOffset(years=1)),
    )
    df.to_hdf(temp_file, key="data")
    result = read_hdf(temp_file)
    tm.assert_frame_equal(result, df)
    assert result.index.freq == df.index.freq


def test_read_hdf_freq_bytes_attr(temp_file):
    # GH#35917: HDF5 files written by old (Py2-era) pandas stored freq as
    #  a Python 2 pickle byte-string that pytables can't unpickle in Py3.
    #  Reading should warn and drop the freq, not propagate the bytes.
    import tables

    old_freq_bytes = (
        b"ccopy_reg\n_reconstructor\np0\n(cpandas.tseries.offsets\nDay\n"
        b"p1\nc__builtin__\nobject\np2\nNtp3\nRp4\n(dp5\nVn\np6\nI1\n"
        b"sVnormalize\np7\nI00\nsV_cache\np8\n(dp9\nsb."
    )
    df = DataFrame({"a": [1, 2, 3]}, index=date_range("2020-01-01", periods=3))
    df.to_hdf(temp_file, key="data", format="fixed")

    with tables.open_file(temp_file, "a") as f:
        f.root.data.axis1._v_attrs.freq = old_freq_bytes
        # pytables returns this as np.bytes_ on read-back; the guard's
        #  isinstance(..., bytes) check covers it because np.bytes_ subclasses
        #  bytes. Pin that assumption so a future pytables change is caught here.
        assert isinstance(f.root.data.axis1._v_attrs["freq"], np.bytes_)

    with tm.assert_produces_warning(
        UserWarning, match="Could not decode freq attribute"
    ):
        result = read_hdf(temp_file)
    assert result.index.freq is None
    tm.assert_index_equal(result.index, df.index._with_freq(None))
