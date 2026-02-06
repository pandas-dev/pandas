import collections
from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DatetimeIndex,
    Index,
    Interval,
    IntervalIndex,
    MultiIndex,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    array,
)
import pandas._testing as tm
from pandas.tests.base.common import allow_na_ops


@pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
def test_value_counts(index_or_series_obj):
    obj = index_or_series_obj
    obj = np.repeat(obj, range(1, len(obj) + 1))
    result = obj.value_counts()

    counter = collections.Counter(obj)
    expected = Series(dict(counter.most_common()), dtype=np.int64, name="count")

    if obj.dtype != np.float16:
        expected.index = expected.index.astype(obj.dtype)
    else:
        with pytest.raises(NotImplementedError, match="float16 indexes are not "):
            expected.index.astype(obj.dtype)
        return
    if isinstance(expected.index, MultiIndex):
        expected.index.names = obj.names
    else:
        expected.index.name = obj.name

    if not isinstance(result.dtype, np.dtype):
        if getattr(obj.dtype, "storage", "") == "pyarrow":
            expected = expected.astype("int64[pyarrow]")
        else:
            # i.e IntegerDtype
            expected = expected.astype("Int64")

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("null_obj", [np.nan, None])
@pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
def test_value_counts_null(null_obj, index_or_series_obj):
    orig = index_or_series_obj

    if not allow_na_ops(orig):
        pytest.skip("type doesn't allow for NA operations")
    elif len(orig) < 1:
        pytest.skip("Test doesn't make sense on empty data")
    elif isinstance(orig, MultiIndex):
        pytest.skip(f"MultiIndex can't hold '{null_obj}'")

    obj = orig.copy(deep=True)
    values = obj._values
    values[0:2] = null_obj

    klass = type(obj)
    repeated_values = np.repeat(values, range(1, len(values) + 1))
    obj = klass(repeated_values, dtype=obj.dtype)

    # because np.nan == np.nan is False, but None == None is True
    # np.nan would be duplicated, whereas None wouldn't
    counter = collections.Counter(obj.dropna())
    expected = Series(dict(counter.most_common()), dtype=np.int64, name="count")

    if obj.dtype != np.float16:
        expected.index = expected.index.astype(obj.dtype)
    else:
        with pytest.raises(NotImplementedError, match="float16 indexes are not "):
            expected.index.astype(obj.dtype)
        return
    expected.index.name = obj.name

    result = obj.value_counts()

    if not isinstance(result.dtype, np.dtype):
        if getattr(obj.dtype, "storage", "") == "pyarrow":
            expected = expected.astype("int64[pyarrow]")
        else:
            # i.e IntegerDtype
            expected = expected.astype("Int64")
    tm.assert_series_equal(result, expected)

    expected[null_obj] = 3

    result = obj.value_counts(dropna=False)
    expected = expected.sort_index()
    result = result.sort_index()
    tm.assert_series_equal(result, expected)


def test_value_counts_inferred(index_or_series, using_infer_string):
    klass = index_or_series
    s_values = ["a", "b", "b", "b", "b", "c", "d", "d", "a", "a"]
    s = klass(s_values)
    expected = Series([4, 3, 2, 1], index=["b", "a", "d", "c"], name="count")
    tm.assert_series_equal(s.value_counts(), expected)

    if isinstance(s, Index):
        exp = Index(np.unique(np.array(s_values, dtype=np.object_)))
        tm.assert_index_equal(s.unique(), exp)
    else:
        exp = np.unique(np.array(s_values, dtype=np.object_))
        if using_infer_string:
            exp = array(exp, dtype="str")
        tm.assert_equal(s.unique(), exp)

    assert s.nunique() == 4
    # don't sort, have to sort after the fact as not sorting is
    # platform-dep
    hist = s.value_counts(sort=False).sort_values()
    expected = Series([3, 1, 4, 2], index=list("acbd"), name="count").sort_values()
    tm.assert_series_equal(hist, expected)

    # sort ascending
    hist = s.value_counts(ascending=True)
    expected = Series([1, 2, 3, 4], index=list("cdab"), name="count")
    tm.assert_series_equal(hist, expected)

    # relative histogram.
    hist = s.value_counts(normalize=True)
    expected = Series(
        [0.4, 0.3, 0.2, 0.1], index=["b", "a", "d", "c"], name="proportion"
    )
    tm.assert_series_equal(hist, expected)


def test_value_counts_bins(index_or_series, using_infer_string):
    klass = index_or_series
    s_values = ["a", "b", "b", "b", "b", "c", "d", "d", "a", "a"]
    s = klass(s_values)

    # bins
    msg = "bins argument only works with numeric data"
    with pytest.raises(TypeError, match=msg):
        s.value_counts(bins=1)

    s1 = Series([1, 1, 2, 3])
    res1 = s1.value_counts(bins=1)
    exp1 = Series({Interval(0.997, 3.0): 4}, name="count")
    tm.assert_series_equal(res1, exp1)
    res1n = s1.value_counts(bins=1, normalize=True)
    exp1n = Series({Interval(0.997, 3.0): 1.0}, name="proportion")
    tm.assert_series_equal(res1n, exp1n)

    if isinstance(s1, Index):
        tm.assert_index_equal(s1.unique(), Index([1, 2, 3]))
    else:
        exp = np.array([1, 2, 3], dtype=np.int64)
        tm.assert_numpy_array_equal(s1.unique(), exp)

    assert s1.nunique() == 3

    # these return the same
    res4 = s1.value_counts(bins=4, dropna=True)
    intervals = IntervalIndex.from_breaks([0.997, 1.5, 2.0, 2.5, 3.0])
    exp4 = Series([2, 1, 1, 0], index=intervals.take([0, 1, 3, 2]), name="count")
    tm.assert_series_equal(res4, exp4)

    res4 = s1.value_counts(bins=4, dropna=False)
    intervals = IntervalIndex.from_breaks([0.997, 1.5, 2.0, 2.5, 3.0])
    exp4 = Series([2, 1, 1, 0], index=intervals.take([0, 1, 3, 2]), name="count")
    tm.assert_series_equal(res4, exp4)

    res4n = s1.value_counts(bins=4, normalize=True)
    exp4n = Series(
        [0.5, 0.25, 0.25, 0], index=intervals.take([0, 1, 3, 2]), name="proportion"
    )
    tm.assert_series_equal(res4n, exp4n)

    # handle NA's properly
    s_values = ["a", "b", "b", "b", np.nan, np.nan, "d", "d", "a", "a", "b"]
    s = klass(s_values)
    expected = Series([4, 3, 2], index=["b", "a", "d"], name="count")
    tm.assert_series_equal(s.value_counts(), expected)

    if isinstance(s, Index):
        exp = Index(["a", "b", np.nan, "d"])
        tm.assert_index_equal(s.unique(), exp)
    else:
        exp = np.array(["a", "b", np.nan, "d"], dtype=object)
        if using_infer_string:
            exp = array(exp, dtype="str")
        tm.assert_equal(s.unique(), exp)
    assert s.nunique() == 3

    s = klass({}) if klass is dict else klass({}, dtype=object)
    expected = Series([], dtype=np.int64, name="count")
    tm.assert_series_equal(s.value_counts(), expected, check_index_type=False)
    # returned dtype differs depending on original
    if isinstance(s, Index):
        tm.assert_index_equal(s.unique(), Index([]), exact=False)
    else:
        tm.assert_numpy_array_equal(s.unique(), np.array([]), check_dtype=False)

    assert s.nunique() == 0


def test_value_counts_datetime64(index_or_series, unit):
    klass = index_or_series

    # GH 3002, datetime64[ns]
    # don't test names though
    df = pd.DataFrame(
        {
            "person_id": ["xxyyzz", "xxyyzz", "xxyyzz", "xxyyww", "foofoo", "foofoo"],
            "dt": pd.to_datetime(
                [
                    "2010-01-01",
                    "2010-01-01",
                    "2010-01-01",
                    "2009-01-01",
                    "2008-09-09",
                    "2008-09-09",
                ]
            ).as_unit(unit),
            "food": ["PIE", "GUM", "EGG", "EGG", "PIE", "GUM"],
        }
    )

    s = klass(df["dt"].copy())
    s.name = None
    idx = pd.to_datetime(
        ["2010-01-01 00:00:00", "2008-09-09 00:00:00", "2009-01-01 00:00:00"]
    ).as_unit(unit)
    expected_s = Series([3, 2, 1], index=idx, name="count")
    tm.assert_series_equal(s.value_counts(), expected_s)

    expected = array(
        np.array(
            ["2010-01-01 00:00:00", "2009-01-01 00:00:00", "2008-09-09 00:00:00"],
            dtype=f"datetime64[{unit}]",
        )
    )
    result = s.unique()
    if isinstance(s, Index):
        tm.assert_index_equal(result, DatetimeIndex(expected))
    else:
        tm.assert_extension_array_equal(result, expected)

    assert s.nunique() == 3

    # with NaT
    s = df["dt"].copy()
    s = klass(list(s.values) + [pd.NaT] * 4)
    if klass is Series:
        s = s.dt.as_unit(unit)
    else:
        s = s.as_unit(unit)

    result = s.value_counts()
    assert result.index.dtype == f"datetime64[{unit}]"
    tm.assert_series_equal(result, expected_s)

    result = s.value_counts(dropna=False)
    expected_s = pd.concat(
        [
            Series([4], index=DatetimeIndex([pd.NaT]).as_unit(unit), name="count"),
            expected_s,
        ]
    )
    tm.assert_series_equal(result, expected_s)

    assert s.dtype == f"datetime64[{unit}]"
    unique = s.unique()
    assert unique.dtype == f"datetime64[{unit}]"

    # numpy_array_equal cannot compare pd.NaT
    if isinstance(s, Index):
        exp_idx = DatetimeIndex([*expected.tolist(), pd.NaT]).as_unit(unit)
        tm.assert_index_equal(unique, exp_idx)
    else:
        tm.assert_extension_array_equal(unique[:3], expected)
        assert pd.isna(unique[3])

    assert s.nunique() == 3
    assert s.nunique(dropna=False) == 4


def test_value_counts_timedelta64(index_or_series, unit):
    # timedelta64[ns]
    klass = index_or_series

    day = Timedelta(timedelta(1)).as_unit(unit)
    tdi = TimedeltaIndex([day], name="dt").as_unit(unit)

    tdvals = np.zeros(6, dtype=f"m8[{unit}]") + day
    td = klass(tdvals, name="dt")

    result = td.value_counts()
    expected_s = Series([6], index=tdi, name="count")
    tm.assert_series_equal(result, expected_s)

    expected = tdi
    result = td.unique()
    if isinstance(td, Index):
        tm.assert_index_equal(result, expected)
    else:
        tm.assert_extension_array_equal(result, expected._values)

    td2 = day + np.zeros(6, dtype=f"m8[{unit}]")
    td2 = klass(td2, name="dt")
    result2 = td2.value_counts()
    tm.assert_series_equal(result2, expected_s)


def test_value_counts_with_nan(dropna, index_or_series):
    # GH31944
    klass = index_or_series
    values = [True, pd.NA, np.nan]
    obj = klass(values)
    res = obj.value_counts(dropna=dropna)
    if dropna is True:
        expected = Series([1], index=Index([True], dtype=obj.dtype), name="count")
    else:
        expected = Series([1, 1, 1], index=[True, pd.NA, np.nan], name="count")
    tm.assert_series_equal(res, expected)


def test_value_counts_object_inference_deprecated():
    # GH#56161
    dti = pd.date_range("2016-01-01", periods=3, tz="UTC")

    idx = dti.astype(object)
    res = idx.value_counts()

    exp = dti.value_counts()
    exp.index = exp.index.astype(object)
    tm.assert_series_equal(res, exp)


@pytest.mark.parametrize(
    ("index", "expected_index"),
    [
        [
            pd.date_range("2016-01-01", periods=5, freq="D", unit="ns"),
            pd.date_range("2016-01-01", periods=5, freq="D", unit="ns"),
        ],
        [
            pd.timedelta_range(Timedelta(0), periods=5, freq="h"),
            pd.timedelta_range(Timedelta(0), periods=5, freq="h"),
        ],
        [
            DatetimeIndex(
                [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(1)]
                + [Timestamp("2016-01-02")]
                + [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(1, 5)]
            ),
            DatetimeIndex(pd.date_range("2016-01-01", periods=5, freq="D", unit="us")),
        ],
        [
            TimedeltaIndex(
                [Timedelta(hours=i) for i in range(1)]
                + [Timedelta(hours=1)]
                + [Timedelta(hours=i) for i in range(1, 5)],
            ),
            TimedeltaIndex(pd.timedelta_range(Timedelta(hours=0), periods=5, freq="h")),
        ],
        [
            DatetimeIndex(
                [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(2)]
                + [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(3, 5)],
            ),
            DatetimeIndex(
                [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(2)]
                + [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(3, 5)],
            ),
        ],
        [
            TimedeltaIndex(
                [Timedelta(hours=i) for i in range(2)]
                + [Timedelta(hours=i) for i in range(3, 5)],
            ),
            TimedeltaIndex(
                [Timedelta(hours=i) for i in range(2)]
                + [Timedelta(hours=i) for i in range(3, 5)],
            ),
        ],
        [
            DatetimeIndex(
                [Timestamp("2016-01-01")]
                + [pd.NaT]
                + [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(1, 5)],
            ),
            DatetimeIndex(
                [Timestamp("2016-01-01")]
                + [pd.NaT]
                + [Timestamp("2016-01-01") + Timedelta(days=i) for i in range(1, 5)],
            ),
        ],
        [
            TimedeltaIndex(
                [Timedelta(hours=0)]
                + [pd.NaT]
                + [Timedelta(hours=i) for i in range(1, 5)],
            ),
            TimedeltaIndex(
                [Timedelta(hours=0)]
                + [pd.NaT]
                + [Timedelta(hours=i) for i in range(1, 5)],
            ),
        ],
    ],
)
def test_value_counts_index_datetimelike(index, expected_index):
    vc = index.value_counts(sort=False, dropna=False)
    tm.assert_index_equal(vc.index, expected_index)
