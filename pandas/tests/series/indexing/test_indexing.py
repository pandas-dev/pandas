""" test get/set & misc """

from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.common import is_scalar

import pandas as pd
from pandas import Categorical, DataFrame, MultiIndex, Series, Timedelta, Timestamp
import pandas._testing as tm

from pandas.tseries.offsets import BDay


def test_basic_indexing():
    s = Series(np.random.randn(5), index=["a", "b", "a", "a", "b"])

    msg = "index 5 is out of bounds for axis 0 with size 5"
    with pytest.raises(IndexError, match=msg):
        s[5]
    with pytest.raises(IndexError, match=msg):
        s[5] = 0

    with pytest.raises(KeyError, match=r"^'c'$"):
        s["c"]

    s = s.sort_index()

    with pytest.raises(IndexError, match=msg):
        s[5]
    msg = r"index 5 is out of bounds for axis (0|1) with size 5|^5$"
    with pytest.raises(IndexError, match=msg):
        s[5] = 0


def test_basic_getitem_with_labels(datetime_series):
    indices = datetime_series.index[[5, 10, 15]]

    result = datetime_series[indices]
    expected = datetime_series.reindex(indices)
    tm.assert_series_equal(result, expected)

    result = datetime_series[indices[0] : indices[2]]
    expected = datetime_series.loc[indices[0] : indices[2]]
    tm.assert_series_equal(result, expected)

    # integer indexes, be careful
    s = Series(np.random.randn(10), index=list(range(0, 20, 2)))
    inds = [0, 2, 5, 7, 8]
    arr_inds = np.array([0, 2, 5, 7, 8])
    with pytest.raises(KeyError, match="with any missing labels"):
        s[inds]

    with pytest.raises(KeyError, match="with any missing labels"):
        s[arr_inds]

    # GH12089
    # with tz for values
    s = Series(
        pd.date_range("2011-01-01", periods=3, tz="US/Eastern"), index=["a", "b", "c"]
    )
    expected = Timestamp("2011-01-01", tz="US/Eastern")
    result = s.loc["a"]
    assert result == expected
    result = s.iloc[0]
    assert result == expected
    result = s["a"]
    assert result == expected


def test_getitem_setitem_ellipsis():
    s = Series(np.random.randn(10))

    np.fix(s)

    result = s[...]
    tm.assert_series_equal(result, s)

    s[...] = 5
    assert (result == 5).all()


def test_getitem_get(datetime_series, string_series, object_series):
    idx1 = string_series.index[5]
    idx2 = object_series.index[5]

    assert string_series[idx1] == string_series.get(idx1)
    assert object_series[idx2] == object_series.get(idx2)

    assert string_series[idx1] == string_series[5]
    assert object_series[idx2] == object_series[5]

    assert string_series.get(-1) == string_series.get(string_series.index[-1])
    assert string_series[5] == string_series.get(string_series.index[5])

    # missing
    d = datetime_series.index[0] - BDay()
    msg = r"Timestamp\('1999-12-31 00:00:00', freq='B'\)"
    with pytest.raises(KeyError, match=msg):
        datetime_series[d]

    # None
    # GH 5652
    s1 = Series(dtype=object)
    s2 = Series(dtype=object, index=list("abc"))
    for s in [s1, s2]:
        result = s.get(None)
        assert result is None


def test_getitem_fancy(string_series, object_series):
    slice1 = string_series[[1, 2, 3]]
    slice2 = object_series[[1, 2, 3]]
    assert string_series.index[2] == slice1.index[1]
    assert object_series.index[2] == slice2.index[1]
    assert string_series[2] == slice1[1]
    assert object_series[2] == slice2[1]


def test_getitem_generator(string_series):
    gen = (x > 0 for x in string_series)
    result = string_series[gen]
    result2 = string_series[iter(string_series > 0)]
    expected = string_series[string_series > 0]
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result2, expected)


def test_type_promotion():
    # GH12599
    s = pd.Series(dtype=object)
    s["a"] = pd.Timestamp("2016-01-01")
    s["b"] = 3.0
    s["c"] = "foo"
    expected = Series([pd.Timestamp("2016-01-01"), 3.0, "foo"], index=["a", "b", "c"])
    tm.assert_series_equal(s, expected)


@pytest.mark.parametrize(
    "result_1, duplicate_item, expected_1",
    [
        [
            pd.Series({1: 12, 2: [1, 2, 2, 3]}),
            pd.Series({1: 313}),
            pd.Series({1: 12}, dtype=object),
        ],
        [
            pd.Series({1: [1, 2, 3], 2: [1, 2, 2, 3]}),
            pd.Series({1: [1, 2, 3]}),
            pd.Series({1: [1, 2, 3]}),
        ],
    ],
)
def test_getitem_with_duplicates_indices(result_1, duplicate_item, expected_1):
    # GH 17610
    result = result_1.append(duplicate_item)
    expected = expected_1.append(duplicate_item)
    tm.assert_series_equal(result[1], expected)
    assert result[2] == result_1[2]


def test_getitem_out_of_bounds(datetime_series):
    # don't segfault, GH #495
    msg = r"index \d+ is out of bounds for axis 0 with size \d+"
    with pytest.raises(IndexError, match=msg):
        datetime_series[len(datetime_series)]

    # GH #917
    msg = r"index -\d+ is out of bounds for axis 0 with size \d+"
    s = Series([], dtype=object)
    with pytest.raises(IndexError, match=msg):
        s[-1]


def test_getitem_setitem_integers():
    # caused bug without test
    s = Series([1, 2, 3], ["a", "b", "c"])

    assert s.iloc[0] == s["a"]
    s.iloc[0] = 5
    tm.assert_almost_equal(s["a"], 5)


def test_getitem_box_float64(datetime_series):
    value = datetime_series[5]
    assert isinstance(value, np.float64)


@pytest.mark.parametrize(
    "arr",
    [np.random.randn(10), tm.makeDateIndex(10, name="a").tz_localize(tz="US/Eastern")],
)
def test_get(arr):
    # GH 21260
    s = Series(arr, index=[2 * i for i in range(len(arr))])
    assert s.get(4) == s.iloc[2]

    result = s.get([4, 6])
    expected = s.iloc[[2, 3]]
    tm.assert_series_equal(result, expected)

    result = s.get(slice(2))
    expected = s.iloc[[0, 1]]
    tm.assert_series_equal(result, expected)

    assert s.get(-1) is None
    assert s.get(s.index.max() + 1) is None

    s = Series(arr[:6], index=list("abcdef"))
    assert s.get("c") == s.iloc[2]

    result = s.get(slice("b", "d"))
    expected = s.iloc[[1, 2, 3]]
    tm.assert_series_equal(result, expected)

    result = s.get("Z")
    assert result is None

    assert s.get(4) == s.iloc[4]
    assert s.get(-1) == s.iloc[-1]
    assert s.get(len(s)) is None

    # GH 21257
    s = pd.Series(arr)
    s2 = s[::2]
    assert s2.get(1) is None


def test_series_box_timestamp():
    rng = pd.date_range("20090415", "20090519", freq="B")
    ser = Series(rng)

    assert isinstance(ser[5], pd.Timestamp)

    rng = pd.date_range("20090415", "20090519", freq="B")
    ser = Series(rng, index=rng)
    assert isinstance(ser[5], pd.Timestamp)

    assert isinstance(ser.iat[5], pd.Timestamp)


def test_series_box_timedelta():
    rng = pd.timedelta_range("1 day 1 s", periods=5, freq="h")
    ser = pd.Series(rng)
    assert isinstance(ser[0], Timedelta)
    assert isinstance(ser.at[1], Timedelta)
    assert isinstance(ser.iat[2], Timedelta)
    assert isinstance(ser.loc[3], Timedelta)
    assert isinstance(ser.iloc[4], Timedelta)


def test_getitem_ambiguous_keyerror():
    s = Series(range(10), index=list(range(0, 20, 2)))
    with pytest.raises(KeyError, match=r"^1$"):
        s[1]
    with pytest.raises(KeyError, match=r"^1$"):
        s.loc[1]


def test_getitem_unordered_dup():
    obj = Series(range(5), index=["c", "a", "a", "b", "b"])
    assert is_scalar(obj["c"])
    assert obj["c"] == 0


def test_getitem_dups_with_missing():
    # breaks reindex, so need to use .loc internally
    # GH 4246
    s = Series([1, 2, 3, 4], ["foo", "bar", "foo", "bah"])
    with pytest.raises(KeyError, match="with any missing labels"):
        s.loc[["foo", "bar", "bah", "bam"]]

    with pytest.raises(KeyError, match="with any missing labels"):
        s[["foo", "bar", "bah", "bam"]]


def test_getitem_dups():
    s = Series(range(5), index=["A", "A", "B", "C", "C"], dtype=np.int64)
    expected = Series([3, 4], index=["C", "C"], dtype=np.int64)
    result = s["C"]
    tm.assert_series_equal(result, expected)


def test_setitem_ambiguous_keyerror():
    s = Series(range(10), index=list(range(0, 20, 2)))

    # equivalent of an append
    s2 = s.copy()
    s2[1] = 5
    expected = s.append(Series([5], index=[1]))
    tm.assert_series_equal(s2, expected)

    s2 = s.copy()
    s2.loc[1] = 5
    expected = s.append(Series([5], index=[1]))
    tm.assert_series_equal(s2, expected)


def test_getitem_dataframe():
    rng = list(range(10))
    s = pd.Series(10, index=rng)
    df = pd.DataFrame(rng, index=rng)
    msg = (
        "Indexing a Series with DataFrame is not supported, "
        "use the appropriate DataFrame column"
    )
    with pytest.raises(TypeError, match=msg):
        s[df > 5]


def test_setitem(datetime_series, string_series):
    datetime_series[datetime_series.index[5]] = np.NaN
    datetime_series[[1, 2, 17]] = np.NaN
    datetime_series[6] = np.NaN
    assert np.isnan(datetime_series[6])
    assert np.isnan(datetime_series[2])
    datetime_series[np.isnan(datetime_series)] = 5
    assert not np.isnan(datetime_series[2])

    # caught this bug when writing tests
    series = Series(tm.makeIntIndex(20).astype(float), index=tm.makeIntIndex(20))

    series[::2] = 0
    assert (series[::2] == 0).all()

    # set item that's not contained
    s = string_series.copy()
    s["foobar"] = 1

    app = Series([1], index=["foobar"], name="series")
    expected = string_series.append(app)
    tm.assert_series_equal(s, expected)

    # Test for issue #10193
    key = pd.Timestamp("2012-01-01")
    series = pd.Series(dtype=object)
    series[key] = 47
    expected = pd.Series(47, [key])
    tm.assert_series_equal(series, expected)

    series = pd.Series([], pd.DatetimeIndex([], freq="D"), dtype=object)
    series[key] = 47
    expected = pd.Series(47, pd.DatetimeIndex([key], freq="D"))
    tm.assert_series_equal(series, expected)


def test_setitem_dtypes():
    # change dtypes
    # GH 4463
    expected = Series([np.nan, 2, 3])

    s = Series([1, 2, 3])
    s.iloc[0] = np.nan
    tm.assert_series_equal(s, expected)

    s = Series([1, 2, 3])
    s.loc[0] = np.nan
    tm.assert_series_equal(s, expected)

    s = Series([1, 2, 3])
    s[0] = np.nan
    tm.assert_series_equal(s, expected)

    s = Series([False])
    s.loc[0] = np.nan
    tm.assert_series_equal(s, Series([np.nan]))

    s = Series([False, True])
    s.loc[0] = np.nan
    tm.assert_series_equal(s, Series([np.nan, 1.0]))


def test_set_value(datetime_series, string_series):
    idx = datetime_series.index[10]
    res = datetime_series._set_value(idx, 0)
    assert res is None
    assert datetime_series[idx] == 0

    # equiv
    s = string_series.copy()
    res = s._set_value("foobar", 0)
    assert res is None
    assert s.index[-1] == "foobar"
    assert s["foobar"] == 0

    s = string_series.copy()
    s.loc["foobar"] = 0
    assert s.index[-1] == "foobar"
    assert s["foobar"] == 0


def test_setslice(datetime_series):
    sl = datetime_series[5:20]
    assert len(sl) == len(sl.index)
    assert sl.index.is_unique is True


def test_2d_to_1d_assignment_raises():
    x = np.random.randn(2, 2)
    y = pd.Series(range(2))

    msg = (
        r"shape mismatch: value array of shape \(2,2\) could not be "
        r"broadcast to indexing result of shape \(2,\)"
    )
    with pytest.raises(ValueError, match=msg):
        y.loc[range(2)] = x

    msg = r"could not broadcast input array from shape \(2,2\) into shape \(2\)"
    with pytest.raises(ValueError, match=msg):
        y.loc[:] = x


# FutureWarning from NumPy about [slice(None, 5).
@pytest.mark.filterwarnings("ignore:Using a non-tuple:FutureWarning")
def test_basic_getitem_setitem_corner(datetime_series):
    # invalid tuples, e.g. td.ts[:, None] vs. td.ts[:, 2]
    msg = "Can only tuple-index with a MultiIndex"
    with pytest.raises(ValueError, match=msg):
        datetime_series[:, 2]
    with pytest.raises(ValueError, match=msg):
        datetime_series[:, 2] = 2

    # weird lists. [slice(0, 5)] will work but not two slices
    with tm.assert_produces_warning(FutureWarning):
        # GH#31299
        result = datetime_series[[slice(None, 5)]]
    expected = datetime_series[:5]
    tm.assert_series_equal(result, expected)

    # OK
    msg = r"unhashable type(: 'slice')?"
    with pytest.raises(TypeError, match=msg):
        datetime_series[[5, slice(None, None)]]
    with pytest.raises(TypeError, match=msg):
        datetime_series[[5, slice(None, None)]] = 2


@pytest.mark.parametrize("tz", ["US/Eastern", "UTC", "Asia/Tokyo"])
def test_setitem_with_tz(tz):
    orig = pd.Series(pd.date_range("2016-01-01", freq="H", periods=3, tz=tz))
    assert orig.dtype == f"datetime64[ns, {tz}]"

    # scalar
    s = orig.copy()
    s[1] = pd.Timestamp("2011-01-01", tz=tz)
    exp = pd.Series(
        [
            pd.Timestamp("2016-01-01 00:00", tz=tz),
            pd.Timestamp("2011-01-01 00:00", tz=tz),
            pd.Timestamp("2016-01-01 02:00", tz=tz),
        ]
    )
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[1] = pd.Timestamp("2011-01-01", tz=tz)
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[1] = pd.Timestamp("2011-01-01", tz=tz)
    tm.assert_series_equal(s, exp)

    # vector
    vals = pd.Series(
        [pd.Timestamp("2011-01-01", tz=tz), pd.Timestamp("2012-01-01", tz=tz)],
        index=[1, 2],
    )
    assert vals.dtype == f"datetime64[ns, {tz}]"

    s[[1, 2]] = vals
    exp = pd.Series(
        [
            pd.Timestamp("2016-01-01 00:00", tz=tz),
            pd.Timestamp("2011-01-01 00:00", tz=tz),
            pd.Timestamp("2012-01-01 00:00", tz=tz),
        ]
    )
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)


def test_setitem_with_tz_dst():
    # GH XXX
    tz = "US/Eastern"
    orig = pd.Series(pd.date_range("2016-11-06", freq="H", periods=3, tz=tz))
    assert orig.dtype == f"datetime64[ns, {tz}]"

    # scalar
    s = orig.copy()
    s[1] = pd.Timestamp("2011-01-01", tz=tz)
    exp = pd.Series(
        [
            pd.Timestamp("2016-11-06 00:00-04:00", tz=tz),
            pd.Timestamp("2011-01-01 00:00-05:00", tz=tz),
            pd.Timestamp("2016-11-06 01:00-05:00", tz=tz),
        ]
    )
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[1] = pd.Timestamp("2011-01-01", tz=tz)
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[1] = pd.Timestamp("2011-01-01", tz=tz)
    tm.assert_series_equal(s, exp)

    # vector
    vals = pd.Series(
        [pd.Timestamp("2011-01-01", tz=tz), pd.Timestamp("2012-01-01", tz=tz)],
        index=[1, 2],
    )
    assert vals.dtype == f"datetime64[ns, {tz}]"

    s[[1, 2]] = vals
    exp = pd.Series(
        [
            pd.Timestamp("2016-11-06 00:00", tz=tz),
            pd.Timestamp("2011-01-01 00:00", tz=tz),
            pd.Timestamp("2012-01-01 00:00", tz=tz),
        ]
    )
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)


def test_categorical_assigning_ops():
    orig = Series(Categorical(["b", "b"], categories=["a", "b"]))
    s = orig.copy()
    s[:] = "a"
    exp = Series(Categorical(["a", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[1] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[s.index > 0] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[[False, True]] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.index = ["x", "y"]
    s["y"] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]), index=["x", "y"])
    tm.assert_series_equal(s, exp)

    # ensure that one can set something to np.nan
    s = Series(Categorical([1, 2, 3]))
    exp = Series(Categorical([1, np.nan, 3], categories=[1, 2, 3]))
    s[1] = np.nan
    tm.assert_series_equal(s, exp)


def test_getitem_categorical_str():
    # GH#31765
    ser = pd.Series(range(5), index=pd.Categorical(["a", "b", "c", "a", "b"]))
    result = ser["a"]
    expected = ser.iloc[[0, 3]]
    tm.assert_series_equal(result, expected)

    # Check the intermediate steps work as expected
    result = ser.index.get_value(ser, "a")
    tm.assert_series_equal(result, expected)


def test_slice(string_series, object_series):
    numSlice = string_series[10:20]
    numSliceEnd = string_series[-10:]
    objSlice = object_series[10:20]

    assert string_series.index[9] not in numSlice.index
    assert object_series.index[9] not in objSlice.index

    assert len(numSlice) == len(numSlice.index)
    assert string_series[numSlice.index[0]] == numSlice[numSlice.index[0]]

    assert numSlice.index[1] == string_series.index[11]
    assert tm.equalContents(numSliceEnd, np.array(string_series)[-10:])

    # Test return view.
    sl = string_series[10:20]
    sl[:] = 0

    assert (string_series[10:20] == 0).all()


def test_slice_can_reorder_not_uniquely_indexed():
    s = Series(1, index=["a", "a", "b", "b", "c"])
    s[::-1]  # it works!


def test_loc_setitem(string_series):
    inds = string_series.index[[3, 4, 7]]

    result = string_series.copy()
    result.loc[inds] = 5

    expected = string_series.copy()
    expected[[3, 4, 7]] = 5
    tm.assert_series_equal(result, expected)

    result.iloc[5:10] = 10
    expected[5:10] = 10
    tm.assert_series_equal(result, expected)

    # set slice with indices
    d1, d2 = string_series.index[[5, 15]]
    result.loc[d1:d2] = 6
    expected[5:16] = 6  # because it's inclusive
    tm.assert_series_equal(result, expected)

    # set index value
    string_series.loc[d1] = 4
    string_series.loc[d2] = 6
    assert string_series[d1] == 4
    assert string_series[d2] == 6


def test_setitem_na():
    # these induce dtype changes
    expected = Series([np.nan, 3, np.nan, 5, np.nan, 7, np.nan, 9, np.nan])
    s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])
    s[::2] = np.nan
    tm.assert_series_equal(s, expected)

    # gets coerced to float, right?
    expected = Series([np.nan, 1, np.nan, 0])
    s = Series([True, True, False, False])
    s[::2] = np.nan
    tm.assert_series_equal(s, expected)

    expected = Series([np.nan, np.nan, np.nan, np.nan, np.nan, 5, 6, 7, 8, 9])
    s = Series(np.arange(10))
    s[:5] = np.nan
    tm.assert_series_equal(s, expected)


def test_timedelta_assignment():
    # GH 8209
    s = Series([], dtype=object)
    s.loc["B"] = timedelta(1)
    tm.assert_series_equal(s, Series(Timedelta("1 days"), index=["B"]))

    s = s.reindex(s.index.insert(0, "A"))
    tm.assert_series_equal(s, Series([np.nan, Timedelta("1 days")], index=["A", "B"]))

    result = s.fillna(timedelta(1))
    expected = Series(Timedelta("1 days"), index=["A", "B"])
    tm.assert_series_equal(result, expected)

    s.loc["A"] = timedelta(1)
    tm.assert_series_equal(s, expected)

    # GH 14155
    s = Series(10 * [np.timedelta64(10, "m")])
    s.loc[[1, 2, 3]] = np.timedelta64(20, "m")
    expected = pd.Series(10 * [np.timedelta64(10, "m")])
    expected.loc[[1, 2, 3]] = pd.Timedelta(np.timedelta64(20, "m"))
    tm.assert_series_equal(s, expected)


@pytest.mark.parametrize(
    "nat_val,should_cast",
    [
        (pd.NaT, True),
        (np.timedelta64("NaT", "ns"), False),
        (np.datetime64("NaT", "ns"), True),
    ],
)
@pytest.mark.parametrize("tz", [None, "UTC"])
def test_dt64_series_assign_nat(nat_val, should_cast, tz):
    # some nat-like values should be cast to datetime64 when inserting
    #  into a datetime64 series.  Others should coerce to object
    #  and retain their dtypes.
    dti = pd.date_range("2016-01-01", periods=3, tz=tz)
    base = pd.Series(dti)
    expected = pd.Series([pd.NaT] + list(dti[1:]), dtype=dti.dtype)
    if not should_cast:
        expected = expected.astype(object)

    ser = base.copy(deep=True)
    ser[0] = nat_val
    tm.assert_series_equal(ser, expected)

    ser = base.copy(deep=True)
    ser.loc[0] = nat_val
    tm.assert_series_equal(ser, expected)

    ser = base.copy(deep=True)
    ser.iloc[0] = nat_val
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize(
    "nat_val,should_cast",
    [
        (pd.NaT, True),
        (np.timedelta64("NaT", "ns"), True),
        (np.datetime64("NaT", "ns"), False),
    ],
)
def test_td64_series_assign_nat(nat_val, should_cast):
    # some nat-like values should be cast to timedelta64 when inserting
    #  into a timedelta64 series.  Others should coerce to object
    #  and retain their dtypes.
    base = pd.Series([0, 1, 2], dtype="m8[ns]")
    expected = pd.Series([pd.NaT, 1, 2], dtype="m8[ns]")
    if not should_cast:
        expected = expected.astype(object)

    ser = base.copy(deep=True)
    ser[0] = nat_val
    tm.assert_series_equal(ser, expected)

    ser = base.copy(deep=True)
    ser.loc[0] = nat_val
    tm.assert_series_equal(ser, expected)

    ser = base.copy(deep=True)
    ser.iloc[0] = nat_val
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize(
    "td",
    [
        pd.Timedelta("9 days"),
        pd.Timedelta("9 days").to_timedelta64(),
        pd.Timedelta("9 days").to_pytimedelta(),
    ],
)
def test_append_timedelta_does_not_cast(td):
    # GH#22717 inserting a Timedelta should _not_ cast to int64
    expected = pd.Series(["x", td], index=[0, "td"], dtype=object)

    ser = pd.Series(["x"])
    ser["td"] = td
    tm.assert_series_equal(ser, expected)
    assert isinstance(ser["td"], pd.Timedelta)

    ser = pd.Series(["x"])
    ser.loc["td"] = pd.Timedelta("9 days")
    tm.assert_series_equal(ser, expected)
    assert isinstance(ser["td"], pd.Timedelta)


def test_underlying_data_conversion():
    # GH 4080
    df = DataFrame({c: [1, 2, 3] for c in ["a", "b", "c"]})
    df.set_index(["a", "b", "c"], inplace=True)
    s = Series([1], index=[(2, 2, 2)])
    df["val"] = 0
    df
    df["val"].update(s)

    expected = DataFrame(dict(a=[1, 2, 3], b=[1, 2, 3], c=[1, 2, 3], val=[0, 1, 0]))
    expected.set_index(["a", "b", "c"], inplace=True)
    tm.assert_frame_equal(df, expected)

    # GH 3970
    # these are chained assignments as well
    pd.set_option("chained_assignment", None)
    df = DataFrame({"aa": range(5), "bb": [2.2] * 5})
    df["cc"] = 0.0

    ck = [True] * len(df)

    df["bb"].iloc[0] = 0.13

    # TODO: unused
    df_tmp = df.iloc[ck]  # noqa

    df["bb"].iloc[0] = 0.15
    assert df["bb"].iloc[0] == 0.15
    pd.set_option("chained_assignment", "raise")

    # GH 3217
    df = DataFrame(dict(a=[1, 3], b=[np.nan, 2]))
    df["c"] = np.nan
    df["c"].update(pd.Series(["foo"], index=[0]))

    expected = DataFrame(dict(a=[1, 3], b=[np.nan, 2], c=["foo", np.nan]))
    tm.assert_frame_equal(df, expected)


def test_preserve_refs(datetime_series):
    seq = datetime_series[[5, 10, 15]]
    seq[1] = np.NaN
    assert not np.isnan(datetime_series[10])


def test_cast_on_putmask():
    # GH 2746

    # need to upcast
    s = Series([1, 2], index=[1, 2], dtype="int64")
    s[[True, False]] = Series([0], index=[1], dtype="int64")
    expected = Series([0, 2], index=[1, 2], dtype="int64")

    tm.assert_series_equal(s, expected)


def test_type_promote_putmask():
    # GH8387: test that changing types does not break alignment
    ts = Series(np.random.randn(100), index=np.arange(100, 0, -1)).round(5)
    left, mask = ts.copy(), ts > 0
    right = ts[mask].copy().map(str)
    left[mask] = right
    tm.assert_series_equal(left, ts.map(lambda t: str(t) if t > 0 else t))

    s = Series([0, 1, 2, 0])
    mask = s > 0
    s2 = s[mask].map(str)
    s[mask] = s2
    tm.assert_series_equal(s, Series([0, "1", "2", 0]))

    s = Series([0, "foo", "bar", 0])
    mask = Series([False, True, True, False])
    s2 = s[mask]
    s[mask] = s2
    tm.assert_series_equal(s, Series([0, "foo", "bar", 0]))


def test_multilevel_preserve_name():
    index = MultiIndex(
        levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
        codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
        names=["first", "second"],
    )
    s = Series(np.random.randn(len(index)), index=index, name="sth")

    result = s["foo"]
    result2 = s.loc["foo"]
    assert result.name == s.name
    assert result2.name == s.name


def test_setitem_scalar_into_readonly_backing_data():
    # GH14359: test that you cannot mutate a read only buffer

    array = np.zeros(5)
    array.flags.writeable = False  # make the array immutable
    series = Series(array)

    for n in range(len(series)):
        msg = "assignment destination is read-only"
        with pytest.raises(ValueError, match=msg):
            series[n] = 1

        assert array[n] == 0


def test_setitem_slice_into_readonly_backing_data():
    # GH14359: test that you cannot mutate a read only buffer

    array = np.zeros(5)
    array.flags.writeable = False  # make the array immutable
    series = Series(array)

    msg = "assignment destination is read-only"
    with pytest.raises(ValueError, match=msg):
        series[1:3] = 1

    assert not array.any()


"""
miscellaneous methods
"""


def test_pop():
    # GH 6600
    df = DataFrame({"A": 0, "B": np.arange(5, dtype="int64"), "C": 0})
    k = df.iloc[4]

    result = k.pop("B")
    assert result == 4

    expected = Series([0, 0], index=["A", "C"], name=4)
    tm.assert_series_equal(k, expected)


def test_uint_drop(any_int_dtype):
    # see GH18311
    # assigning series.loc[0] = 4 changed series.dtype to int
    series = pd.Series([1, 2, 3], dtype=any_int_dtype)
    series.loc[0] = 4
    expected = pd.Series([4, 2, 3], dtype=any_int_dtype)
    tm.assert_series_equal(series, expected)


def test_getitem_2d_no_warning():
    # https://github.com/pandas-dev/pandas/issues/30867
    # Don't want to support this long-term, but
    # for now ensure that the warning from Index
    # doesn't comes through via Series.__getitem__.
    series = pd.Series([1, 2, 3], index=[1, 2, 3])
    with tm.assert_produces_warning(None):
        series[:, None]


def test_getitem_unrecognized_scalar():
    # GH#32684 a scalar key that is not recognized by lib.is_scalar

    # a series that might be produced via `frame.dtypes`
    ser = pd.Series([1, 2], index=[np.dtype("O"), np.dtype("i8")])

    key = ser.index[1]

    result = ser[key]
    assert result == 2
