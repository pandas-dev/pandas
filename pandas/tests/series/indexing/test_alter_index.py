from datetime import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import Categorical, Series, date_range, isna
import pandas._testing as tm


@pytest.mark.parametrize(
    "first_slice,second_slice",
    [
        [[2, None], [None, -5]],
        [[None, 0], [None, -5]],
        [[None, -5], [None, 0]],
        [[None, 0], [None, 0]],
    ],
)
@pytest.mark.parametrize("fill", [None, -1])
def test_align(datetime_series, first_slice, second_slice, join_type, fill):
    a = datetime_series[slice(*first_slice)]
    b = datetime_series[slice(*second_slice)]

    aa, ab = a.align(b, join=join_type, fill_value=fill)

    join_index = a.index.join(b.index, how=join_type)
    if fill is not None:
        diff_a = aa.index.difference(join_index)
        diff_b = ab.index.difference(join_index)
        if len(diff_a) > 0:
            assert (aa.reindex(diff_a) == fill).all()
        if len(diff_b) > 0:
            assert (ab.reindex(diff_b) == fill).all()

    ea = a.reindex(join_index)
    eb = b.reindex(join_index)

    if fill is not None:
        ea = ea.fillna(fill)
        eb = eb.fillna(fill)

    tm.assert_series_equal(aa, ea)
    tm.assert_series_equal(ab, eb)
    assert aa.name == "ts"
    assert ea.name == "ts"
    assert ab.name == "ts"
    assert eb.name == "ts"


@pytest.mark.parametrize(
    "first_slice,second_slice",
    [
        [[2, None], [None, -5]],
        [[None, 0], [None, -5]],
        [[None, -5], [None, 0]],
        [[None, 0], [None, 0]],
    ],
)
@pytest.mark.parametrize("method", ["pad", "bfill"])
@pytest.mark.parametrize("limit", [None, 1])
def test_align_fill_method(
    datetime_series, first_slice, second_slice, join_type, method, limit
):
    a = datetime_series[slice(*first_slice)]
    b = datetime_series[slice(*second_slice)]

    aa, ab = a.align(b, join=join_type, method=method, limit=limit)

    join_index = a.index.join(b.index, how=join_type)
    ea = a.reindex(join_index)
    eb = b.reindex(join_index)

    ea = ea.fillna(method=method, limit=limit)
    eb = eb.fillna(method=method, limit=limit)

    tm.assert_series_equal(aa, ea)
    tm.assert_series_equal(ab, eb)


def test_align_nocopy(datetime_series):
    b = datetime_series[:5].copy()

    # do copy
    a = datetime_series.copy()
    ra, _ = a.align(b, join="left")
    ra[:5] = 5
    assert not (a[:5] == 5).any()

    # do not copy
    a = datetime_series.copy()
    ra, _ = a.align(b, join="left", copy=False)
    ra[:5] = 5
    assert (a[:5] == 5).all()

    # do copy
    a = datetime_series.copy()
    b = datetime_series[:5].copy()
    _, rb = a.align(b, join="right")
    rb[:3] = 5
    assert not (b[:3] == 5).any()

    # do not copy
    a = datetime_series.copy()
    b = datetime_series[:5].copy()
    _, rb = a.align(b, join="right", copy=False)
    rb[:2] = 5
    assert (b[:2] == 5).all()


def test_align_same_index(datetime_series):
    a, b = datetime_series.align(datetime_series, copy=False)
    assert a.index is datetime_series.index
    assert b.index is datetime_series.index

    a, b = datetime_series.align(datetime_series, copy=True)
    assert a.index is not datetime_series.index
    assert b.index is not datetime_series.index


def test_align_multiindex():
    # GH 10665

    midx = pd.MultiIndex.from_product(
        [range(2), range(3), range(2)], names=("a", "b", "c")
    )
    idx = pd.Index(range(2), name="b")
    s1 = pd.Series(np.arange(12, dtype="int64"), index=midx)
    s2 = pd.Series(np.arange(2, dtype="int64"), index=idx)

    # these must be the same results (but flipped)
    res1l, res1r = s1.align(s2, join="left")
    res2l, res2r = s2.align(s1, join="right")

    expl = s1
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = pd.Series([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)

    res1l, res1r = s1.align(s2, join="right")
    res2l, res2r = s2.align(s1, join="left")

    exp_idx = pd.MultiIndex.from_product(
        [range(2), range(2), range(2)], names=("a", "b", "c")
    )
    expl = pd.Series([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = pd.Series([0, 0, 1, 1] * 2, index=exp_idx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)


@pytest.mark.parametrize("method", ["backfill", "bfill", "pad", "ffill", None])
def test_align_method(method):
    # GH31788
    ser = pd.Series(range(3), index=range(3))
    df = pd.DataFrame(0.0, index=range(3), columns=range(3))

    result_ser, result_df = ser.align(df, method=method)
    tm.assert_series_equal(result_ser, ser)
    tm.assert_frame_equal(result_df, df)


def test_reindex(datetime_series, string_series):
    identity = string_series.reindex(string_series.index)

    # __array_interface__ is not defined for older numpies
    # and on some pythons
    try:
        assert np.may_share_memory(string_series.index, identity.index)
    except AttributeError:
        pass

    assert identity.index.is_(string_series.index)
    assert identity.index.identical(string_series.index)

    subIndex = string_series.index[10:20]
    subSeries = string_series.reindex(subIndex)

    for idx, val in subSeries.items():
        assert val == string_series[idx]

    subIndex2 = datetime_series.index[10:20]
    subTS = datetime_series.reindex(subIndex2)

    for idx, val in subTS.items():
        assert val == datetime_series[idx]
    stuffSeries = datetime_series.reindex(subIndex)

    assert np.isnan(stuffSeries).all()

    # This is extremely important for the Cython code to not screw up
    nonContigIndex = datetime_series.index[::2]
    subNonContig = datetime_series.reindex(nonContigIndex)
    for idx, val in subNonContig.items():
        assert val == datetime_series[idx]

    # return a copy the same index here
    result = datetime_series.reindex()
    assert not (result is datetime_series)


def test_reindex_nan():
    ts = Series([2, 3, 5, 7], index=[1, 4, np.nan, 8])

    i, j = [np.nan, 1, np.nan, 8, 4, np.nan], [2, 0, 2, 3, 1, 2]
    tm.assert_series_equal(ts.reindex(i), ts.iloc[j])

    ts.index = ts.index.astype("object")

    # reindex coerces index.dtype to float, loc/iloc doesn't
    tm.assert_series_equal(ts.reindex(i), ts.iloc[j], check_index_type=False)


def test_reindex_series_add_nat():
    rng = date_range("1/1/2000 00:00:00", periods=10, freq="10s")
    series = Series(rng)

    result = series.reindex(range(15))
    assert np.issubdtype(result.dtype, np.dtype("M8[ns]"))

    mask = result.isna()
    assert mask[-5:].all()
    assert not mask[:-5].any()


def test_reindex_with_datetimes():
    rng = date_range("1/1/2000", periods=20)
    ts = Series(np.random.randn(20), index=rng)

    result = ts.reindex(list(ts.index[5:10]))
    expected = ts[5:10]
    tm.assert_series_equal(result, expected)

    result = ts[list(ts.index[5:10])]
    tm.assert_series_equal(result, expected)


def test_reindex_corner(datetime_series):
    # (don't forget to fix this) I think it's fixed
    empty = Series(dtype=object)
    empty.reindex(datetime_series.index, method="pad")  # it works

    # corner case: pad empty series
    reindexed = empty.reindex(datetime_series.index, method="pad")

    # pass non-Index
    reindexed = datetime_series.reindex(list(datetime_series.index))
    tm.assert_series_equal(datetime_series, reindexed)

    # bad fill method
    ts = datetime_series[::2]
    msg = (
        r"Invalid fill method\. Expecting pad \(ffill\), backfill"
        r" \(bfill\) or nearest\. Got foo"
    )
    with pytest.raises(ValueError, match=msg):
        ts.reindex(datetime_series.index, method="foo")


def test_reindex_pad():
    s = Series(np.arange(10), dtype="int64")
    s2 = s[::2]

    reindexed = s2.reindex(s.index, method="pad")
    reindexed2 = s2.reindex(s.index, method="ffill")
    tm.assert_series_equal(reindexed, reindexed2)

    expected = Series([0, 0, 2, 2, 4, 4, 6, 6, 8, 8], index=np.arange(10))
    tm.assert_series_equal(reindexed, expected)

    # GH4604
    s = Series([1, 2, 3, 4, 5], index=["a", "b", "c", "d", "e"])
    new_index = ["a", "g", "c", "f"]
    expected = Series([1, 1, 3, 3], index=new_index)

    # this changes dtype because the ffill happens after
    result = s.reindex(new_index).ffill()
    tm.assert_series_equal(result, expected.astype("float64"))

    result = s.reindex(new_index).ffill(downcast="infer")
    tm.assert_series_equal(result, expected)

    expected = Series([1, 5, 3, 5], index=new_index)
    result = s.reindex(new_index, method="ffill")
    tm.assert_series_equal(result, expected)

    # inference of new dtype
    s = Series([True, False, False, True], index=list("abcd"))
    new_index = "agc"
    result = s.reindex(list(new_index)).ffill()
    expected = Series([True, True, False], index=list(new_index))
    tm.assert_series_equal(result, expected)

    # GH4618 shifted series downcasting
    s = Series(False, index=range(0, 5))
    result = s.shift(1).fillna(method="bfill")
    expected = Series(False, index=range(0, 5))
    tm.assert_series_equal(result, expected)


def test_reindex_nearest():
    s = Series(np.arange(10, dtype="int64"))
    target = [0.1, 0.9, 1.5, 2.0]
    actual = s.reindex(target, method="nearest")
    expected = Series(np.around(target).astype("int64"), target)
    tm.assert_series_equal(expected, actual)

    actual = s.reindex_like(actual, method="nearest")
    tm.assert_series_equal(expected, actual)

    actual = s.reindex_like(actual, method="nearest", tolerance=1)
    tm.assert_series_equal(expected, actual)
    actual = s.reindex_like(actual, method="nearest", tolerance=[1, 2, 3, 4])
    tm.assert_series_equal(expected, actual)

    actual = s.reindex(target, method="nearest", tolerance=0.2)
    expected = Series([0, 1, np.nan, 2], target)
    tm.assert_series_equal(expected, actual)

    actual = s.reindex(target, method="nearest", tolerance=[0.3, 0.01, 0.4, 3])
    expected = Series([0, np.nan, np.nan, 2], target)
    tm.assert_series_equal(expected, actual)


def test_reindex_backfill():
    pass


def test_reindex_int(datetime_series):
    ts = datetime_series[::2]
    int_ts = Series(np.zeros(len(ts), dtype=int), index=ts.index)

    # this should work fine
    reindexed_int = int_ts.reindex(datetime_series.index)

    # if NaNs introduced
    assert reindexed_int.dtype == np.float_

    # NO NaNs introduced
    reindexed_int = int_ts.reindex(int_ts.index[::2])
    assert reindexed_int.dtype == np.int_


def test_reindex_bool(datetime_series):
    # A series other than float, int, string, or object
    ts = datetime_series[::2]
    bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)

    # this should work fine
    reindexed_bool = bool_ts.reindex(datetime_series.index)

    # if NaNs introduced
    assert reindexed_bool.dtype == np.object_

    # NO NaNs introduced
    reindexed_bool = bool_ts.reindex(bool_ts.index[::2])
    assert reindexed_bool.dtype == np.bool_


def test_reindex_bool_pad(datetime_series):
    # fail
    ts = datetime_series[5:]
    bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)
    filled_bool = bool_ts.reindex(datetime_series.index, method="pad")
    assert isna(filled_bool[:5]).all()


def test_reindex_categorical():
    index = date_range("20000101", periods=3)

    # reindexing to an invalid Categorical
    s = Series(["a", "b", "c"], dtype="category")
    result = s.reindex(index)
    expected = Series(
        Categorical(values=[np.nan, np.nan, np.nan], categories=["a", "b", "c"])
    )
    expected.index = index
    tm.assert_series_equal(result, expected)

    # partial reindexing
    expected = Series(Categorical(values=["b", "c"], categories=["a", "b", "c"]))
    expected.index = [1, 2]
    result = s.reindex([1, 2])
    tm.assert_series_equal(result, expected)

    expected = Series(Categorical(values=["c", np.nan], categories=["a", "b", "c"]))
    expected.index = [2, 3]
    result = s.reindex([2, 3])
    tm.assert_series_equal(result, expected)


def test_reindex_like(datetime_series):
    other = datetime_series[::2]
    tm.assert_series_equal(
        datetime_series.reindex(other.index), datetime_series.reindex_like(other)
    )

    # GH 7179
    day1 = datetime(2013, 3, 5)
    day2 = datetime(2013, 5, 5)
    day3 = datetime(2014, 3, 5)

    series1 = Series([5, None, None], [day1, day2, day3])
    series2 = Series([None, None], [day1, day3])

    result = series1.reindex_like(series2, method="pad")
    expected = Series([5, np.nan], index=[day1, day3])
    tm.assert_series_equal(result, expected)


def test_reindex_fill_value():
    # -----------------------------------------------------------
    # floats
    floats = Series([1.0, 2.0, 3.0])
    result = floats.reindex([1, 2, 3])
    expected = Series([2.0, 3.0, np.nan], index=[1, 2, 3])
    tm.assert_series_equal(result, expected)

    result = floats.reindex([1, 2, 3], fill_value=0)
    expected = Series([2.0, 3.0, 0], index=[1, 2, 3])
    tm.assert_series_equal(result, expected)

    # -----------------------------------------------------------
    # ints
    ints = Series([1, 2, 3])

    result = ints.reindex([1, 2, 3])
    expected = Series([2.0, 3.0, np.nan], index=[1, 2, 3])
    tm.assert_series_equal(result, expected)

    # don't upcast
    result = ints.reindex([1, 2, 3], fill_value=0)
    expected = Series([2, 3, 0], index=[1, 2, 3])
    assert issubclass(result.dtype.type, np.integer)
    tm.assert_series_equal(result, expected)

    # -----------------------------------------------------------
    # objects
    objects = Series([1, 2, 3], dtype=object)

    result = objects.reindex([1, 2, 3])
    expected = Series([2, 3, np.nan], index=[1, 2, 3], dtype=object)
    tm.assert_series_equal(result, expected)

    result = objects.reindex([1, 2, 3], fill_value="foo")
    expected = Series([2, 3, "foo"], index=[1, 2, 3], dtype=object)
    tm.assert_series_equal(result, expected)

    # ------------------------------------------------------------
    # bools
    bools = Series([True, False, True])

    result = bools.reindex([1, 2, 3])
    expected = Series([False, True, np.nan], index=[1, 2, 3], dtype=object)
    tm.assert_series_equal(result, expected)

    result = bools.reindex([1, 2, 3], fill_value=False)
    expected = Series([False, True, False], index=[1, 2, 3])
    tm.assert_series_equal(result, expected)


def test_reindex_datetimeindexes_tz_naive_and_aware():
    # GH 8306
    idx = date_range("20131101", tz="America/Chicago", periods=7)
    newidx = date_range("20131103", periods=10, freq="H")
    s = Series(range(7), index=idx)
    with pytest.raises(TypeError):
        s.reindex(newidx, method="ffill")


def test_reindex_empty_series_tz_dtype():
    # GH 20869
    result = Series(dtype="datetime64[ns, UTC]").reindex([0, 1])
    expected = Series([pd.NaT] * 2, dtype="datetime64[ns, UTC]")
    tm.assert_equal(result, expected)


def test_rename():
    # GH 17407
    s = Series(range(1, 6), index=pd.Index(range(2, 7), name="IntIndex"))
    result = s.rename(str)
    expected = s.rename(lambda i: str(i))
    tm.assert_series_equal(result, expected)

    assert result.name == expected.name


@pytest.mark.parametrize(
    "data, index, drop_labels, axis, expected_data, expected_index",
    [
        # Unique Index
        ([1, 2], ["one", "two"], ["two"], 0, [1], ["one"]),
        ([1, 2], ["one", "two"], ["two"], "rows", [1], ["one"]),
        ([1, 1, 2], ["one", "two", "one"], ["two"], 0, [1, 2], ["one", "one"]),
        # GH 5248 Non-Unique Index
        ([1, 1, 2], ["one", "two", "one"], "two", 0, [1, 2], ["one", "one"]),
        ([1, 1, 2], ["one", "two", "one"], ["one"], 0, [1], ["two"]),
        ([1, 1, 2], ["one", "two", "one"], "one", 0, [1], ["two"]),
    ],
)
def test_drop_unique_and_non_unique_index(
    data, index, axis, drop_labels, expected_data, expected_index
):

    s = Series(data=data, index=index)
    result = s.drop(drop_labels, axis=axis)
    expected = Series(data=expected_data, index=expected_index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "data, index, drop_labels, axis, error_type, error_desc",
    [
        # single string/tuple-like
        (range(3), list("abc"), "bc", 0, KeyError, "not found in axis"),
        # bad axis
        (range(3), list("abc"), ("a",), 0, KeyError, "not found in axis"),
        (range(3), list("abc"), "one", "columns", ValueError, "No axis named columns"),
    ],
)
def test_drop_exception_raised(data, index, drop_labels, axis, error_type, error_desc):

    with pytest.raises(error_type, match=error_desc):
        Series(data, index=index).drop(drop_labels, axis=axis)


def test_drop_with_ignore_errors():
    # errors='ignore'
    s = Series(range(3), index=list("abc"))
    result = s.drop("bc", errors="ignore")
    tm.assert_series_equal(result, s)
    result = s.drop(["a", "d"], errors="ignore")
    expected = s.iloc[1:]
    tm.assert_series_equal(result, expected)

    # GH 8522
    s = Series([2, 3], index=[True, False])
    assert s.index.is_object()
    result = s.drop(True)
    expected = Series([3], index=[False])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("index", [[1, 2, 3], [1, 1, 3]])
@pytest.mark.parametrize("drop_labels", [[], [1], [3]])
def test_drop_empty_list(index, drop_labels):
    # GH 21494
    expected_index = [i for i in index if i not in drop_labels]
    series = pd.Series(index=index, dtype=object).drop(drop_labels)
    expected = pd.Series(index=expected_index, dtype=object)
    tm.assert_series_equal(series, expected)


@pytest.mark.parametrize(
    "data, index, drop_labels",
    [
        (None, [1, 2, 3], [1, 4]),
        (None, [1, 2, 2], [1, 4]),
        ([2, 3], [0, 1], [False, True]),
    ],
)
def test_drop_non_empty_list(data, index, drop_labels):
    # GH 21494 and GH 16877
    with pytest.raises(KeyError, match="not found in axis"):
        dtype = object if data is None else None
        pd.Series(data=data, index=index, dtype=dtype).drop(drop_labels)
