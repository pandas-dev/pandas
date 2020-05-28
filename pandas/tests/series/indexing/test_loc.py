import numpy as np
import pytest

import pandas as pd
from pandas import Series, Timestamp
import pandas._testing as tm


@pytest.mark.parametrize("val,expected", [(2 ** 63 - 1, 3), (2 ** 63, 4)])
def test_loc_uint64(val, expected):
    # see gh-19399
    s = Series({2 ** 63 - 1: 3, 2 ** 63: 4})
    assert s.loc[val] == expected


def test_loc_getitem(string_series, datetime_series):
    inds = string_series.index[[3, 4, 7]]
    tm.assert_series_equal(string_series.loc[inds], string_series.reindex(inds))
    tm.assert_series_equal(string_series.iloc[5::2], string_series[5::2])

    # slice with indices
    d1, d2 = datetime_series.index[[5, 15]]
    result = datetime_series.loc[d1:d2]
    expected = datetime_series.truncate(d1, d2)
    tm.assert_series_equal(result, expected)

    # boolean
    mask = string_series > string_series.median()
    tm.assert_series_equal(string_series.loc[mask], string_series[mask])

    # ask for index value
    assert datetime_series.loc[d1] == datetime_series[d1]
    assert datetime_series.loc[d2] == datetime_series[d2]


def test_loc_getitem_not_monotonic(datetime_series):
    d1, d2 = datetime_series.index[[5, 15]]

    ts2 = datetime_series[::2][[1, 2, 0]]

    msg = r"Timestamp\('2000-01-10 00:00:00'\)"
    with pytest.raises(KeyError, match=msg):
        ts2.loc[d1:d2]
    with pytest.raises(KeyError, match=msg):
        ts2.loc[d1:d2] = 0


def test_loc_getitem_setitem_integer_slice_keyerrors():
    s = Series(np.random.randn(10), index=list(range(0, 20, 2)))

    # this is OK
    cp = s.copy()
    cp.iloc[4:10] = 0
    assert (cp.iloc[4:10] == 0).all()

    # so is this
    cp = s.copy()
    cp.iloc[3:11] = 0
    assert (cp.iloc[3:11] == 0).values.all()

    result = s.iloc[2:6]
    result2 = s.loc[3:11]
    expected = s.reindex([4, 6, 8, 10])

    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result2, expected)

    # non-monotonic, raise KeyError
    s2 = s.iloc[list(range(5)) + list(range(9, 4, -1))]
    with pytest.raises(KeyError, match=r"^3$"):
        s2.loc[3:11]
    with pytest.raises(KeyError, match=r"^3$"):
        s2.loc[3:11] = 0


def test_loc_getitem_iterator(string_series):
    idx = iter(string_series.index[:10])
    result = string_series.loc[idx]
    tm.assert_series_equal(result, string_series[:10])


def test_loc_setitem_boolean(string_series):
    mask = string_series > string_series.median()

    result = string_series.copy()
    result.loc[mask] = 0
    expected = string_series
    expected[mask] = 0
    tm.assert_series_equal(result, expected)


def test_loc_setitem_corner(string_series):
    inds = list(string_series.index[[5, 8, 12]])
    string_series.loc[inds] = 5
    msg = r"\['foo'\] not in index"
    with pytest.raises(KeyError, match=msg):
        string_series.loc[inds + ["foo"]] = 5


def test_basic_setitem_with_labels(datetime_series):
    indices = datetime_series.index[[5, 10, 15]]

    cp = datetime_series.copy()
    exp = datetime_series.copy()
    cp[indices] = 0
    exp.loc[indices] = 0
    tm.assert_series_equal(cp, exp)

    cp = datetime_series.copy()
    exp = datetime_series.copy()
    cp[indices[0] : indices[2]] = 0
    exp.loc[indices[0] : indices[2]] = 0
    tm.assert_series_equal(cp, exp)

    # integer indexes, be careful
    s = Series(np.random.randn(10), index=list(range(0, 20, 2)))
    inds = [0, 4, 6]
    arr_inds = np.array([0, 4, 6])

    cp = s.copy()
    exp = s.copy()
    s[inds] = 0
    s.loc[inds] = 0
    tm.assert_series_equal(cp, exp)

    cp = s.copy()
    exp = s.copy()
    s[arr_inds] = 0
    s.loc[arr_inds] = 0
    tm.assert_series_equal(cp, exp)

    inds_notfound = [0, 4, 5, 6]
    arr_inds_notfound = np.array([0, 4, 5, 6])
    msg = r"\[5\] not in index"
    with pytest.raises(KeyError, match=msg):
        s[inds_notfound] = 0
    with pytest.raises(Exception, match=msg):
        s[arr_inds_notfound] = 0

    # GH12089
    # with tz for values
    s = Series(
        pd.date_range("2011-01-01", periods=3, tz="US/Eastern"), index=["a", "b", "c"]
    )
    s2 = s.copy()
    expected = Timestamp("2011-01-03", tz="US/Eastern")
    s2.loc["a"] = expected
    result = s2.loc["a"]
    assert result == expected

    s2 = s.copy()
    s2.iloc[0] = expected
    result = s2.iloc[0]
    assert result == expected

    s2 = s.copy()
    s2["a"] = expected
    result = s2["a"]
    assert result == expected
