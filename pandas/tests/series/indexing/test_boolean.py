import numpy as np
import pytest

from pandas import Index, Series
import pandas._testing as tm
from pandas.core.indexing import IndexingError

from pandas.tseries.offsets import BDay


def test_getitem_boolean(string_series):
    s = string_series
    mask = s > s.median()

    # passing list is OK
    result = s[list(mask)]
    expected = s[mask]
    tm.assert_series_equal(result, expected)
    tm.assert_index_equal(result.index, s.index[mask])


def test_getitem_boolean_empty():
    s = Series([], dtype=np.int64)
    s.index.name = "index_name"
    s = s[s.isna()]
    assert s.index.name == "index_name"
    assert s.dtype == np.int64

    # GH5877
    # indexing with empty series
    s = Series(["A", "B"])
    expected = Series(np.nan, index=["C"], dtype=object)
    result = s[Series(["C"], dtype=object)]
    tm.assert_series_equal(result, expected)

    s = Series(["A", "B"])
    expected = Series(dtype=object, index=Index([], dtype="int64"))
    result = s[Series([], dtype=object)]
    tm.assert_series_equal(result, expected)

    # invalid because of the boolean indexer
    # that's empty or not-aligned
    msg = (
        r"Unalignable boolean Series provided as indexer \(index of "
        r"the boolean Series and of the indexed object do not match"
    )
    with pytest.raises(IndexingError, match=msg):
        s[Series([], dtype=bool)]

    with pytest.raises(IndexingError, match=msg):
        s[Series([True], dtype=bool)]


def test_getitem_boolean_object(string_series):
    # using column from DataFrame

    s = string_series
    mask = s > s.median()
    omask = mask.astype(object)

    # getitem
    result = s[omask]
    expected = s[mask]
    tm.assert_series_equal(result, expected)

    # setitem
    s2 = s.copy()
    cop = s.copy()
    cop[omask] = 5
    s2[mask] = 5
    tm.assert_series_equal(cop, s2)

    # nans raise exception
    omask[5:10] = np.nan
    msg = "Cannot mask with non-boolean array containing NA / NaN values"
    with pytest.raises(ValueError, match=msg):
        s[omask]
    with pytest.raises(ValueError, match=msg):
        s[omask] = 5


def test_getitem_setitem_boolean_corner(datetime_series):
    ts = datetime_series
    mask_shifted = ts.shift(1, freq=BDay()) > ts.median()

    # these used to raise...??

    msg = (
        r"Unalignable boolean Series provided as indexer \(index of "
        r"the boolean Series and of the indexed object do not match"
    )
    with pytest.raises(IndexingError, match=msg):
        ts[mask_shifted]
    with pytest.raises(IndexingError, match=msg):
        ts[mask_shifted] = 1

    with pytest.raises(IndexingError, match=msg):
        ts.loc[mask_shifted]
    with pytest.raises(IndexingError, match=msg):
        ts.loc[mask_shifted] = 1


def test_setitem_boolean(string_series):
    mask = string_series > string_series.median()

    # similar indexed series
    result = string_series.copy()
    result[mask] = string_series * 2
    expected = string_series * 2
    tm.assert_series_equal(result[mask], expected[mask])

    # needs alignment
    result = string_series.copy()
    result[mask] = (string_series * 2)[0:5]
    expected = (string_series * 2)[0:5].reindex_like(string_series)
    expected[-mask] = string_series[mask]
    tm.assert_series_equal(result[mask], expected[mask])


def test_get_set_boolean_different_order(string_series):
    ordered = string_series.sort_values()

    # setting
    copy = string_series.copy()
    copy[ordered > 0] = 0

    expected = string_series.copy()
    expected[expected > 0] = 0

    tm.assert_series_equal(copy, expected)

    # getting
    sel = string_series[ordered > 0]
    exp = string_series[string_series > 0]
    tm.assert_series_equal(sel, exp)
