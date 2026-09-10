import numpy as np

import pandas as pd
import pandas._testing as tm


class TestGetLevelValues:
    def test_get_level_values_box_datetime64(self):
        dates = pd.date_range("1/1/2000", periods=4)
        levels = [dates, [0, 1]]
        codes = [[0, 0, 1, 1, 2, 2, 3, 3], [0, 1, 0, 1, 0, 1, 0, 1]]

        index = pd.MultiIndex(levels=levels, codes=codes)

        assert isinstance(index.get_level_values(0)[0], pd.Timestamp)


def test_get_level_values(idx):
    result = idx.get_level_values(0)
    expected = pd.Index(["foo", "foo", "bar", "baz", "qux", "qux"], name="first")
    tm.assert_index_equal(result, expected)
    assert result.name == "first"

    result = idx.get_level_values("first")
    expected = idx.get_level_values(0)
    tm.assert_index_equal(result, expected)

    # GH 10460
    index = pd.MultiIndex(
        levels=[pd.CategoricalIndex(["A", "B"]), pd.CategoricalIndex([1, 2, 3])],
        codes=[np.array([0, 0, 0, 1, 1, 1]), np.array([0, 1, 2, 0, 1, 2])],
    )

    exp = pd.CategoricalIndex(["A", "A", "A", "B", "B", "B"])
    tm.assert_index_equal(index.get_level_values(0), exp)
    exp = pd.CategoricalIndex([1, 2, 3, 1, 2, 3])
    tm.assert_index_equal(index.get_level_values(1), exp)


def test_get_level_values_all_na():
    # GH#17924 when level entirely consists of nan
    arrays = [[np.nan, np.nan, np.nan], ["a", np.nan, 1]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([np.nan, np.nan, np.nan], dtype=np.float64)
    tm.assert_index_equal(result, expected)

    result = index.get_level_values(1)
    expected = pd.Index(["a", np.nan, 1], dtype=object)
    tm.assert_index_equal(result, expected)


def test_get_level_values_int_with_na():
    # GH#17924
    arrays = [["a", "b", "b"], [1, np.nan, 2]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = pd.Index([1, np.nan, 2])
    tm.assert_index_equal(result, expected)

    arrays = [["a", "b", "b"], [np.nan, np.nan, 2]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = pd.Index([np.nan, np.nan, 2])
    tm.assert_index_equal(result, expected)


def test_get_level_values_na():
    arrays = [[np.nan, np.nan, np.nan], ["a", np.nan, 1]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([np.nan, np.nan, np.nan])
    tm.assert_index_equal(result, expected)

    result = index.get_level_values(1)
    expected = pd.Index(["a", np.nan, 1])
    tm.assert_index_equal(result, expected)

    arrays = [["a", "b", "b"], pd.DatetimeIndex([0, 1, pd.NaT])]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = pd.DatetimeIndex([0, 1, pd.NaT])
    tm.assert_index_equal(result, expected)

    arrays = [[], []]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([], dtype=object)
    tm.assert_index_equal(result, expected)


def test_get_level_values_when_periods():
    # GH33131. See also discussion in GH32669.
    # This test can probably be removed when PeriodIndex._engine is removed.
    idx = pd.MultiIndex.from_arrays(
        [pd.PeriodIndex([pd.Period("2019Q1"), pd.Period("2019Q2")], name="b")]
    )
    idx2 = pd.MultiIndex.from_arrays(
        [idx._get_level_values(level) for level in range(idx.nlevels)]
    )
    assert all(x.is_monotonic_increasing for x in idx2.levels)


def test_values_loses_freq_of_underlying_index():
    # GH#49054
    idx = pd.DatetimeIndex(pd.date_range("20200101", periods=3, freq="BME"))
    expected = idx.copy(deep=True)
    idx2 = pd.Index([1, 2, 3])
    midx = pd.MultiIndex(levels=[idx, idx2], codes=[[0, 1, 2], [0, 1, 2]])
    midx.values
    assert idx.freq is not None
    tm.assert_index_equal(idx, expected)


def test_get_level_values_gets_frequency_correctly():
    # GH#57949 GH#58327
    datetime_index = pd.date_range(
        start=pd.to_datetime("1/1/2018"), periods=4, freq="YS"
    )
    other_index = ["A"]
    multi_index = pd.MultiIndex.from_product([datetime_index, other_index])

    assert multi_index.get_level_values(0).freq == datetime_index.freq
