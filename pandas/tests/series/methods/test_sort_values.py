import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestSeriesSortValues:
    @pytest.mark.filterwarnings("ignore:The inplace keyword in Series.sort_values")
    def test_sort_values(self, datetime_series):
        # check indexes are reordered corresponding with the values
        ser = pd.Series([3, 2, 4, 1], ["A", "B", "C", "D"])
        expected = pd.Series([1, 2, 3, 4], ["D", "B", "A", "C"])
        result = ser.sort_values()
        tm.assert_series_equal(expected, result)

        ts = datetime_series.copy()
        ts[:5] = np.nan
        vals = ts.values

        result = ts.sort_values()
        assert np.isnan(result[-5:]).all()
        tm.assert_numpy_array_equal(result[:-5].values, np.sort(vals[5:]))

        # na_position
        result = ts.sort_values(na_position="first")
        assert np.isnan(result[:5]).all()
        tm.assert_numpy_array_equal(result[5:].values, np.sort(vals[5:]))

        # something object-type
        ser = pd.Series(["A", "B"], [1, 2])
        # no failure
        ser.sort_values()

        # ascending=False
        ordered = ts.sort_values(ascending=False)
        expected = np.sort(ts.dropna().values)[::-1]
        tm.assert_almost_equal(expected, ordered.dropna().values)
        ordered = ts.sort_values(ascending=False, na_position="first")
        tm.assert_almost_equal(expected, ordered.dropna().values)

        # ascending=[False] should behave the same as ascending=False
        ordered = ts.sort_values(ascending=[False])
        expected = ts.sort_values(ascending=False)
        tm.assert_series_equal(expected, ordered)
        ordered = ts.sort_values(ascending=[False], na_position="first")
        expected = ts.sort_values(ascending=False, na_position="first")
        tm.assert_series_equal(expected, ordered)

        msg = 'For argument "ascending" expected type bool, received type NoneType.'
        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending=None)
        msg = r"Length of ascending \(0\) must be 1 for Series"
        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending=[])
        msg = r"Length of ascending \(3\) must be 1 for Series"
        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending=[1, 2, 3])
        msg = r"Length of ascending \(2\) must be 1 for Series"
        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending=[False, False])
        msg = 'For argument "ascending" expected type bool, received type str.'
        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending="foobar")

        # inplace=True
        ts = datetime_series.copy()
        return_value = ts.sort_values(ascending=False, inplace=True)
        assert return_value is None
        tm.assert_series_equal(ts, datetime_series.sort_values(ascending=False))
        tm.assert_index_equal(
            ts.index, datetime_series.sort_values(ascending=False).index
        )

        # GH#5856/5853
        # Series.sort_values operating on a view
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((10, 4)))
        s = df.iloc[:, 0]
        s.sort_values(inplace=True)
        tm.assert_series_equal(s, df.iloc[:, 0].sort_values())

    def test_sort_values_categorical(self):
        cat = pd.Series(pd.Categorical(["a", "b", "b", "a"], ordered=False))

        # sort in the categories order
        expected = pd.Series(
            pd.Categorical(["a", "a", "b", "b"], ordered=False), index=[0, 3, 1, 2]
        )
        result = cat.sort_values()
        tm.assert_series_equal(result, expected)

        cat = pd.Series(pd.Categorical(["a", "c", "b", "d"], ordered=True))
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        cat = pd.Series(
            pd.Categorical(
                ["a", "c", "b", "d"], categories=["a", "b", "c", "d"], ordered=True
            )
        )
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        res = cat.sort_values(ascending=False)
        exp = np.array(["d", "c", "b", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        raw_cat1 = pd.Categorical(
            ["a", "b", "c", "d"], categories=["a", "b", "c", "d"], ordered=False
        )
        raw_cat2 = pd.Categorical(
            ["a", "b", "c", "d"], categories=["d", "c", "b", "a"], ordered=True
        )
        s = ["a", "b", "c", "d"]
        df = pd.DataFrame(
            {"unsort": raw_cat1, "sort": raw_cat2, "string": s, "values": [1, 2, 3, 4]}
        )

        # Cats must be sorted in a dataframe
        res = df.sort_values(by=["string"], ascending=False)
        exp = np.array(["d", "c", "b", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(res["sort"].values.__array__(), exp)
        assert res["sort"].dtype == "category"

        res = df.sort_values(by=["sort"], ascending=False)
        exp = df.sort_values(by=["string"], ascending=True)
        tm.assert_series_equal(res["values"], exp["values"])
        assert res["sort"].dtype == "category"
        assert res["unsort"].dtype == "category"

        # unordered cat, but we allow this
        df.sort_values(by=["unsort"], ascending=False)

        # multi-columns sort
        # GH#7848
        df = pd.DataFrame(
            {"id": [6, 5, 4, 3, 2, 1], "raw_grade": ["a", "b", "b", "a", "a", "e"]}
        )
        df["grade"] = pd.Categorical(df["raw_grade"], ordered=True)
        df["grade"] = df["grade"].cat.set_categories(["b", "e", "a"])

        # sorts 'grade' according to the order of the categories
        result = df.sort_values(by=["grade"])
        expected = df.iloc[[1, 2, 5, 0, 3, 4]]
        tm.assert_frame_equal(result, expected)

        # multi
        result = df.sort_values(by=["grade", "id"])
        expected = df.iloc[[2, 1, 5, 4, 3, 0]]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.filterwarnings("ignore:The inplace keyword in Series.sort_values")
    @pytest.mark.parametrize("inplace", [True, False])
    @pytest.mark.parametrize(
        "original_list, sorted_list, ignore_index, output_index",
        [
            ([2, 3, 6, 1], [6, 3, 2, 1], True, [0, 1, 2, 3]),
            ([2, 3, 6, 1], [6, 3, 2, 1], False, [2, 1, 0, 3]),
        ],
    )
    def test_sort_values_ignore_index(
        self, inplace, original_list, sorted_list, ignore_index, output_index
    ):
        # GH 30114
        ser = pd.Series(original_list)
        expected = pd.Series(sorted_list, index=output_index)
        kwargs = {"ignore_index": ignore_index, "inplace": inplace}

        if inplace:
            result_ser = ser.copy()
            result_ser.sort_values(ascending=False, **kwargs)
        else:
            result_ser = ser.sort_values(ascending=False, **kwargs)

        tm.assert_series_equal(result_ser, expected)
        tm.assert_series_equal(ser, pd.Series(original_list))

    @pytest.mark.filterwarnings("ignore:The inplace keyword in Series.sort_values")
    def test_sort_values_ignore_index_on_already_sorted(self):
        # GH 65833 - ignore_index had no effect on already sorted Series
        ser = pd.Series([1, 2, 3], index=[2, 3, 4])
        expected = pd.Series([1, 2, 3], index=[0, 1, 2])
        result = ser.sort_values(ignore_index=True)
        tm.assert_series_equal(result, expected)

        ser2 = pd.Series([1, 2, 3], index=[2, 3, 4])
        ser2.sort_values(ignore_index=True, inplace=True)
        tm.assert_series_equal(ser2, expected)

        ser3 = pd.Series([1, 2, 3], index=[2, 3, 4])
        result3 = ser3.sort_values(ignore_index=False)
        tm.assert_series_equal(result3, pd.Series([1, 2, 3], index=[2, 3, 4]))

    def test_mergesort_descending_stability(self):
        # GH 28697
        s = pd.Series([1, 2, 1, 3], ["first", "b", "second", "c"])
        result = s.sort_values(ascending=False, kind="mergesort")
        expected = pd.Series([3, 2, 1, 1], ["c", "b", "first", "second"])
        tm.assert_series_equal(result, expected)

    def test_sort_values_validate_ascending_for_value_error(self):
        # GH41634
        ser = pd.Series([23, 7, 21])

        msg = 'For argument "ascending" expected type bool, received type str.'
        with pytest.raises(ValueError, match=msg):
            ser.sort_values(ascending="False")

    def test_sort_values_validate_ascending_functional(self, ascending):
        # GH41634
        ser = pd.Series([23, 7, 21])
        expected = np.sort(ser.values)

        sorted_ser = ser.sort_values(ascending=ascending)
        if not ascending:
            expected = expected[::-1]

        result = sorted_ser.values
        tm.assert_numpy_array_equal(result, expected)


class TestSeriesSortingKey:
    def test_sort_values_key(self):
        series = pd.Series(np.array(["Hello", "goodbye"]))

        result = series.sort_values(axis=0)
        expected = series
        tm.assert_series_equal(result, expected)

        result = series.sort_values(axis=0, key=lambda x: x.str.lower())
        expected = series[::-1]
        tm.assert_series_equal(result, expected)

    def test_sort_values_key_nan(self):
        series = pd.Series(np.array([0, 5, np.nan, 3, 2, np.nan]))

        result = series.sort_values(axis=0)
        expected = series.iloc[[0, 4, 3, 1, 2, 5]]
        tm.assert_series_equal(result, expected)

        result = series.sort_values(axis=0, key=lambda x: x + 5)
        expected = series.iloc[[0, 4, 3, 1, 2, 5]]
        tm.assert_series_equal(result, expected)

        result = series.sort_values(axis=0, key=lambda x: -x, ascending=False)
        expected = series.iloc[[0, 4, 3, 1, 2, 5]]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "ascending, expected_index", [(True, [0, 2, 1]), (False, [1, 2, 0])]
    )
    def test_sort_values_key_map_unordered_categorical(self, ascending, expected_index):
        # GH#58153
        ser = pd.Series(["March", "Dec", "April"], dtype="category")
        mapping = {"March": 0, "April": 1, "Dec": 3}

        result = ser.sort_values(key=lambda x: x.map(mapping), ascending=ascending)

        tm.assert_series_equal(result, ser.loc[expected_index])

    def test_sort_values_key_map_swaps_categories(self):
        # GH#58153 the new labels are the same set as the old ones
        ser = pd.Series(["a", "b", "c"], dtype="category")
        mapping = {"a": "c", "b": "b", "c": "a"}

        result = ser.sort_values(key=lambda x: x.map(mapping))

        tm.assert_series_equal(result, ser.loc[[2, 1, 0]])

    @pytest.mark.parametrize(
        "na_position, expected_index",
        [("last", [0, 2, 1, 3]), ("first", [3, 0, 2, 1])],
    )
    def test_sort_values_key_map_categorical_with_na(self, na_position, expected_index):
        # GH#58153
        ser = pd.Series(["March", "Dec", "April", None], dtype="category")
        mapping = {"March": 0, "April": 1, "Dec": 3}

        result = ser.sort_values(key=lambda x: x.map(mapping), na_position=na_position)

        tm.assert_series_equal(result, ser.loc[expected_index])

    def test_sort_values_key_map_categorical_unsortable_labels(self):
        # GH#58153 labels that cannot be compared keep the category order
        ser = pd.Series(["b", "a"], dtype="category")

        result = ser.sort_values(key=lambda x: x.map({"a": 1, "b": "z"}))

        tm.assert_series_equal(result, ser.loc[[1, 0]])

    @pytest.mark.parametrize(
        "key",
        [
            lambda x: x.cat.reorder_categories(["c", "a", "b", "d"]),
            lambda x: x.cat.set_categories(["c", "a", "b", "d", "e"]),
            lambda x: x.cat.reorder_categories(["c", "b", "a", "d"]),
        ],
    )
    def test_sort_values_key_keeps_explicit_category_order(self, key):
        # GH#58153 a key that only sets the order must still be respected
        ser = pd.Series(
            pd.Categorical(["a", "b", "c"], categories=["c", "b", "a", "d"])
        )

        result = ser.sort_values(key=key)
        expected = ser.iloc[np.argsort(key(ser).cat.codes.to_numpy(), kind="stable")]

        tm.assert_series_equal(result, expected)

    def test_sort_values_key_moves_unused_category_only(self):
        # GH#58153 codes are unchanged and the values are the same, so the
        # category order is still what the key asks for
        ser = pd.Series(pd.Categorical(["a", "b"], categories=["b", "a", "c", "d"]))

        result = ser.sort_values(
            key=lambda x: x.cat.reorder_categories(["b", "a", "d", "c"])
        )

        tm.assert_series_equal(result, ser.loc[[1, 0]])

    def test_sort_values_key_map_ordered_categorical_keeps_order(self):
        # GH#58153 category order is meaningful for ordered categoricals
        ser = pd.Series(
            pd.Categorical(["low", "HIGH"], categories=["low", "HIGH"], ordered=True)
        )

        result = ser.sort_values(key=lambda x: x.map(str.lower))

        tm.assert_series_equal(result, ser)


def test_sort_values_inplace_depr():
    msg = "The inplace keyword in Series.sort_values is deprecated"

    ser = pd.Series([3, 1, 2])
    ser_orig = ser.copy()
    expected = pd.Series([1, 2, 3], index=[1, 2, 0])

    # does not use keyword, no warning
    with tm.assert_produces_warning(False):
        result = ser.sort_values()
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(ser, ser_orig)

    # uses keyword, set to false, warning
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.sort_values(inplace=False)
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(ser, ser_orig)

    # uses keyword, set to true, warning
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        ser.sort_values(inplace=True)
    tm.assert_series_equal(ser, expected)
