import numpy as np
import pytest

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import Categorical, DataFrame, Series
import pandas._testing as tm


class TestSeriesDtypes:
    def test_dtype(self, datetime_series):

        assert datetime_series.dtype == np.dtype("float64")
        assert datetime_series.dtypes == np.dtype("float64")

    def test_astype_from_categorical(self):
        items = ["a", "b", "c", "a"]
        s = Series(items)
        exp = Series(Categorical(items))
        res = s.astype("category")
        tm.assert_series_equal(res, exp)

        items = [1, 2, 3, 1]
        s = Series(items)
        exp = Series(Categorical(items))
        res = s.astype("category")
        tm.assert_series_equal(res, exp)

        df = DataFrame({"cats": [1, 2, 3, 4, 5, 6], "vals": [1, 2, 3, 4, 5, 6]})
        cats = Categorical([1, 2, 3, 4, 5, 6])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        df = DataFrame(
            {"cats": ["a", "b", "b", "a", "a", "d"], "vals": [1, 2, 3, 4, 5, 6]}
        )
        cats = Categorical(["a", "b", "b", "a", "a", "d"])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        # with keywords
        lst = ["a", "b", "c", "a"]
        s = Series(lst)
        exp = Series(Categorical(lst, ordered=True))
        res = s.astype(CategoricalDtype(None, ordered=True))
        tm.assert_series_equal(res, exp)

        exp = Series(Categorical(lst, categories=list("abcdef"), ordered=True))
        res = s.astype(CategoricalDtype(list("abcdef"), ordered=True))
        tm.assert_series_equal(res, exp)

    def test_astype_categorical_to_other(self):
        cat = Categorical([f"{i} - {i + 499}" for i in range(0, 10000, 500)])
        ser = Series(np.random.RandomState(0).randint(0, 10000, 100)).sort_values()
        ser = pd.cut(ser, range(0, 10500, 500), right=False, labels=cat)

        expected = ser
        tm.assert_series_equal(ser.astype("category"), expected)
        tm.assert_series_equal(ser.astype(CategoricalDtype()), expected)
        msg = r"Cannot cast object dtype to float64"
        with pytest.raises(ValueError, match=msg):
            ser.astype("float64")

        cat = Series(Categorical(["a", "b", "b", "a", "a", "c", "c", "c"]))
        exp = Series(["a", "b", "b", "a", "a", "c", "c", "c"])
        tm.assert_series_equal(cat.astype("str"), exp)
        s2 = Series(Categorical(["1", "2", "3", "4"]))
        exp2 = Series([1, 2, 3, 4]).astype("int")
        tm.assert_series_equal(s2.astype("int"), exp2)

        # object don't sort correctly, so just compare that we have the same
        # values
        def cmp(a, b):
            tm.assert_almost_equal(np.sort(np.unique(a)), np.sort(np.unique(b)))

        expected = Series(np.array(ser.values), name="value_group")
        cmp(ser.astype("object"), expected)
        cmp(ser.astype(np.object_), expected)

        # array conversion
        tm.assert_almost_equal(np.array(ser), np.array(ser.values))

        tm.assert_series_equal(ser.astype("category"), ser)
        tm.assert_series_equal(ser.astype(CategoricalDtype()), ser)

        roundtrip_expected = ser.cat.set_categories(
            ser.cat.categories.sort_values()
        ).cat.remove_unused_categories()
        result = ser.astype("object").astype("category")
        tm.assert_series_equal(result, roundtrip_expected)
        result = ser.astype("object").astype(CategoricalDtype())
        tm.assert_series_equal(result, roundtrip_expected)

    def test_astype_categorical_invalid_conversions(self):
        # invalid conversion (these are NOT a dtype)
        cat = Categorical([f"{i} - {i + 499}" for i in range(0, 10000, 500)])
        ser = Series(np.random.RandomState(0).randint(0, 10000, 100)).sort_values()
        ser = pd.cut(ser, range(0, 10500, 500), right=False, labels=cat)

        msg = (
            "dtype '<class 'pandas.core.arrays.categorical.Categorical'>' "
            "not understood"
        )
        with pytest.raises(TypeError, match=msg):
            ser.astype(Categorical)
        with pytest.raises(TypeError, match=msg):
            ser.astype("object").astype(Categorical)

    def test_categorical_astype_to_int(self, any_int_or_nullable_int_dtype):
        # GH 39402

        df = DataFrame(data={"col1": pd.array([2.0, 1.0, 3.0])})
        df.col1 = df.col1.astype("category")
        df.col1 = df.col1.astype(any_int_or_nullable_int_dtype)
        expected = DataFrame(
            {"col1": pd.array([2, 1, 3], dtype=any_int_or_nullable_int_dtype)}
        )
        tm.assert_frame_equal(df, expected)

    def test_series_to_categorical(self):
        # see gh-16524: test conversion of Series to Categorical
        series = Series(["a", "b", "c"])

        result = Series(series, dtype="category")
        expected = Series(["a", "b", "c"], dtype="category")

        tm.assert_series_equal(result, expected)

    def test_reindex_astype_order_consistency(self):
        # GH 17444
        s = Series([1, 2, 3], index=[2, 0, 1])
        new_index = [0, 1, 2]
        temp_dtype = "category"
        new_dtype = str
        s1 = s.reindex(new_index).astype(temp_dtype).astype(new_dtype)
        s2 = s.astype(temp_dtype).reindex(new_index).astype(new_dtype)
        tm.assert_series_equal(s1, s2)
