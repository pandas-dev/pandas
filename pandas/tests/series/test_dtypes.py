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

        value = np.random.RandomState(0).randint(0, 10000, 100)
        df = DataFrame({"value": value})
        labels = [f"{i} - {i + 499}" for i in range(0, 10000, 500)]
        cat_labels = Categorical(labels, labels)

        df = df.sort_values(by=["value"], ascending=True)
        df["value_group"] = pd.cut(
            df.value, range(0, 10500, 500), right=False, labels=cat_labels
        )

        s = df["value_group"]
        expected = s
        tm.assert_series_equal(s.astype("category"), expected)
        tm.assert_series_equal(s.astype(CategoricalDtype()), expected)
        msg = r"could not convert string to float|invalid literal for float\(\)"
        with pytest.raises(ValueError, match=msg):
            s.astype("float64")

        cat = Series(Categorical(["a", "b", "b", "a", "a", "c", "c", "c"]))
        exp = Series(["a", "b", "b", "a", "a", "c", "c", "c"])
        tm.assert_series_equal(cat.astype("str"), exp)
        s2 = Series(Categorical(["1", "2", "3", "4"]))
        exp2 = Series([1, 2, 3, 4]).astype(int)
        tm.assert_series_equal(s2.astype("int"), exp2)

        # object don't sort correctly, so just compare that we have the same
        # values
        def cmp(a, b):
            tm.assert_almost_equal(np.sort(np.unique(a)), np.sort(np.unique(b)))

        expected = Series(np.array(s.values), name="value_group")
        cmp(s.astype("object"), expected)
        cmp(s.astype(np.object_), expected)

        # array conversion
        tm.assert_almost_equal(np.array(s), np.array(s.values))

        tm.assert_series_equal(s.astype("category"), s)
        tm.assert_series_equal(s.astype(CategoricalDtype()), s)

        roundtrip_expected = s.cat.set_categories(
            s.cat.categories.sort_values()
        ).cat.remove_unused_categories()
        tm.assert_series_equal(
            s.astype("object").astype("category"), roundtrip_expected
        )
        tm.assert_series_equal(
            s.astype("object").astype(CategoricalDtype()), roundtrip_expected
        )

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
