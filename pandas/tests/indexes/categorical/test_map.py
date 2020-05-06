import numpy as np
import pytest

import pandas as pd
from pandas import CategoricalIndex, Index
import pandas._testing as tm


class TestMap:
    @pytest.mark.parametrize(
        "data, categories",
        [
            (list("abcbca"), list("cab")),
            (pd.interval_range(0, 3).repeat(3), pd.interval_range(0, 3)),
        ],
        ids=["string", "interval"],
    )
    def test_map_str(self, data, categories, ordered):
        # GH 31202 - override base class since we want to maintain categorical/ordered
        index = CategoricalIndex(data, categories=categories, ordered=ordered)
        result = index.map(str)
        expected = CategoricalIndex(
            map(str, data), categories=map(str, categories), ordered=ordered
        )
        tm.assert_index_equal(result, expected)

    def test_map(self):
        ci = pd.CategoricalIndex(list("ABABC"), categories=list("CBA"), ordered=True)
        result = ci.map(lambda x: x.lower())
        exp = pd.CategoricalIndex(list("ababc"), categories=list("cba"), ordered=True)
        tm.assert_index_equal(result, exp)

        ci = pd.CategoricalIndex(
            list("ABABC"), categories=list("BAC"), ordered=False, name="XXX"
        )
        result = ci.map(lambda x: x.lower())
        exp = pd.CategoricalIndex(
            list("ababc"), categories=list("bac"), ordered=False, name="XXX"
        )
        tm.assert_index_equal(result, exp)

        # GH 12766: Return an index not an array
        tm.assert_index_equal(
            ci.map(lambda x: 1), Index(np.array([1] * 5, dtype=np.int64), name="XXX")
        )

        # change categories dtype
        ci = pd.CategoricalIndex(list("ABABC"), categories=list("BAC"), ordered=False)

        def f(x):
            return {"A": 10, "B": 20, "C": 30}.get(x)

        result = ci.map(f)
        exp = pd.CategoricalIndex(
            [10, 20, 10, 20, 30], categories=[20, 10, 30], ordered=False
        )
        tm.assert_index_equal(result, exp)

        result = ci.map(pd.Series([10, 20, 30], index=["A", "B", "C"]))
        tm.assert_index_equal(result, exp)

        result = ci.map({"A": 10, "B": 20, "C": 30})
        tm.assert_index_equal(result, exp)

    def test_map_with_categorical_series(self):
        # GH 12756
        a = pd.Index([1, 2, 3, 4])
        b = pd.Series(["even", "odd", "even", "odd"], dtype="category")
        c = pd.Series(["even", "odd", "even", "odd"])

        exp = CategoricalIndex(["odd", "even", "odd", np.nan])
        tm.assert_index_equal(a.map(b), exp)
        exp = pd.Index(["odd", "even", "odd", np.nan])
        tm.assert_index_equal(a.map(c), exp)

    @pytest.mark.parametrize(
        ("data", "f"),
        (
            ([1, 1, np.nan], pd.isna),
            ([1, 2, np.nan], pd.isna),
            ([1, 1, np.nan], {1: False}),
            ([1, 2, np.nan], {1: False, 2: False}),
            ([1, 1, np.nan], pd.Series([False, False])),
            ([1, 2, np.nan], pd.Series([False, False, False])),
        ),
    )
    def test_map_with_nan(self, data, f):  # GH 24241
        values = pd.Categorical(data)
        result = values.map(f)
        if data[1] == 1:
            expected = pd.Categorical([False, False, np.nan])
            tm.assert_categorical_equal(result, expected)
        else:
            expected = pd.Index([False, False, np.nan])
            tm.assert_index_equal(result, expected)
