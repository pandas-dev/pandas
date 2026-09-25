from datetime import datetime

import numpy as np

from pandas.errors import Pandas4Warning

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
import pandas._testing as tm


class TestCategoricalConcat:
    def test_categorical_concat(self, sort):
        # See GH 10177
        df1 = pd.DataFrame(
            np.arange(18, dtype="int64").reshape(6, 3), columns=["a", "b", "c"]
        )

        df2 = pd.DataFrame(
            np.arange(14, dtype="int64").reshape(7, 2), columns=["a", "c"]
        )

        cat_values = ["one", "one", "two", "one", "two", "two", "one"]
        df2["h"] = pd.Series(pd.Categorical(cat_values))

        res = pd.concat((df1, df2), axis=0, ignore_index=True, sort=sort)
        exp = pd.DataFrame(
            {
                "a": [0, 3, 6, 9, 12, 15, 0, 2, 4, 6, 8, 10, 12],
                "b": [
                    1,
                    4,
                    7,
                    10,
                    13,
                    16,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "c": [2, 5, 8, 11, 14, 17, 1, 3, 5, 7, 9, 11, 13],
                "h": [None] * 6 + cat_values,
            }
        )
        exp["h"] = exp["h"].astype(df2["h"].dtype)
        tm.assert_frame_equal(res, exp)

    def test_categorical_concat_dtypes(self, using_infer_string):
        # GH8143
        index = ["cat", "obj", "num"]
        cat = pd.Categorical(["a", "b", "c"])
        obj = pd.Series(["a", "b", "c"])
        num = pd.Series([1, 2, 3])
        df = pd.concat([pd.Series(cat), obj, num], axis=1, keys=index)

        result = df.dtypes == (object if not using_infer_string else "str")
        expected = pd.Series([False, True, False], index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == "int64"
        expected = pd.Series([False, False, True], index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == "category"
        expected = pd.Series([True, False, False], index=index)
        tm.assert_series_equal(result, expected)

    def test_concat_categoricalindex(self):
        # GH 16111, categories that aren't lexsorted
        categories = [9, 0, 1, 2, 3]

        a = pd.Series(1, index=pd.CategoricalIndex([9, 0], categories=categories))
        b = pd.Series(2, index=pd.CategoricalIndex([0, 1], categories=categories))
        c = pd.Series(3, index=pd.CategoricalIndex([1, 2], categories=categories))

        result = pd.concat([a, b, c], axis=1)

        exp_idx = pd.CategoricalIndex([9, 0, 1, 2], categories=categories)
        exp = pd.DataFrame(
            {
                0: [1, 1, np.nan, np.nan],
                1: [np.nan, 2, 2, np.nan],
                2: [np.nan, np.nan, 3, 3],
            },
            columns=[0, 1, 2],
            index=exp_idx,
        )
        tm.assert_frame_equal(result, exp)

    def test_categorical_concat_preserve(self):
        # GH 8641  series concat not preserving category dtype
        # GH 13524 can concat different categories
        s = pd.Series(list("abc"), dtype="category")
        s2 = pd.Series(list("abd"), dtype="category")

        exp = pd.Series(list("abcabd"))
        res = pd.concat([s, s2], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = pd.Series(list("abcabc"), dtype="category")
        res = pd.concat([s, s], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = pd.Series(list("abcabc"), index=[0, 1, 2, 0, 1, 2], dtype="category")
        res = pd.concat([s, s])
        tm.assert_series_equal(res, exp)

        a = pd.Series(np.arange(6, dtype="int64"))
        b = pd.Series(list("aabbca"))

        df2 = pd.DataFrame({"A": a, "B": b.astype(CategoricalDtype(list("cab")))})
        res = pd.concat([df2, df2])
        exp = pd.DataFrame(
            {
                "A": pd.concat([a, a]),
                "B": pd.concat([b, b]).astype(CategoricalDtype(list("cab"))),
            }
        )
        tm.assert_frame_equal(res, exp)

    def test_categorical_index_preserver(self):
        a = pd.Series(np.arange(6, dtype="int64"))
        b = pd.Series(list("aabbca"))

        df2 = pd.DataFrame(
            {"A": a, "B": b.astype(CategoricalDtype(list("cab")))}
        ).set_index("B")
        result = pd.concat([df2, df2])
        expected = pd.DataFrame(
            {
                "A": pd.concat([a, a]),
                "B": pd.concat([b, b]).astype(CategoricalDtype(list("cab"))),
            }
        ).set_index("B")
        tm.assert_frame_equal(result, expected)

        # wrong categories -> uses concat_compat, which casts to object
        msg = "Constructing a Categorical with a dtype and values containing"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df3 = pd.DataFrame(
                {"A": a, "B": pd.Categorical(b, categories=list("abe"))}
            ).set_index("B")
        result = pd.concat([df2, df3])
        expected = pd.concat(
            [
                df2.set_axis(df2.index.astype(object), axis=0),
                df3.set_axis(df3.index.astype(object), axis=0),
            ]
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_categorical_tz(self):
        # GH-23816
        a = pd.Series(pd.date_range("2017-01-01", periods=2, tz="US/Pacific"))
        b = pd.Series(["a", "b"], dtype="category")
        result = pd.concat([a, b], ignore_index=True)
        expected = pd.Series(
            [
                pd.Timestamp("2017-01-01", tz="US/Pacific"),
                pd.Timestamp("2017-01-02", tz="US/Pacific"),
                "a",
                "b",
            ]
        )
        tm.assert_series_equal(result, expected)

    def test_concat_categorical_datetime(self):
        # GH-39443
        df1 = pd.DataFrame(
            {"x": pd.Series(datetime(2021, 1, 1), index=[0], dtype="category")}
        )
        df2 = pd.DataFrame(
            {"x": pd.Series(datetime(2021, 1, 2), index=[1], dtype="category")}
        )

        result = pd.concat([df1, df2])
        expected = pd.DataFrame(
            {
                "x": pd.Series(
                    [datetime(2021, 1, 1), datetime(2021, 1, 2)],
                    index=[0, 1],
                )
            }
        )

        tm.assert_equal(result, expected)

    def test_concat_categorical_unchanged(self):
        # GH-12007
        # test fix for when concat on categorical and float
        # coerces dtype categorical -> float
        df = pd.DataFrame(pd.Series(["a", "b", "c"], dtype="category", name="A"))
        ser = pd.Series([0, 1, 2], index=[0, 1, 3], name="B")
        result = pd.concat([df, ser], axis=1)
        expected = pd.DataFrame(
            {
                "A": pd.Series(["a", "b", "c", np.nan], dtype="category"),
                "B": pd.Series([0, 1, np.nan, 2], dtype="float"),
            }
        )
        tm.assert_equal(result, expected)

    def test_categorical_concat_gh7864(self):
        # GH 7864
        # make sure ordering is preserved
        df = pd.DataFrame({"id": [1, 2, 3, 4, 5, 6], "raw_grade": list("abbaae")})
        df["grade"] = pd.Categorical(df["raw_grade"])
        df["grade"].cat.set_categories(["e", "a", "b"])

        df1 = df[0:3]
        df2 = df[3:]

        tm.assert_index_equal(df["grade"].cat.categories, df1["grade"].cat.categories)
        tm.assert_index_equal(df["grade"].cat.categories, df2["grade"].cat.categories)

        dfx = pd.concat([df1, df2])
        tm.assert_index_equal(df["grade"].cat.categories, dfx["grade"].cat.categories)

    def test_categorical_index_upcast(self):
        # GH 17629
        # test upcasting to object when concatenating on categorical indexes
        # with non-identical categories

        a = pd.DataFrame({"foo": [1, 2]}, index=pd.Categorical(["foo", "bar"]))
        b = pd.DataFrame({"foo": [4, 3]}, index=pd.Categorical(["baz", "bar"]))

        res = pd.concat([a, b])
        exp = pd.DataFrame({"foo": [1, 2, 4, 3]}, index=["foo", "bar", "baz", "bar"])

        tm.assert_equal(res, exp)

        a = pd.Series([1, 2], index=pd.Categorical(["foo", "bar"]))
        b = pd.Series([4, 3], index=pd.Categorical(["baz", "bar"]))

        res = pd.concat([a, b])
        exp = pd.Series([1, 2, 4, 3], index=["foo", "bar", "baz", "bar"])

        tm.assert_equal(res, exp)

    def test_categorical_missing_from_one_frame(self):
        # GH 25412
        df1 = pd.DataFrame({"f1": [1, 2, 3]})
        df2 = pd.DataFrame(
            {"f1": [2, 3, 1], "f2": pd.Series([4, 4, 4]).astype("category")}
        )
        result = pd.concat([df1, df2], sort=True)
        dtype = CategoricalDtype([4])
        expected = pd.DataFrame(
            {
                "f1": [1, 2, 3, 2, 3, 1],
                "f2": pd.Categorical.from_codes([-1, -1, -1, 0, 0, 0], dtype=dtype),
            },
            index=[0, 1, 2, 0, 1, 2],
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_categorical_same_categories_different_order(self):
        # https://github.com/pandas-dev/pandas/issues/24845

        c1 = pd.CategoricalIndex(["a", "a"], categories=["a", "b"], ordered=False)
        c2 = pd.CategoricalIndex(["b", "b"], categories=["b", "a"], ordered=False)
        c3 = pd.CategoricalIndex(
            ["a", "a", "b", "b"], categories=["a", "b"], ordered=False
        )

        df1 = pd.DataFrame({"A": [1, 2]}, index=c1)
        df2 = pd.DataFrame({"A": [3, 4]}, index=c2)

        result = pd.concat((df1, df2))
        expected = pd.DataFrame({"A": [1, 2, 3, 4]}, index=c3)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_categorical_index_preserves_dtype(self):
        # GH#14016 keys given as a CategoricalIndex should stay categorical
        #  in the resulting MultiIndex level, not be cast to object
        cidx = pd.CategoricalIndex(["y", "x"], categories=list("xyz"), ordered=True)
        df = pd.DataFrame([[10, 11, 12]])
        result = pd.concat([df, df], keys=cidx).index.levels[0]
        tm.assert_index_equal(result, cidx)


# ---------------------------------------------------------------------
# pd.concat with union_categories=True (GH#14177)


def test_union_categories_series_different_categories():
    # GH#14177
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    s2 = pd.Series(pd.Categorical(["b", "c"], categories=["b", "c"]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "b", "c"], categories=["a", "b", "c"]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_series_same_categories():
    # Same categories should still work with union_categories=True
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    s2 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "a", "b"], categories=["a", "b"]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_dataframe_different_categories():
    # GH#14177
    df1 = pd.DataFrame({"x": pd.Categorical(["a", "b"], categories=["a", "b"])})
    df2 = pd.DataFrame({"x": pd.Categorical(["b", "c"], categories=["b", "c"])})
    result = pd.concat([df1, df2], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {"x": pd.Categorical(["a", "b", "b", "c"], categories=["a", "b", "c"])},
    )
    tm.assert_frame_equal(result, expected)


def test_union_categories_default_false_preserves_existing_behavior():
    # union_categories=False (default) should not preserve categorical
    # dtype when categories differ
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    s2 = pd.Series(pd.Categorical(["b", "c"], categories=["b", "c"]))
    result = pd.concat([s1, s2], ignore_index=True)
    expected = pd.Series(["a", "b", "b", "c"])
    tm.assert_series_equal(result, expected)


def test_union_categories_mixed_categorical_noncategorical():
    # When mixing categorical and non-categorical, keyword has no effect
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    s2 = pd.Series(["b", "c"])
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    # Should fall through to existing behavior
    expected = pd.Series(["a", "b", "b", "c"])
    tm.assert_series_equal(result, expected)


def test_union_categories_ordered_same_order_subset_categories():
    # Ordered categoricals with compatible categories
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"], ordered=True))
    s2 = pd.Series(pd.Categorical(["a"], categories=["a", "b"], ordered=True))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "a"], categories=["a", "b"], ordered=True),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_ordered_incompatible_categories_drops_order():
    # Ordered categoricals with different categories union to an unordered
    # result rather than raising (GH#14177, dev-call decision)
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"], ordered=True))
    s2 = pd.Series(pd.Categorical(["b", "c"], categories=["b", "c"], ordered=True))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "b", "c"], categories=["a", "b", "c"]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_mixed_ordered_unordered_drops_order():
    # Mixed ordered/unordered with the same categories: dtypes differ, so
    # the result is unordered rather than raising
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"], ordered=True))
    s2 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "a", "b"], categories=["a", "b"]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_ordered_same_categories_different_order_drops_order():
    # Same category set in a different order means the dtypes differ, so
    # the result is the unordered union rather than an ordered categorical
    s1 = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"], ordered=True))
    s2 = pd.Series(pd.Categorical(["a", "b"], categories=["b", "a"], ordered=True))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(["a", "b", "a", "b"], categories=["a", "b"]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_incompatible_category_dtypes():
    # Categories with different dtypes are cast to a common dtype instead of
    # falling back to a non-categorical result
    s1 = pd.Series(pd.Categorical(["a", "b"]))
    s2 = pd.Series(pd.Categorical([1, 2]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(
            np.array(["a", "b", 1, 2], dtype=object),
            categories=np.array(["a", "b", 1, 2], dtype=object),
        )
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_numeric_category_dtypes():
    # int and float categories are unioned under the common float64 dtype
    s1 = pd.Series(pd.Categorical([1, 2]))
    s2 = pd.Series(pd.Categorical([2.5, 3.5]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical([1.0, 2.0, 2.5, 3.5], categories=[1.0, 2.0, 2.5, 3.5])
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_lossy_category_cast():
    # int64 categories that are not exactly representable as float64 collapse
    # under the common dtype, matching what a non-categorical concat gives
    s1 = pd.Series(pd.Categorical([2**53, 2**53 + 1]))
    s2 = pd.Series(pd.Categorical([2.5]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(
            [float(2**53), float(2**53), 2.5], categories=[float(2**53), 2.5]
        )
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_empty_categorical():
    # GH#14177 an empty Categorical has object-dtype categories, which must not
    # drag the result off of categorical dtype
    ser = pd.Series(pd.Categorical(["a", "b"]))
    result = pd.concat(
        [ser, pd.Series(dtype="category")], ignore_index=True, union_categories=True
    )
    expected = pd.Series(pd.Categorical(["a", "b"], categories=["a", "b"]))
    tm.assert_series_equal(result, expected)


def test_union_categories_all_nan_categorical():
    # an all-NaN Categorical has no categories, so it contributes nothing but
    # the NaN itself
    ser = pd.Series(pd.Categorical(["a", "b"]))
    result = pd.concat(
        [ser, pd.Series(pd.Categorical([np.nan]))],
        ignore_index=True,
        union_categories=True,
    )
    expected = pd.Series(pd.Categorical(["a", "b", np.nan], categories=["a", "b"]))
    tm.assert_series_equal(result, expected)


def test_union_categories_dataframe_column_missing_from_one_frame():
    # Column present in only one frame still keeps categorical dtype,
    # with NaN for the missing rows
    df1 = pd.DataFrame({"x": pd.Categorical(["a", "b"]), "y": [1, 2]})
    df2 = pd.DataFrame({"y": [3]})
    result = pd.concat([df1, df2], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "x": pd.Categorical(["a", "b", np.nan], categories=["a", "b"]),
            "y": [1, 2, 3],
        }
    )
    tm.assert_frame_equal(result, expected)


def test_union_categories_different_categories_column_missing_from_one_frame():
    # Differing categories AND a frame missing the column: the all-NA filler
    # must not degrade the result to object
    df1 = pd.DataFrame({"x": pd.Categorical(["a", "b"]), "y": [1, 2]})
    df2 = pd.DataFrame({"x": pd.Categorical(["b", "c"]), "y": [3, 4]})
    df3 = pd.DataFrame({"y": [5]})
    result = pd.concat([df1, df2, df3], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "x": pd.Categorical(
                ["a", "b", "b", "c", np.nan], categories=["a", "b", "c"]
            ),
            "y": [1, 2, 3, 4, 5],
        }
    )
    tm.assert_frame_equal(result, expected)

    # missing-column frame first
    result = pd.concat([df3, df1, df2], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "y": [5, 1, 2, 3, 4],
            "x": pd.Categorical(
                [np.nan, "a", "b", "b", "c"], categories=["a", "b", "c"]
            ),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_union_categories_ordered_column_missing_from_one_frame():
    # Identical ordered dtypes stay ordered even when another frame lacks
    # the column
    df1 = pd.DataFrame(
        {"x": pd.Categorical(["a"], categories=["a", "b"], ordered=True), "y": [1]}
    )
    df2 = pd.DataFrame(
        {"x": pd.Categorical(["b"], categories=["a", "b"], ordered=True), "y": [2]}
    )
    df3 = pd.DataFrame({"y": [3]})
    result = pd.concat([df1, df2, df3], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "x": pd.Categorical(
                ["a", "b", np.nan], categories=["a", "b"], ordered=True
            ),
            "y": [1, 2, 3],
        }
    )
    tm.assert_frame_equal(result, expected)


def test_union_categories_incompatible_category_dtypes_column_missing():
    # Categories with differing dtypes are unioned also when a frame lacks
    # the column
    df1 = pd.DataFrame({"x": pd.Categorical(["a"]), "y": [1]})
    df2 = pd.DataFrame({"x": pd.Categorical([2]), "y": [2]})
    df3 = pd.DataFrame({"y": [3]})
    result = pd.concat([df1, df2, df3], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "x": pd.Categorical(
                np.array(["a", 2, np.nan], dtype=object),
                categories=np.array(["a", 2], dtype=object),
            ),
            "y": [1, 2, 3],
        }
    )
    tm.assert_frame_equal(result, expected)


def test_union_categories_integer_categories():
    s1 = pd.Series(pd.Categorical([1, 2], categories=[1, 2]))
    s2 = pd.Series(pd.Categorical([2, 3], categories=[2, 3]))
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical([1, 2, 2, 3], categories=[1, 2, 3]),
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_datetime_categories():
    cat1 = pd.Categorical(
        pd.to_datetime(["2020-01-01", "2020-01-02"]),
        categories=pd.to_datetime(["2020-01-01", "2020-01-02"]),
    )
    cat2 = pd.Categorical(
        pd.to_datetime(["2020-01-02", "2020-01-03"]),
        categories=pd.to_datetime(["2020-01-02", "2020-01-03"]),
    )
    s1 = pd.Series(cat1)
    s2 = pd.Series(cat2)
    result = pd.concat([s1, s2], ignore_index=True, union_categories=True)
    expected = pd.Series(
        pd.Categorical(
            pd.to_datetime(["2020-01-01", "2020-01-02", "2020-01-02", "2020-01-03"]),
            categories=pd.to_datetime(["2020-01-01", "2020-01-02", "2020-01-03"]),
        )
    )
    tm.assert_series_equal(result, expected)


def test_union_categories_dataframe_multiple_categorical_columns():
    # Multiple columns, each categorical with different categories
    df1 = pd.DataFrame(
        {
            "x": pd.Categorical(["a", "b"], categories=["a", "b"]),
            "y": pd.Categorical([1, 2], categories=[1, 2]),
        }
    )
    df2 = pd.DataFrame(
        {
            "x": pd.Categorical(["c"], categories=["b", "c"]),
            "y": pd.Categorical([3], categories=[2, 3]),
        }
    )
    result = pd.concat([df1, df2], ignore_index=True, union_categories=True)
    expected = pd.DataFrame(
        {
            "x": pd.Categorical(["a", "b", "c"], categories=["a", "b", "c"]),
            "y": pd.Categorical([1, 2, 3], categories=[1, 2, 3]),
        }
    )
    tm.assert_frame_equal(result, expected)
