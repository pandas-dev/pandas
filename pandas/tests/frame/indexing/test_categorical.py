import numpy as np
import pytest

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import Categorical, DataFrame, Index, Series
import pandas._testing as tm

msg1 = "Cannot setitem on a Categorical with a new category, set the categories first"
msg2 = "Cannot set a Categorical with another, without identical categories"


class TestDataFrameIndexingCategorical:
    def test_assignment(self):
        # assignment
        df = DataFrame(
            {"value": np.array(np.random.randint(0, 10000, 100), dtype="int32")}
        )
        labels = Categorical([f"{i} - {i + 499}" for i in range(0, 10000, 500)])

        df = df.sort_values(by=["value"], ascending=True)
        s = pd.cut(df.value, range(0, 10500, 500), right=False, labels=labels)
        d = s.values
        df["D"] = d
        str(df)

        result = df.dtypes
        expected = Series(
            [np.dtype("int32"), CategoricalDtype(categories=labels, ordered=False)],
            index=["value", "D"],
        )
        tm.assert_series_equal(result, expected)

        df["E"] = s
        str(df)

        result = df.dtypes
        expected = Series(
            [
                np.dtype("int32"),
                CategoricalDtype(categories=labels, ordered=False),
                CategoricalDtype(categories=labels, ordered=False),
            ],
            index=["value", "D", "E"],
        )
        tm.assert_series_equal(result, expected)

        result1 = df["D"]
        result2 = df["E"]
        tm.assert_categorical_equal(result1._mgr._block.values, d)

        # sorting
        s.name = "E"
        tm.assert_series_equal(result2.sort_index(), s.sort_index())

        cat = Categorical([1, 2, 3, 10], categories=[1, 2, 3, 4, 10])
        df = DataFrame(Series(cat))

    @pytest.fixture
    def orig(self):
        cats = Categorical(["a", "a", "a", "a", "a", "a", "a"], categories=["a", "b"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n"])
        values = [1, 1, 1, 1, 1, 1, 1]
        orig = DataFrame({"cats": cats, "values": values}, index=idx)
        return orig

    @pytest.fixture
    def exp_single_row(self):
        # The expected values if we change a single row
        cats1 = Categorical(["a", "a", "b", "a", "a", "a", "a"], categories=["a", "b"])
        idx1 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values1 = [1, 1, 2, 1, 1, 1, 1]
        exp_single_row = DataFrame({"cats": cats1, "values": values1}, index=idx1)
        return exp_single_row

    @pytest.fixture
    def exp_multi_row(self):
        # assign multiple rows (mixed values) (-> array) -> exp_multi_row
        # changed multiple rows
        cats2 = Categorical(["a", "a", "b", "b", "a", "a", "a"], categories=["a", "b"])
        idx2 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values2 = [1, 1, 2, 2, 1, 1, 1]
        exp_multi_row = DataFrame({"cats": cats2, "values": values2}, index=idx2)
        return exp_multi_row

    @pytest.fixture
    def exp_parts_cats_col(self):
        # changed part of the cats column
        cats3 = Categorical(["a", "a", "b", "b", "a", "a", "a"], categories=["a", "b"])
        idx3 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values3 = [1, 1, 1, 1, 1, 1, 1]
        exp_parts_cats_col = DataFrame({"cats": cats3, "values": values3}, index=idx3)
        return exp_parts_cats_col

    @pytest.fixture
    def exp_single_cats_value(self):
        # changed single value in cats col
        cats4 = Categorical(["a", "a", "b", "a", "a", "a", "a"], categories=["a", "b"])
        idx4 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values4 = [1, 1, 1, 1, 1, 1, 1]
        exp_single_cats_value = DataFrame(
            {"cats": cats4, "values": values4}, index=idx4
        )
        return exp_single_cats_value

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc])
    def test_loc_iloc_setitem_list_of_lists(self, orig, exp_multi_row, indexer):
        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()

        key = slice(2, 4)
        if indexer is tm.loc:
            key = slice("j", "k")

        indexer(df)[key, :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        df = orig.copy()
        with pytest.raises(ValueError, match=msg1):
            indexer(df)[key, :] = [["c", 2], ["c", 2]]

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc, tm.at, tm.iat])
    def test_loc_iloc_at_iat_setitem_single_value_in_categories(
        self, orig, exp_single_cats_value, indexer
    ):
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()

        key = (2, 0)
        if indexer in [tm.loc, tm.at]:
            key = (df.index[2], df.columns[0])

        # "b" is among the categories for df["cat"}]
        indexer(df)[key] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        # "c" is not among the categories for df["cat"]
        with pytest.raises(ValueError, match=msg1):
            indexer(df)[key] = "c"

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc])
    def test_loc_iloc_setitem_mask_single_value_in_categories(
        self, orig, exp_single_cats_value, indexer
    ):
        # mask with single True
        df = orig.copy()

        mask = df.index == "j"
        key = 0
        if indexer is tm.loc:
            key = df.columns[key]

        indexer(df)[mask, key] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc])
    def test_iloc_setitem_full_row_non_categorical_rhs(
        self, orig, exp_single_row, indexer
    ):
        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()

        key = 2
        if indexer is tm.loc:
            key = df.index[2]

        # not categorical dtype, but "b" _is_ among the categories for df["cat"]
        indexer(df)[key, :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        # "c" is not among the categories for df["cat"]
        with pytest.raises(ValueError, match=msg1):
            indexer(df)[key, :] = ["c", 2]

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc])
    def test_loc_iloc_setitem_partial_col_categorical_rhs(
        self, orig, exp_parts_cats_col, indexer
    ):
        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()

        key = (slice(2, 4), 0)
        if indexer is tm.loc:
            key = (slice("j", "k"), df.columns[0])

        # same categories as we currently have in df["cats"]
        compat = Categorical(["b", "b"], categories=["a", "b"])
        indexer(df)[key] = compat
        tm.assert_frame_equal(df, exp_parts_cats_col)

        # categories do not match df["cat"]'s, but "b" is among them
        semi_compat = Categorical(list("bb"), categories=list("abc"))
        with pytest.raises(ValueError, match=msg2):
            # different categories but holdable values
            #  -> not sure if this should fail or pass
            indexer(df)[key] = semi_compat

        # categories do not match df["cat"]'s, and "c" is not among them
        incompat = Categorical(list("cc"), categories=list("abc"))
        with pytest.raises(ValueError, match=msg2):
            # different values
            indexer(df)[key] = incompat

    @pytest.mark.parametrize("indexer", [tm.loc, tm.iloc])
    def test_loc_iloc_setitem_non_categorical_rhs(
        self, orig, exp_parts_cats_col, indexer
    ):
        # assign a part of a column with dtype != categorical -> exp_parts_cats_col
        df = orig.copy()

        key = (slice(2, 4), 0)
        if indexer is tm.loc:
            key = (slice("j", "k"), df.columns[0])

        # "b" is among the categories for df["cat"]
        indexer(df)[key] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        # "c" not part of the categories
        with pytest.raises(ValueError, match=msg1):
            indexer(df)[key] = ["c", "c"]

    def test_setitem_mask_categorical(self, exp_multi_row):
        # fancy indexing

        catsf = Categorical(
            ["a", "a", "c", "c", "a", "a", "a"], categories=["a", "b", "c"]
        )
        idxf = Index(["h", "i", "j", "k", "l", "m", "n"])
        valuesf = [1, 1, 3, 3, 1, 1, 1]
        df = DataFrame({"cats": catsf, "values": valuesf}, index=idxf)

        exp_fancy = exp_multi_row.copy()
        return_value = exp_fancy["cats"].cat.set_categories(
            ["a", "b", "c"], inplace=True
        )
        assert return_value is None

        mask = df["cats"] == "c"
        df[mask] = ["b", 2]
        # category c is kept in .categories
        tm.assert_frame_equal(df, exp_fancy)

    def test_loc_setitem_categorical_values_partial_column_slice(
        self, using_array_manager
    ):
        # Assigning a Category to parts of a int/... column uses the values of
        # the Categorical
        df = DataFrame({"a": [1, 1, 1, 1, 1], "b": list("aaaaa")})
        exp = DataFrame({"a": [1, "b", "b", 1, 1], "b": list("aabba")})
        if using_array_manager:
            with pytest.raises(ValueError, match=""):
                df.loc[1:2, "a"] = Categorical(["b", "b"], categories=["a", "b"])
        else:
            df.loc[1:2, "a"] = Categorical(["b", "b"], categories=["a", "b"])
            df.loc[2:3, "b"] = Categorical(["b", "b"], categories=["a", "b"])
            tm.assert_frame_equal(df, exp)

    def test_loc_setitem_single_row_categorical(self, using_array_manager):
        # GH 25495
        df = DataFrame({"Alpha": ["a"], "Numeric": [0]})
        categories = Categorical(df["Alpha"], categories=["a", "b", "c"])
        df.loc[:, "Alpha"] = categories

        result = df["Alpha"]
        expected = Series(categories, index=df.index, name="Alpha")
        if using_array_manager:
            # with ArrayManager the object dtype is preserved
            expected = expected.astype(object)
        tm.assert_series_equal(result, expected)

    def test_loc_indexing_preserves_index_category_dtype(self):
        # GH 15166
        df = DataFrame(
            data=np.arange(2, 22, 2),
            index=pd.MultiIndex(
                levels=[pd.CategoricalIndex(["a", "b"]), range(10)],
                codes=[[0] * 5 + [1] * 5, range(10)],
                names=["Index1", "Index2"],
            ),
        )

        expected = pd.CategoricalIndex(
            ["a", "b"],
            categories=["a", "b"],
            ordered=False,
            name="Index1",
            dtype="category",
        )

        result = df.index.levels[0]
        tm.assert_index_equal(result, expected)

        result = df.loc[["a"]].index.levels[0]
        tm.assert_index_equal(result, expected)
