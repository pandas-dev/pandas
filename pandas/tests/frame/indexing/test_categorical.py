import numpy as np
import pytest

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import Categorical, DataFrame, Index, Series
import pandas._testing as tm


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

    def test_assigning_ops(self):
        # systematically test the assigning operations:
        # for all slicing ops:
        #  for value in categories and value not in categories:

        #   - assign a single value -> exp_single_cats_value

        #   - assign a complete row (mixed values) -> exp_single_row

        # assign multiple rows (mixed values) (-> array) -> exp_multi_row

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col

        cats = Categorical(["a", "a", "a", "a", "a", "a", "a"], categories=["a", "b"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n"])
        values = [1, 1, 1, 1, 1, 1, 1]
        orig = DataFrame({"cats": cats, "values": values}, index=idx)

        # the expected values
        # changed single row
        cats1 = Categorical(["a", "a", "b", "a", "a", "a", "a"], categories=["a", "b"])
        idx1 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values1 = [1, 1, 2, 1, 1, 1, 1]
        exp_single_row = DataFrame({"cats": cats1, "values": values1}, index=idx1)

        # changed multiple rows
        cats2 = Categorical(["a", "a", "b", "b", "a", "a", "a"], categories=["a", "b"])
        idx2 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values2 = [1, 1, 2, 2, 1, 1, 1]
        exp_multi_row = DataFrame({"cats": cats2, "values": values2}, index=idx2)

        # changed part of the cats column
        cats3 = Categorical(["a", "a", "b", "b", "a", "a", "a"], categories=["a", "b"])
        idx3 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values3 = [1, 1, 1, 1, 1, 1, 1]
        exp_parts_cats_col = DataFrame({"cats": cats3, "values": values3}, index=idx3)

        # changed single value in cats col
        cats4 = Categorical(["a", "a", "b", "a", "a", "a", "a"], categories=["a", "b"])
        idx4 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values4 = [1, 1, 1, 1, 1, 1, 1]
        exp_single_cats_value = DataFrame(
            {"cats": cats4, "values": values4}, index=idx4
        )

        #  iloc
        # ###############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.iloc[2, 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.iloc[df.index == "j", 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        msg1 = (
            "Cannot setitem on a Categorical with a new category, "
            "set the categories first"
        )
        msg2 = "Cannot set a Categorical with another, without identical categories"
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.iloc[2, 0] = "c"

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.iloc[2, :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.iloc[2, :] = ["c", 2]

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.iloc[2:4, :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.iloc[2:4, :] = [["c", 2], ["c", 2]]

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4, 0] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg2):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.iloc[2:4, 0] = Categorical(list("bb"), categories=list("abc"))

        with pytest.raises(ValueError, match=msg2):
            # different values
            df = orig.copy()
            df.iloc[2:4, 0] = Categorical(list("cc"), categories=list("abc"))

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4, 0] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg1):
            df.iloc[2:4, 0] = ["c", "c"]

        #  loc
        # ##############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.loc["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.loc[df.index == "j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j", "cats"] = "c"

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.loc["j", :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j", :] = ["c", 2]

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.loc["j":"k", :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j":"k", :] = [["c", 2], ["c", 2]]

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", "cats"] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg2):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.loc["j":"k", "cats"] = Categorical(
                ["b", "b"], categories=["a", "b", "c"]
            )

        with pytest.raises(ValueError, match=msg2):
            # different values
            df = orig.copy()
            df.loc["j":"k", "cats"] = Categorical(
                ["c", "c"], categories=["a", "b", "c"]
            )

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", "cats"] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg1):
            df.loc["j":"k", "cats"] = ["c", "c"]

        #  loc
        # ##############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.loc["j", df.columns[0]] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.loc[df.index == "j", df.columns[0]] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j", df.columns[0]] = "c"

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.loc["j", :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j", :] = ["c", 2]

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.loc["j":"k", :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.loc["j":"k", :] = [["c", 2], ["c", 2]]

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", df.columns[0]] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg2):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.loc["j":"k", df.columns[0]] = Categorical(
                ["b", "b"], categories=["a", "b", "c"]
            )

        with pytest.raises(ValueError, match=msg2):
            # different values
            df = orig.copy()
            df.loc["j":"k", df.columns[0]] = Categorical(
                ["c", "c"], categories=["a", "b", "c"]
            )

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", df.columns[0]] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError, match=msg1):
            df.loc["j":"k", df.columns[0]] = ["c", "c"]

        # iat
        df = orig.copy()
        df.iat[2, 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.iat[2, 0] = "c"

        # at
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.at["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.at["j", "cats"] = "c"

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

        df[df["cats"] == "c"] = ["b", 2]
        # category c is kept in .categories
        tm.assert_frame_equal(df, exp_fancy)

        # set_value
        df = orig.copy()
        df.at["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        with pytest.raises(ValueError, match=msg1):
            df = orig.copy()
            df.at["j", "cats"] = "c"

        # Assigning a Category to parts of a int/... column uses the values of
        # the Categorical
        df = DataFrame({"a": [1, 1, 1, 1, 1], "b": list("aaaaa")})
        exp = DataFrame({"a": [1, "b", "b", 1, 1], "b": list("aabba")})
        df.loc[1:2, "a"] = Categorical(["b", "b"], categories=["a", "b"])
        df.loc[2:3, "b"] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp)

    def test_loc_setitem_single_row_categorical(self):
        # GH 25495
        df = DataFrame({"Alpha": ["a"], "Numeric": [0]})
        categories = Categorical(df["Alpha"], categories=["a", "b", "c"])
        df.loc[:, "Alpha"] = categories

        result = df["Alpha"]
        expected = Series(categories, index=df.index, name="Alpha")
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
