import numpy as np
import pytest

from pandas import DataFrame, MultiIndex
import pandas.util.testing as tm


class TestMultiIndexPartial:
    def test_getitem_partial_int(self):
        # GH 12416
        # with single item
        l1 = [10, 20]
        l2 = ["a", "b"]
        df = DataFrame(index=range(2), columns=MultiIndex.from_product([l1, l2]))
        expected = DataFrame(index=range(2), columns=l2)
        result = df[20]
        tm.assert_frame_equal(result, expected)

        # with list
        expected = DataFrame(
            index=range(2), columns=MultiIndex.from_product([l1[1:], l2])
        )
        result = df[[20]]
        tm.assert_frame_equal(result, expected)

        # missing item:
        with pytest.raises(KeyError, match="1"):
            df[1]
        with pytest.raises(KeyError, match=r"'\[1\] not in index'"):
            df[[1]]

    def test_series_slice_partial(self):
        pass

    def test_xs_partial(
        self,
        multiindex_dataframe_random_data,
        multiindex_year_month_day_dataframe_random_data,
    ):
        frame = multiindex_dataframe_random_data
        ymd = multiindex_year_month_day_dataframe_random_data
        result = frame.xs("foo")
        result2 = frame.loc["foo"]
        expected = frame.T["foo"].T
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result, result2)

        result = ymd.xs((2000, 4))
        expected = ymd.loc[2000, 4]
        tm.assert_frame_equal(result, expected)

        # ex from #1796
        index = MultiIndex(
            levels=[["foo", "bar"], ["one", "two"], [-1, 1]],
            codes=[
                [0, 0, 0, 0, 1, 1, 1, 1],
                [0, 0, 1, 1, 0, 0, 1, 1],
                [0, 1, 0, 1, 0, 1, 0, 1],
            ],
        )
        df = DataFrame(np.random.randn(8, 4), index=index, columns=list("abcd"))

        result = df.xs(["foo", "one"])
        expected = df.loc["foo", "one"]
        tm.assert_frame_equal(result, expected)

    def test_getitem_partial(self, multiindex_year_month_day_dataframe_random_data):
        ymd = multiindex_year_month_day_dataframe_random_data
        ymd = ymd.T
        result = ymd[2000, 2]

        expected = ymd.reindex(columns=ymd.columns[ymd.columns.codes[1] == 1])
        expected.columns = expected.columns.droplevel(0).droplevel(0)
        tm.assert_frame_equal(result, expected)

    def test_fancy_slice_partial(
        self,
        multiindex_dataframe_random_data,
        multiindex_year_month_day_dataframe_random_data,
    ):
        frame = multiindex_dataframe_random_data
        result = frame.loc["bar":"baz"]
        expected = frame[3:7]
        tm.assert_frame_equal(result, expected)

        ymd = multiindex_year_month_day_dataframe_random_data
        result = ymd.loc[(2000, 2):(2000, 4)]
        lev = ymd.index.codes[1]
        expected = ymd[(lev >= 1) & (lev <= 3)]
        tm.assert_frame_equal(result, expected)

    def test_getitem_partial_column_select(self):
        idx = MultiIndex(
            codes=[[0, 0, 0], [0, 1, 1], [1, 0, 1]],
            levels=[["a", "b"], ["x", "y"], ["p", "q"]],
        )
        df = DataFrame(np.random.rand(3, 2), index=idx)

        result = df.loc[("a", "y"), :]
        expected = df.loc[("a", "y")]
        tm.assert_frame_equal(result, expected)

        result = df.loc[("a", "y"), [1, 0]]
        expected = df.loc[("a", "y")][[1, 0]]
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError, match=r"\('a', 'foo'\)"):
            df.loc[("a", "foo"), :]

    def test_partial_set(self, multiindex_year_month_day_dataframe_random_data):
        # GH #397
        ymd = multiindex_year_month_day_dataframe_random_data
        df = ymd.copy()
        exp = ymd.copy()
        df.loc[2000, 4] = 0
        exp.loc[2000, 4].values[:] = 0
        tm.assert_frame_equal(df, exp)

        df["A"].loc[2000, 4] = 1
        exp["A"].loc[2000, 4].values[:] = 1
        tm.assert_frame_equal(df, exp)

        df.loc[2000] = 5
        exp.loc[2000].values[:] = 5
        tm.assert_frame_equal(df, exp)

        # this works...for now
        df["A"].iloc[14] = 5
        assert df["A"][14] == 5

    # ---------------------------------------------------------------------
    # AMBIGUOUS CASES!

    def test_partial_loc_missing(self, multiindex_year_month_day_dataframe_random_data):
        pytest.skip("skipping for now")

        ymd = multiindex_year_month_day_dataframe_random_data
        result = ymd.loc[2000, 0]
        expected = ymd.loc[2000]["A"]
        tm.assert_series_equal(result, expected)

        # need to put in some work here

        # self.ymd.loc[2000, 0] = 0
        # assert (self.ymd.loc[2000]['A'] == 0).all()

        # Pretty sure the second (and maybe even the first) is already wrong.
        with pytest.raises(Exception):
            ymd.loc[(2000, 6)]
        with pytest.raises(Exception):
            ymd.loc[(2000, 6), 0]

    # ---------------------------------------------------------------------

    def test_setitem_multiple_partial(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        expected = frame.copy()
        result = frame.copy()
        result.loc[["foo", "bar"]] = 0
        expected.loc["foo"] = 0
        expected.loc["bar"] = 0
        tm.assert_frame_equal(result, expected)

        expected = frame.copy()
        result = frame.copy()
        result.loc["foo":"bar"] = 0
        expected.loc["foo"] = 0
        expected.loc["bar"] = 0
        tm.assert_frame_equal(result, expected)

        expected = frame["A"].copy()
        result = frame["A"].copy()
        result.loc[["foo", "bar"]] = 0
        expected.loc["foo"] = 0
        expected.loc["bar"] = 0
        tm.assert_series_equal(result, expected)

        expected = frame["A"].copy()
        result = frame["A"].copy()
        result.loc["foo":"bar"] = 0
        expected.loc["foo"] = 0
        expected.loc["bar"] = 0
        tm.assert_series_equal(result, expected)


def test_loc_getitem_partial_both_axis():
    # gh-12660
    iterables = [["a", "b"], [2, 1]]
    columns = MultiIndex.from_product(iterables, names=["col1", "col2"])
    rows = MultiIndex.from_product(iterables, names=["row1", "row2"])
    df = DataFrame(np.random.randn(4, 4), index=rows, columns=columns)
    expected = df.iloc[:2, 2:].droplevel("row1").droplevel("col1", axis=1)
    result = df.loc["a", "b"]
    tm.assert_frame_equal(result, expected)
