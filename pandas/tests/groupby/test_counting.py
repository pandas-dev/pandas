from itertools import product

import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Period, Series, Timedelta, Timestamp
import pandas.util.testing as tm


class TestCounting:
    def test_cumcount(self):
        df = DataFrame([["a"], ["a"], ["a"], ["b"], ["a"]], columns=["A"])
        g = df.groupby("A")
        sg = g.A

        expected = Series([0, 1, 2, 0, 3])

        tm.assert_series_equal(expected, g.cumcount())
        tm.assert_series_equal(expected, sg.cumcount())

    def test_cumcount_empty(self):
        ge = DataFrame().groupby(level=0)
        se = Series(dtype=object).groupby(level=0)

        # edge case, as this is usually considered float
        e = Series(dtype="int64")

        tm.assert_series_equal(e, ge.cumcount())
        tm.assert_series_equal(e, se.cumcount())

    def test_cumcount_dupe_index(self):
        df = DataFrame(
            [["a"], ["a"], ["a"], ["b"], ["a"]], columns=["A"], index=[0] * 5
        )
        g = df.groupby("A")
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        tm.assert_series_equal(expected, g.cumcount())
        tm.assert_series_equal(expected, sg.cumcount())

    def test_cumcount_mi(self):
        mi = MultiIndex.from_tuples([[0, 1], [1, 2], [2, 2], [2, 2], [1, 0]])
        df = DataFrame([["a"], ["a"], ["a"], ["b"], ["a"]], columns=["A"], index=mi)
        g = df.groupby("A")
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=mi)

        tm.assert_series_equal(expected, g.cumcount())
        tm.assert_series_equal(expected, sg.cumcount())

    def test_cumcount_groupby_not_col(self):
        df = DataFrame(
            [["a"], ["a"], ["a"], ["b"], ["a"]], columns=["A"], index=[0] * 5
        )
        g = df.groupby([0, 0, 0, 1, 0])
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        tm.assert_series_equal(expected, g.cumcount())
        tm.assert_series_equal(expected, sg.cumcount())

    def test_ngroup(self):
        df = DataFrame({"A": list("aaaba")})
        g = df.groupby("A")
        sg = g.A

        expected = Series([0, 0, 0, 1, 0])

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_distinct(self):
        df = DataFrame({"A": list("abcde")})
        g = df.groupby("A")
        sg = g.A

        expected = Series(range(5), dtype="int64")

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_one_group(self):
        df = DataFrame({"A": [0] * 5})
        g = df.groupby("A")
        sg = g.A

        expected = Series([0] * 5)

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_empty(self):
        ge = DataFrame().groupby(level=0)
        se = Series(dtype=object).groupby(level=0)

        # edge case, as this is usually considered float
        e = Series(dtype="int64")

        tm.assert_series_equal(e, ge.ngroup())
        tm.assert_series_equal(e, se.ngroup())

    def test_ngroup_series_matches_frame(self):
        df = DataFrame({"A": list("aaaba")})
        s = Series(list("aaaba"))

        tm.assert_series_equal(df.groupby(s).ngroup(), s.groupby(s).ngroup())

    def test_ngroup_dupe_index(self):
        df = DataFrame({"A": list("aaaba")}, index=[0] * 5)
        g = df.groupby("A")
        sg = g.A

        expected = Series([0, 0, 0, 1, 0], index=[0] * 5)

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_mi(self):
        mi = MultiIndex.from_tuples([[0, 1], [1, 2], [2, 2], [2, 2], [1, 0]])
        df = DataFrame({"A": list("aaaba")}, index=mi)
        g = df.groupby("A")
        sg = g.A
        expected = Series([0, 0, 0, 1, 0], index=mi)

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_groupby_not_col(self):
        df = DataFrame({"A": list("aaaba")}, index=[0] * 5)
        g = df.groupby([0, 0, 0, 1, 0])
        sg = g.A

        expected = Series([0, 0, 0, 1, 0], index=[0] * 5)

        tm.assert_series_equal(expected, g.ngroup())
        tm.assert_series_equal(expected, sg.ngroup())

    def test_ngroup_descending(self):
        df = DataFrame(["a", "a", "b", "a", "b"], columns=["A"])
        g = df.groupby(["A"])

        ascending = Series([0, 0, 1, 0, 1])
        descending = Series([1, 1, 0, 1, 0])

        tm.assert_series_equal(descending, (g.ngroups - 1) - ascending)
        tm.assert_series_equal(ascending, g.ngroup(ascending=True))
        tm.assert_series_equal(descending, g.ngroup(ascending=False))

    def test_ngroup_matches_cumcount(self):
        # verify one manually-worked out case works
        df = DataFrame(
            [["a", "x"], ["a", "y"], ["b", "x"], ["a", "x"], ["b", "y"]],
            columns=["A", "X"],
        )
        g = df.groupby(["A", "X"])
        g_ngroup = g.ngroup()
        g_cumcount = g.cumcount()
        expected_ngroup = Series([0, 1, 2, 0, 3])
        expected_cumcount = Series([0, 0, 0, 1, 0])

        tm.assert_series_equal(g_ngroup, expected_ngroup)
        tm.assert_series_equal(g_cumcount, expected_cumcount)

    def test_ngroup_cumcount_pair(self):
        # brute force comparison for all small series
        for p in product(range(3), repeat=4):
            df = DataFrame({"a": p})
            g = df.groupby(["a"])

            order = sorted(set(p))
            ngroupd = [order.index(val) for val in p]
            cumcounted = [p[:i].count(val) for i, val in enumerate(p)]

            tm.assert_series_equal(g.ngroup(), Series(ngroupd))
            tm.assert_series_equal(g.cumcount(), Series(cumcounted))

    def test_ngroup_respects_groupby_order(self):
        np.random.seed(0)
        df = DataFrame({"a": np.random.choice(list("abcdef"), 100)})
        for sort_flag in (False, True):
            g = df.groupby(["a"], sort=sort_flag)
            df["group_id"] = -1
            df["group_index"] = -1

            for i, (_, group) in enumerate(g):
                df.loc[group.index, "group_id"] = i
                for j, ind in enumerate(group.index):
                    df.loc[ind, "group_index"] = j

            tm.assert_series_equal(Series(df["group_id"].values), g.ngroup())
            tm.assert_series_equal(Series(df["group_index"].values), g.cumcount())

    @pytest.mark.parametrize(
        "datetimelike",
        [
            [
                Timestamp("2016-05-{i:02d} 20:09:25+00:00".format(i=i))
                for i in range(1, 4)
            ],
            [Timestamp("2016-05-{i:02d} 20:09:25".format(i=i)) for i in range(1, 4)],
            [Timedelta(x, unit="h") for x in range(1, 4)],
            [Period(freq="2W", year=2017, month=x) for x in range(1, 4)],
        ],
    )
    def test_count_with_datetimelike(self, datetimelike):
        # test for #13393, where DataframeGroupBy.count() fails
        # when counting a datetimelike column.

        df = DataFrame({"x": ["a", "a", "b"], "y": datetimelike})
        res = df.groupby("x").count()
        expected = DataFrame({"y": [2, 1]}, index=["a", "b"])
        expected.index.name = "x"
        tm.assert_frame_equal(expected, res)

    def test_count_with_only_nans_in_first_group(self):
        # GH21956
        df = DataFrame({"A": [np.nan, np.nan], "B": ["a", "b"], "C": [1, 2]})
        result = df.groupby(["A", "B"]).C.count()
        mi = MultiIndex(levels=[[], ["a", "b"]], codes=[[], []], names=["A", "B"])
        expected = Series([], index=mi, dtype=np.int64, name="C")
        tm.assert_series_equal(result, expected, check_index_type=False)
