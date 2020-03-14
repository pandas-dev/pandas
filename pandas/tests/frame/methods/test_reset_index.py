from datetime import datetime

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    IntervalIndex,
    MultiIndex,
    RangeIndex,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm


class TestResetIndex:
    def test_reset_index_tz(self, tz_aware_fixture):
        # GH 3950
        # reset_index with single level
        tz = tz_aware_fixture
        idx = date_range("1/1/2011", periods=5, freq="D", tz=tz, name="idx")
        df = DataFrame({"a": range(5), "b": ["A", "B", "C", "D", "E"]}, index=idx)

        expected = DataFrame(
            {
                "idx": [
                    datetime(2011, 1, 1),
                    datetime(2011, 1, 2),
                    datetime(2011, 1, 3),
                    datetime(2011, 1, 4),
                    datetime(2011, 1, 5),
                ],
                "a": range(5),
                "b": ["A", "B", "C", "D", "E"],
            },
            columns=["idx", "a", "b"],
        )
        expected["idx"] = expected["idx"].apply(lambda d: Timestamp(d, tz=tz))
        tm.assert_frame_equal(df.reset_index(), expected)

    def test_reset_index_with_intervals(self):
        idx = IntervalIndex.from_breaks(np.arange(11), name="x")
        original = DataFrame({"x": idx, "y": np.arange(10)})[["x", "y"]]

        result = original.set_index("x")
        expected = DataFrame({"y": np.arange(10)}, index=idx)
        tm.assert_frame_equal(result, expected)

        result2 = result.reset_index()
        tm.assert_frame_equal(result2, original)

    def test_reset_index(self, float_frame):
        stacked = float_frame.stack()[::2]
        stacked = DataFrame({"foo": stacked, "bar": stacked})

        names = ["first", "second"]
        stacked.index.names = names
        deleveled = stacked.reset_index()
        for i, (lev, level_codes) in enumerate(
            zip(stacked.index.levels, stacked.index.codes)
        ):
            values = lev.take(level_codes)
            name = names[i]
            tm.assert_index_equal(values, Index(deleveled[name]))

        stacked.index.names = [None, None]
        deleveled2 = stacked.reset_index()
        tm.assert_series_equal(
            deleveled["first"], deleveled2["level_0"], check_names=False
        )
        tm.assert_series_equal(
            deleveled["second"], deleveled2["level_1"], check_names=False
        )

        # default name assigned
        rdf = float_frame.reset_index()
        exp = Series(float_frame.index.values, name="index")
        tm.assert_series_equal(rdf["index"], exp)

        # default name assigned, corner case
        df = float_frame.copy()
        df["index"] = "foo"
        rdf = df.reset_index()
        exp = Series(float_frame.index.values, name="level_0")
        tm.assert_series_equal(rdf["level_0"], exp)

        # but this is ok
        float_frame.index.name = "index"
        deleveled = float_frame.reset_index()
        tm.assert_series_equal(deleveled["index"], Series(float_frame.index))
        tm.assert_index_equal(deleveled.index, Index(np.arange(len(deleveled))))

        # preserve column names
        float_frame.columns.name = "columns"
        resetted = float_frame.reset_index()
        assert resetted.columns.name == "columns"

        # only remove certain columns
        df = float_frame.reset_index().set_index(["index", "A", "B"])
        rs = df.reset_index(["A", "B"])

        # TODO should reset_index check_names ?
        tm.assert_frame_equal(rs, float_frame, check_names=False)

        rs = df.reset_index(["index", "A", "B"])
        tm.assert_frame_equal(rs, float_frame.reset_index(), check_names=False)

        rs = df.reset_index(["index", "A", "B"])
        tm.assert_frame_equal(rs, float_frame.reset_index(), check_names=False)

        rs = df.reset_index("A")
        xp = float_frame.reset_index().set_index(["index", "B"])
        tm.assert_frame_equal(rs, xp, check_names=False)

        # test resetting in place
        df = float_frame.copy()
        resetted = float_frame.reset_index()
        df.reset_index(inplace=True)
        tm.assert_frame_equal(df, resetted, check_names=False)

        df = float_frame.reset_index().set_index(["index", "A", "B"])
        rs = df.reset_index("A", drop=True)
        xp = float_frame.copy()
        del xp["A"]
        xp = xp.set_index(["B"], append=True)
        tm.assert_frame_equal(rs, xp, check_names=False)

    def test_reset_index_name(self):
        df = DataFrame(
            [[1, 2, 3, 4], [5, 6, 7, 8]],
            columns=["A", "B", "C", "D"],
            index=Index(range(2), name="x"),
        )
        assert df.reset_index().index.name is None
        assert df.reset_index(drop=True).index.name is None
        df.reset_index(inplace=True)
        assert df.index.name is None

    def test_reset_index_level(self):
        df = DataFrame([[1, 2, 3, 4], [5, 6, 7, 8]], columns=["A", "B", "C", "D"])

        for levels in ["A", "B"], [0, 1]:
            # With MultiIndex
            result = df.set_index(["A", "B"]).reset_index(level=levels[0])
            tm.assert_frame_equal(result, df.set_index("B"))

            result = df.set_index(["A", "B"]).reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df.set_index("B"))

            result = df.set_index(["A", "B"]).reset_index(level=levels)
            tm.assert_frame_equal(result, df)

            result = df.set_index(["A", "B"]).reset_index(level=levels, drop=True)
            tm.assert_frame_equal(result, df[["C", "D"]])

            # With single-level Index (GH 16263)
            result = df.set_index("A").reset_index(level=levels[0])
            tm.assert_frame_equal(result, df)

            result = df.set_index("A").reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df)

            result = df.set_index(["A"]).reset_index(level=levels[0], drop=True)
            tm.assert_frame_equal(result, df[["B", "C", "D"]])

        # Missing levels - for both MultiIndex and single-level Index:
        for idx_lev in ["A", "B"], ["A"]:
            with pytest.raises(KeyError, match=r"(L|l)evel \(?E\)?"):
                df.set_index(idx_lev).reset_index(level=["A", "E"])
            with pytest.raises(IndexError, match="Too many levels"):
                df.set_index(idx_lev).reset_index(level=[0, 1, 2])

    def test_reset_index_right_dtype(self):
        time = np.arange(0.0, 10, np.sqrt(2) / 2)
        s1 = Series(
            (9.81 * time ** 2) / 2, index=Index(time, name="time"), name="speed"
        )
        df = DataFrame(s1)

        resetted = s1.reset_index()
        assert resetted["time"].dtype == np.float64

        resetted = df.reset_index()
        assert resetted["time"].dtype == np.float64

    def test_reset_index_multiindex_col(self):
        vals = np.random.randn(3, 3).astype(object)
        idx = ["x", "y", "z"]
        full = np.hstack(([[x] for x in idx], vals))
        df = DataFrame(
            vals,
            Index(idx, name="a"),
            columns=[["b", "b", "c"], ["mean", "median", "mean"]],
        )
        rs = df.reset_index()
        xp = DataFrame(
            full, columns=[["a", "b", "b", "c"], ["", "mean", "median", "mean"]]
        )
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index(col_fill=None)
        xp = DataFrame(
            full, columns=[["a", "b", "b", "c"], ["a", "mean", "median", "mean"]]
        )
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index(col_level=1, col_fill="blah")
        xp = DataFrame(
            full, columns=[["blah", "b", "b", "c"], ["a", "mean", "median", "mean"]]
        )
        tm.assert_frame_equal(rs, xp)

        df = DataFrame(
            vals,
            MultiIndex.from_arrays([[0, 1, 2], ["x", "y", "z"]], names=["d", "a"]),
            columns=[["b", "b", "c"], ["mean", "median", "mean"]],
        )
        rs = df.reset_index("a")
        xp = DataFrame(
            full,
            Index([0, 1, 2], name="d"),
            columns=[["a", "b", "b", "c"], ["", "mean", "median", "mean"]],
        )
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index("a", col_fill=None)
        xp = DataFrame(
            full,
            Index(range(3), name="d"),
            columns=[["a", "b", "b", "c"], ["a", "mean", "median", "mean"]],
        )
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index("a", col_fill="blah", col_level=1)
        xp = DataFrame(
            full,
            Index(range(3), name="d"),
            columns=[["blah", "b", "b", "c"], ["a", "mean", "median", "mean"]],
        )
        tm.assert_frame_equal(rs, xp)

    def test_reset_index_multiindex_nan(self):
        # GH#6322, testing reset_index on MultiIndexes
        # when we have a nan or all nan
        df = DataFrame(
            {"A": ["a", "b", "c"], "B": [0, 1, np.nan], "C": np.random.rand(3)}
        )
        rs = df.set_index(["A", "B"]).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame(
            {"A": [np.nan, "b", "c"], "B": [0, 1, 2], "C": np.random.rand(3)}
        )
        rs = df.set_index(["A", "B"]).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame({"A": ["a", "b", "c"], "B": [0, 1, 2], "C": [np.nan, 1.1, 2.2]})
        rs = df.set_index(["A", "B"]).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame(
            {
                "A": ["a", "b", "c"],
                "B": [np.nan, np.nan, np.nan],
                "C": np.random.rand(3),
            }
        )
        rs = df.set_index(["A", "B"]).reset_index()
        tm.assert_frame_equal(rs, df)

    def test_reset_index_with_datetimeindex_cols(self):
        # GH#5818
        df = DataFrame(
            [[1, 2], [3, 4]],
            columns=date_range("1/1/2013", "1/2/2013"),
            index=["A", "B"],
        )

        result = df.reset_index()
        expected = DataFrame(
            [["A", 1, 2], ["B", 3, 4]],
            columns=["index", datetime(2013, 1, 1), datetime(2013, 1, 2)],
        )
        tm.assert_frame_equal(result, expected)

    def test_reset_index_range(self):
        # GH#12071
        df = DataFrame([[0, 0], [1, 1]], columns=["A", "B"], index=RangeIndex(stop=2))
        result = df.reset_index()
        assert isinstance(result.index, RangeIndex)
        expected = DataFrame(
            [[0, 0, 0], [1, 1, 1]],
            columns=["index", "A", "B"],
            index=RangeIndex(stop=2),
        )
        tm.assert_frame_equal(result, expected)
