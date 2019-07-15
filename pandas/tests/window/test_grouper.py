import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas.util.testing as tm


class TestGrouperGrouping:
    def setup_method(self, method):
        self.series = Series(np.arange(10))
        self.frame = DataFrame({"A": [1] * 20 + [2] * 12 + [3] * 8, "B": np.arange(40)})

    def test_mutated(self):

        msg = r"group\(\) got an unexpected keyword argument 'foo'"
        with pytest.raises(TypeError, match=msg):
            self.frame.groupby("A", foo=1)

        g = self.frame.groupby("A")
        assert not g.mutated
        g = self.frame.groupby("A", mutated=True)
        assert g.mutated

    def test_getitem(self):
        g = self.frame.groupby("A")
        g_mutated = self.frame.groupby("A", mutated=True)

        expected = g_mutated.B.apply(lambda x: x.rolling(2).mean())

        result = g.rolling(2).mean().B
        tm.assert_series_equal(result, expected)

        result = g.rolling(2).B.mean()
        tm.assert_series_equal(result, expected)

        result = g.B.rolling(2).mean()
        tm.assert_series_equal(result, expected)

        result = self.frame.B.groupby(self.frame.A).rolling(2).mean()
        tm.assert_series_equal(result, expected)

    def test_getitem_multiple(self):

        # GH 13174
        g = self.frame.groupby("A")
        r = g.rolling(2)
        g_mutated = self.frame.groupby("A", mutated=True)
        expected = g_mutated.B.apply(lambda x: x.rolling(2).count())

        result = r.B.count()
        tm.assert_series_equal(result, expected)

        result = r.B.count()
        tm.assert_series_equal(result, expected)

    def test_rolling(self):
        g = self.frame.groupby("A")
        r = g.rolling(window=4)

        for f in ["sum", "mean", "min", "max", "count", "kurt", "skew"]:

            result = getattr(r, f)()
            expected = g.apply(lambda x: getattr(x.rolling(4), f)())
            tm.assert_frame_equal(result, expected)

        for f in ["std", "var"]:
            result = getattr(r, f)(ddof=1)
            expected = g.apply(lambda x: getattr(x.rolling(4), f)(ddof=1))
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = g.apply(lambda x: x.rolling(4).quantile(0.5))
        tm.assert_frame_equal(result, expected)

    def test_rolling_corr_cov(self):
        g = self.frame.groupby("A")
        r = g.rolling(window=4)

        for f in ["corr", "cov"]:
            result = getattr(r, f)(self.frame)

            def func(x):
                return getattr(x.rolling(4), f)(self.frame)

            expected = g.apply(func)
            tm.assert_frame_equal(result, expected)

            result = getattr(r.B, f)(pairwise=True)

            def func(x):
                return getattr(x.B.rolling(4), f)(pairwise=True)

            expected = g.apply(func)
            tm.assert_series_equal(result, expected)

    def test_rolling_apply(self, raw):
        g = self.frame.groupby("A")
        r = g.rolling(window=4)

        # reduction
        result = r.apply(lambda x: x.sum(), raw=raw)
        expected = g.apply(lambda x: x.rolling(4).apply(lambda y: y.sum(), raw=raw))
        tm.assert_frame_equal(result, expected)

    def test_rolling_apply_mutability(self):
        # GH 14013
        df = pd.DataFrame({"A": ["foo"] * 3 + ["bar"] * 3, "B": [1] * 6})
        g = df.groupby("A")

        mi = pd.MultiIndex.from_tuples(
            [("bar", 3), ("bar", 4), ("bar", 5), ("foo", 0), ("foo", 1), ("foo", 2)]
        )

        mi.names = ["A", None]
        # Grouped column should not be a part of the output
        expected = pd.DataFrame([np.nan, 2.0, 2.0] * 2, columns=["B"], index=mi)

        result = g.rolling(window=2).sum()
        tm.assert_frame_equal(result, expected)

        # Call an arbitrary function on the groupby
        g.sum()

        # Make sure nothing has been mutated
        result = g.rolling(window=2).sum()
        tm.assert_frame_equal(result, expected)

    def test_expanding(self):
        g = self.frame.groupby("A")
        r = g.expanding()

        for f in ["sum", "mean", "min", "max", "count", "kurt", "skew"]:

            result = getattr(r, f)()
            expected = g.apply(lambda x: getattr(x.expanding(), f)())
            tm.assert_frame_equal(result, expected)

        for f in ["std", "var"]:
            result = getattr(r, f)(ddof=0)
            expected = g.apply(lambda x: getattr(x.expanding(), f)(ddof=0))
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = g.apply(lambda x: x.expanding().quantile(0.5))
        tm.assert_frame_equal(result, expected)

    def test_expanding_corr_cov(self):
        g = self.frame.groupby("A")
        r = g.expanding()

        for f in ["corr", "cov"]:
            result = getattr(r, f)(self.frame)

            def func(x):
                return getattr(x.expanding(), f)(self.frame)

            expected = g.apply(func)
            tm.assert_frame_equal(result, expected)

            result = getattr(r.B, f)(pairwise=True)

            def func(x):
                return getattr(x.B.expanding(), f)(pairwise=True)

            expected = g.apply(func)
            tm.assert_series_equal(result, expected)

    def test_expanding_apply(self, raw):
        g = self.frame.groupby("A")
        r = g.expanding()

        # reduction
        result = r.apply(lambda x: x.sum(), raw=raw)
        expected = g.apply(lambda x: x.expanding().apply(lambda y: y.sum(), raw=raw))
        tm.assert_frame_equal(result, expected)
