import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall

import pandas as pd
from pandas import DataFrame, Series
from pandas.core.window import Expanding
from pandas.tests.window.common import Base
import pandas.util.testing as tm


class TestExpanding(Base):
    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
        df
        df.expanding(2).sum()

    @pytest.mark.parametrize("which", ["series", "frame"])
    def test_constructor(self, which):
        # GH 12669

        o = getattr(self, which)
        c = o.expanding

        # valid
        c(min_periods=1)
        c(min_periods=1, center=True)
        c(min_periods=1, center=False)

        # not valid
        for w in [2.0, "foo", np.array([2])]:
            with pytest.raises(ValueError):
                c(min_periods=w)
            with pytest.raises(ValueError):
                c(min_periods=1, center=w)

    @pytest.mark.parametrize("method", ["std", "mean", "sum", "max", "min", "var"])
    def test_numpy_compat(self, method):
        # see gh-12811
        e = Expanding(Series([2, 4, 6]), window=2)

        msg = "numpy operations are not valid with window objects"

        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(e, method)(1, 2, 3)
        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(e, method)(dtype=np.float64)

    @pytest.mark.parametrize(
        "expander",
        [
            1,
            pytest.param(
                "ls",
                marks=pytest.mark.xfail(
                    reason="GH#16425 expanding with offset not supported"
                ),
            ),
        ],
    )
    def test_empty_df_expanding(self, expander):
        # GH 15819 Verifies that datetime and integer expanding windows can be
        # applied to empty DataFrames

        expected = DataFrame()
        result = DataFrame().expanding(expander).sum()
        tm.assert_frame_equal(result, expected)

        # Verifies that datetime and integer expanding windows can be applied
        # to empty DataFrames with datetime index
        expected = DataFrame(index=pd.DatetimeIndex([]))
        result = DataFrame(index=pd.DatetimeIndex([])).expanding(expander).sum()
        tm.assert_frame_equal(result, expected)

    def test_missing_minp_zero(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        # minp=0
        x = pd.Series([np.nan])
        result = x.expanding(min_periods=0).sum()
        expected = pd.Series([0.0])
        tm.assert_series_equal(result, expected)

        # minp=1
        result = x.expanding(min_periods=1).sum()
        expected = pd.Series([np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("klass", [pd.Series, pd.DataFrame])
    def test_iter_raises(self, klass):
        # https://github.com/pandas-dev/pandas/issues/11704
        # Iteration over a Window
        obj = klass([1, 2, 3, 4])
        with pytest.raises(NotImplementedError):
            iter(obj.expanding(2))

    def test_expanding_axis(self, axis_frame):
        # see gh-23372.
        df = DataFrame(np.ones((10, 20)))
        axis = df._get_axis_number(axis_frame)

        if axis == 0:
            expected = DataFrame(
                {i: [np.nan] * 2 + [float(j) for j in range(3, 11)] for i in range(20)}
            )
        else:
            # axis == 1
            expected = DataFrame([[np.nan] * 2 + [float(i) for i in range(3, 21)]] * 10)

        result = df.expanding(3, axis=axis_frame).sum()
        tm.assert_frame_equal(result, expected)
