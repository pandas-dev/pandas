import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
from pandas.tests.plotting.common import TestPlotBase, _check_plot_works


@td.skip_if_no_mpl
class TestCommon(TestPlotBase):
    def test__check_ticks_props(self):
        # GH 34768
        df = DataFrame({"b": [0, 1, 0], "a": [1, 2, 3]})
        ax = _check_plot_works(df.plot, rot=30)
        ax.yaxis.set_tick_params(rotation=30)
        msg = "expected 0.00000 but got "
        with pytest.raises(AssertionError, match=msg):
            self._check_ticks_props(ax, xrot=0)
        with pytest.raises(AssertionError, match=msg):
            self._check_ticks_props(ax, xlabelsize=0)
        with pytest.raises(AssertionError, match=msg):
            self._check_ticks_props(ax, yrot=0)
        with pytest.raises(AssertionError, match=msg):
            self._check_ticks_props(ax, ylabelsize=0)
