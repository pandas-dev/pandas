import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall
import pandas.util._test_decorators as td

import pandas as pd
from pandas import Series
import pandas.core.window as rwindow
from pandas.tests.window.common import Base


@pytest.mark.filterwarnings("ignore:can't resolve package:ImportWarning")
class TestWindow(Base):
    def setup_method(self, method):
        self._create_data()

    @td.skip_if_no_scipy
    @pytest.mark.parametrize("which", ["series", "frame"])
    def test_constructor(self, which):
        # GH 12669

        o = getattr(self, which)
        c = o.rolling

        # valid
        c(win_type="boxcar", window=2, min_periods=1)
        c(win_type="boxcar", window=2, min_periods=1, center=True)
        c(win_type="boxcar", window=2, min_periods=1, center=False)

        # not valid
        for w in [2.0, "foo", np.array([2])]:
            with pytest.raises(ValueError):
                c(win_type="boxcar", window=2, min_periods=w)
            with pytest.raises(ValueError):
                c(win_type="boxcar", window=2, min_periods=1, center=w)

        for wt in ["foobar", 1]:
            with pytest.raises(ValueError):
                c(win_type=wt, window=2)

    @td.skip_if_no_scipy
    @pytest.mark.parametrize("which", ["series", "frame"])
    def test_constructor_with_win_type(self, which, win_types):
        # GH 12669
        o = getattr(self, which)
        c = o.rolling
        c(win_type=win_types, window=2)

    @pytest.mark.parametrize("method", ["sum", "mean"])
    def test_numpy_compat(self, method):
        # see gh-12811
        w = rwindow.Window(Series([2, 4, 6]), window=[0, 2])

        msg = "numpy operations are not valid with window objects"

        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(w, method)(1, 2, 3)
        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(w, method)(dtype=np.float64)

    @td.skip_if_no_scipy
    @pytest.mark.parametrize("arg", ["median", "var", "std", "kurt", "skew"])
    def test_agg_function_support(self, arg):
        df = pd.DataFrame({"A": np.arange(5)})
        roll = df.rolling(2, win_type="triang")

        msg = "'{arg}' is not a valid function for " "'Window' object".format(arg=arg)
        with pytest.raises(AttributeError, match=msg):
            roll.agg(arg)

        with pytest.raises(AttributeError, match=msg):
            roll.agg([arg])

        with pytest.raises(AttributeError, match=msg):
            roll.agg({"A": arg})
