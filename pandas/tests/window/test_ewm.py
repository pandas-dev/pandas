import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall

from pandas import DataFrame, Series
import pandas.core.window as rwindow
from pandas.tests.window.common import Base


class TestEWM(Base):
    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
        df
        df.ewm(com=0.5).mean()

    @pytest.mark.parametrize("which", ["series", "frame"])
    def test_constructor(self, which):
        o = getattr(self, which)
        c = o.ewm

        # valid
        c(com=0.5)
        c(span=1.5)
        c(alpha=0.5)
        c(halflife=0.75)
        c(com=0.5, span=None)
        c(alpha=0.5, com=None)
        c(halflife=0.75, alpha=None)

        # not valid: mutually exclusive
        with pytest.raises(ValueError):
            c(com=0.5, alpha=0.5)
        with pytest.raises(ValueError):
            c(span=1.5, halflife=0.75)
        with pytest.raises(ValueError):
            c(alpha=0.5, span=1.5)

        # not valid: com < 0
        with pytest.raises(ValueError):
            c(com=-0.5)

        # not valid: span < 1
        with pytest.raises(ValueError):
            c(span=0.5)

        # not valid: halflife <= 0
        with pytest.raises(ValueError):
            c(halflife=0)

        # not valid: alpha <= 0 or alpha > 1
        for alpha in (-0.5, 1.5):
            with pytest.raises(ValueError):
                c(alpha=alpha)

    @pytest.mark.parametrize("method", ["std", "mean", "var"])
    def test_numpy_compat(self, method):
        # see gh-12811
        e = rwindow.EWM(Series([2, 4, 6]), alpha=0.5)

        msg = "numpy operations are not valid with window objects"

        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(e, method)(1, 2, 3)
        with pytest.raises(UnsupportedFunctionCall, match=msg):
            getattr(e, method)(dtype=np.float64)
