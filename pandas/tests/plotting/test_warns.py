import pytest

import pandas as pd
import pandas.util.testing as tm

from pandas.plotting import _converter
from pandas.tseries import converter


plt = pytest.importorskip('matplotlib.pyplot')


def test_warns():
    s = pd.Series(range(12), index=pd.date_range('2017', periods=12))
    fig, ax = plt.subplots()

    _converter._WARN = True
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False) as w:
        ax.plot(s.index, s.values)
        plt.close()

    assert len(w) == 1


def test_registering_no_warning():
    s = pd.Series(range(12), index=pd.date_range('2017', periods=12))
    fig, ax = plt.subplots()

    _converter._WARN = True
    converter.register()
    with tm.assert_produces_warning(None) as w:
        ax.plot(s.index, s.values)

    assert len(w) == 0
