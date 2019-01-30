import numpy as np
import pytest

from pandas.core.arrays.numpy_ import PandasArray


@pytest.fixture
def allow_in_pandas(monkeypatch):
    """
    A monkeypatch to tell pandas to let us in.

    By default, passing a PandasArray to an index / series / frame
    constructor will unbox that PandasArray to an ndarray, and treat
    it as a non-EA column. We don't want people using EAs without
    reason.

    The mechanism for this is a check against ABCPandasArray
    in each constructor.

    But, for testing, we need to allow them in pandas. So we patch
    the _typ of PandasArray, so that we evade the ABCPandasArray
    check.
    """
    with monkeypatch.context() as m:
        m.setattr(PandasArray, '_typ', 'extension')
        yield


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def na_cmp():
    def cmp(a, b):
        return np.isnan(a) and np.isnan(b)
    return cmp
