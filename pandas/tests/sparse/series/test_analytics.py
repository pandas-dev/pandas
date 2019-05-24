import numpy as np
import pytest

from pandas import Series, SparseDtype
from pandas.util import testing as tm


@pytest.mark.parametrize('func', [np.exp, np.sqrt], ids=lambda x: x.__name__)
def test_ufunc(func):
    # GH 23743
    # assert we preserve the incoming dtype on ufunc operation
    s = Series([1, np.nan, 3], dtype=SparseDtype('float64', np.nan))
    result = func(s)
    expected = Series(func([1, np.nan, 3]),
                      dtype=SparseDtype('float64', np.nan))
    tm.assert_series_equal(result, expected)
