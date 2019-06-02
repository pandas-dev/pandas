import numpy as np
import pytest

from pandas import Series, SparseDtype
from pandas.util import testing as tm


@pytest.mark.parametrize(
    'data, dtype',
    [([1, np.nan, 3], SparseDtype('float64', np.nan)),
     ([1, 2, 3], SparseDtype('int'))])
@pytest.mark.parametrize('func', [np.exp, np.sqrt], ids=str)
def test_ufunc(data, dtype, func):
    # GH 23743
    # assert we preserve the incoming dtype on ufunc operation
    s = Series(data, dtype=dtype)
    result = func(s)
    expected = Series(func(data),
                      dtype=dtype)
    tm.assert_series_equal(result, expected)
