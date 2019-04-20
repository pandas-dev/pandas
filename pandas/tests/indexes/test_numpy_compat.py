import numpy as np
import pytest

from pandas import Float64Index, Index, Int64Index, UInt64Index
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
from pandas.util import testing as tm


@pytest.mark.parametrize(
    'func', [np.exp, np.exp2, np.expm1, np.log, np.log2, np.log10,
             np.log1p, np.sqrt, np.sin, np.cos, np.tan, np.arcsin,
             np.arccos, np.arctan, np.sinh, np.cosh, np.tanh,
             np.arcsinh, np.arccosh, np.arctanh, np.deg2rad,
             np.rad2deg],
    ids=lambda x: x.__name__)
def test_numpy_ufuncs_basic(indices, func):
    # test ufuncs of numpy, see:
    # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    idx = indices
    if isinstance(idx, DatetimeIndexOpsMixin):
        # raise TypeError or ValueError (PeriodIndex)
        # PeriodIndex behavior should be changed in future version
        with pytest.raises(Exception):
            with np.errstate(all='ignore'):
                func(idx)
    elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
        # coerces to float (e.g. np.sin)
        with np.errstate(all='ignore'):
            result = func(idx)
            exp = Index(func(idx.values), name=idx.name)

        tm.assert_index_equal(result, exp)
        assert isinstance(result, Float64Index)
    else:
        # raise AttributeError or TypeError
        if len(idx) == 0:
            pass
        else:
            with pytest.raises(Exception):
                with np.errstate(all='ignore'):
                    func(idx)


@pytest.mark.parametrize(
    'func', [np.isfinite, np.isinf, np.isnan, np.signbit],
    ids=lambda x: x.__name__)
def test_numpy_ufuncs_other(indices, func):
    # test ufuncs of numpy, see:
    # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    idx = indices
    if isinstance(idx, DatetimeIndexOpsMixin):
        # raise TypeError or ValueError (PeriodIndex)
        with pytest.raises(Exception):
            func(idx)
    elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
        # Results in bool array
        result = func(idx)
        assert isinstance(result, np.ndarray)
        assert not isinstance(result, Index)
    else:
        if len(idx) == 0:
            pass
        else:
            with pytest.raises(Exception):
                func(idx)
