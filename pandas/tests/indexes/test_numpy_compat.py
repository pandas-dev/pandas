import numpy as np
import pytest

from pandas import (
    DatetimeIndex,
    Float64Index,
    Index,
    Int64Index,
    PeriodIndex,
    TimedeltaIndex,
    UInt64Index,
    _np_version_under1p17,
    _np_version_under1p18,
)
import pandas._testing as tm
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin


@pytest.mark.parametrize(
    "func",
    [
        np.exp,
        np.exp2,
        np.expm1,
        np.log,
        np.log2,
        np.log10,
        np.log1p,
        np.sqrt,
        np.sin,
        np.cos,
        np.tan,
        np.arcsin,
        np.arccos,
        np.arctan,
        np.sinh,
        np.cosh,
        np.tanh,
        np.arcsinh,
        np.arccosh,
        np.arctanh,
        np.deg2rad,
        np.rad2deg,
    ],
    ids=lambda x: x.__name__,
)
def test_numpy_ufuncs_basic(index, func):
    # test ufuncs of numpy, see:
    # https://numpy.org/doc/stable/reference/ufuncs.html

    if isinstance(index, DatetimeIndexOpsMixin):
        # raise TypeError or ValueError (PeriodIndex)
        with pytest.raises(Exception):
            with np.errstate(all="ignore"):
                func(index)
    elif isinstance(index, (Float64Index, Int64Index, UInt64Index)):
        # coerces to float (e.g. np.sin)
        with np.errstate(all="ignore"):
            result = func(index)
            exp = Index(func(index.values), name=index.name)

        tm.assert_index_equal(result, exp)
        assert isinstance(result, Float64Index)
    else:
        # raise AttributeError or TypeError
        if len(index) == 0:
            pass
        else:
            with pytest.raises(Exception):
                with np.errstate(all="ignore"):
                    func(index)


@pytest.mark.parametrize(
    "func", [np.isfinite, np.isinf, np.isnan, np.signbit], ids=lambda x: x.__name__
)
def test_numpy_ufuncs_other(index, func):
    # test ufuncs of numpy, see:
    # https://numpy.org/doc/stable/reference/ufuncs.html

    if isinstance(index, (DatetimeIndex, TimedeltaIndex)):
        if isinstance(index, DatetimeIndex) and index.tz is not None:
            if func in [np.isfinite, np.isnan, np.isinf]:
                pytest.xfail(reason="__array_ufunc__ is not defined")

        if not _np_version_under1p18 and func in [np.isfinite, np.isinf, np.isnan]:
            # numpy 1.18(dev) changed isinf and isnan to not raise on dt64/tfd64
            result = func(index)
            assert isinstance(result, np.ndarray)

        elif not _np_version_under1p17 and func in [np.isfinite]:
            # ok under numpy >= 1.17
            # Results in bool array
            result = func(index)
            assert isinstance(result, np.ndarray)
        else:
            # raise TypeError or ValueError (PeriodIndex)
            with pytest.raises(Exception):
                func(index)

    elif isinstance(index, PeriodIndex):
        # raise TypeError or ValueError (PeriodIndex)
        with pytest.raises(Exception):
            func(index)

    elif isinstance(index, (Float64Index, Int64Index, UInt64Index)):
        # Results in bool array
        result = func(index)
        assert isinstance(result, np.ndarray)
        assert not isinstance(result, Index)
    else:
        if len(index) == 0:
            pass
        else:
            with pytest.raises(Exception):
                func(index)
