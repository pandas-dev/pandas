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
def test_numpy_ufuncs_basic(indices, func):
    # test ufuncs of numpy, see:
    # https://docs.scipy.org/doc/numpy/reference/ufuncs.html

    idx = indices
    if isinstance(idx, DatetimeIndexOpsMixin):
        # raise TypeError or ValueError (PeriodIndex)
        with pytest.raises(Exception):
            with np.errstate(all="ignore"):
                func(idx)
    elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
        # coerces to float (e.g. np.sin)
        with np.errstate(all="ignore"):
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
                with np.errstate(all="ignore"):
                    func(idx)


@pytest.mark.parametrize(
    "func", [np.isfinite, np.isinf, np.isnan, np.signbit], ids=lambda x: x.__name__
)
def test_numpy_ufuncs_other(indices, func):
    # test ufuncs of numpy, see:
    # https://docs.scipy.org/doc/numpy/reference/ufuncs.html

    idx = indices
    if isinstance(idx, (DatetimeIndex, TimedeltaIndex)):
        if isinstance(idx, DatetimeIndex) and idx.tz is not None:
            if func in [np.isfinite, np.isnan, np.isinf]:
                pytest.xfail(reason="__array_ufunc__ is not defined")

        if not _np_version_under1p18 and func in [np.isfinite, np.isinf, np.isnan]:
            # numpy 1.18(dev) changed isinf and isnan to not raise on dt64/tfd64
            result = func(idx)
            assert isinstance(result, np.ndarray)

        elif not _np_version_under1p17 and func in [np.isfinite]:
            # ok under numpy >= 1.17
            # Results in bool array
            result = func(idx)
            assert isinstance(result, np.ndarray)
        else:
            # raise TypeError or ValueError (PeriodIndex)
            with pytest.raises(Exception):
                func(idx)

    elif isinstance(idx, PeriodIndex):
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


def test_elementwise_comparison_warning():
    # https://github.com/pandas-dev/pandas/issues/22698#issuecomment-458968300
    # np.array([1, 2]) == 'a' returns False, and produces a
    # FutureWarning that it'll be [False, False] in the future.
    # We just want to ensure that comes through.
    # When NumPy dev actually enforces this change, we'll need to skip
    # this test.
    idx = Index([1, 2])
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = idx == "a"

    expected = np.array([False, False])
    tm.assert_numpy_array_equal(result, expected)
