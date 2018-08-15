# -*- coding: utf-8 -*-
import pytest

import numpy as np
import pandas as pd

from pandas.compat import long


@pytest.fixture(params=[1, np.array(1, dtype=np.int64)])
def one(request):
    # zero-dim integer array behaves like an integer
    return request.param


zeros = [box_cls([0] * 5, dtype=dtype)
         for box_cls in [pd.Index, np.array]
         for dtype in [np.int64, np.uint64, np.float64]]
zeros.extend([np.array(0, dtype=dtype)
              for dtype in [np.int64, np.uint64, np.float64]])
zeros.extend([0, 0.0, long(0)])


@pytest.fixture(params=zeros)
def zero(request):
    # For testing division by (or of) zero for Index with length 5, this
    # gives several scalar-zeros and length-5 vector-zeros
    return request.param


@pytest.fixture(params=[pd.Float64Index(np.arange(5, dtype='float64')),
                        pd.Int64Index(np.arange(5, dtype='int64')),
                        pd.UInt64Index(np.arange(5, dtype='uint64')),
                        pd.RangeIndex(5)],
                ids=lambda x: type(x).__name__)
def numeric_index(request):
    return request.param


@pytest.fixture
def td_series():
    """
    Return a Series with dtype='timedelta64[ns]', including a NaT.
    """
    return pd.Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')


@pytest.fixture(params=[pd.Timedelta('5m4s').to_pytimedelta(),
                        pd.Timedelta('5m4s'),
                        pd.Timedelta('5m4s').to_timedelta64()],
                ids=lambda x: type(x).__name__)
def scalar_td(request):
    """
    Several variants of Timedelta scalars representing 5 minutes and 4 seconds
    """
    return request.param


# ------------------------------------------------------------------
# DateOffset Fixtures

_common_mismatch = [pd.offsets.YearBegin(2),
                    pd.offsets.MonthBegin(1),
                    pd.offsets.Minute()]


@pytest.fixture(params=[pd.Timedelta(minutes=30).to_pytimedelta(),
                        np.timedelta64(30, 's'),
                        pd.Timedelta(seconds=30)] + _common_mismatch)
def not_hourly(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Hourly frequencies.
    """
    return request.param


@pytest.fixture(params=[np.timedelta64(4, 'h'),
                        pd.Timedelta(hours=23).to_pytimedelta(),
                        pd.Timedelta('23:00:00')] + _common_mismatch)
def not_daily(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Daily frequencies.
    """
    return request.param


@pytest.fixture(params=[np.timedelta64(365, 'D'),
                        pd.Timedelta(365).to_pytimedelta(),
                        pd.Timedelta(days=365)] + _common_mismatch)
def mismatched(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Monthly or Annual frequencies.
    """
    return request.param


@pytest.fixture(params=[pd.offsets.Day(3),
                        pd.Timedelta(days=3).to_pytimedelta(),
                        np.timedelta64(3, 'D'),
                        pd.offsets.Hour(72),
                        pd.Timedelta(minutes=60 * 24 * 3).to_pytimedelta(),
                        np.timedelta64(72, 'h'),
                        pd.Timedelta('72:00:00')])
def three_days(request):
    """
    Several timedelta-like and DateOffset objects that each represent
    a 3-day timedelta
    """
    return request.param


@pytest.fixture(params=[pd.offsets.Hour(2),
                        pd.Timedelta(hours=2),
                        np.timedelta64(2, 'h'),
                        pd.offsets.Minute(120),
                        pd.Timedelta(minutes=120).to_pytimedelta(),
                        np.timedelta64(120, 'm')])
def two_hours(request):
    """
    Several timedelta-like and DateOffset objects that each represent
    a 2-hour timedelta
    """
    return request.param


# ------------------------------------------------------------------

@pytest.fixture(params=[pd.Index, pd.Series, pd.DataFrame],
                ids=lambda x: x.__name__)
def box(request):
    """
    Several array-like containers that should have effectively identical
    behavior with respect to arithmetic operations.
    """
    return request.param


@pytest.fixture(params=[
    pd.Index,
    pd.Series,
    pytest.param(pd.DataFrame,
                 marks=pytest.mark.xfail(reason="Tries to broadcast "
                                                "incorrectly",
                                         strict=True, raises=ValueError))
], ids=lambda x: x.__name__)
def box_df_broadcast_failure(request):
    """
    Fixture equivalent to `box` but with the common failing case where
    the DataFrame operation tries to broadcast incorrectly.
    """
    return request.param
