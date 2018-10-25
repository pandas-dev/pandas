"""
Tests for pd.NaT that are logically dependent only on pandas._libs.tslibs;
any tests that depend on Index/Series/etc are in test_vector_compat.py
"""
from datetime import datetime, timedelta

import pytest
import pytz
import numpy as np

from pandas import NaT, Timestamp, Timedelta, Period, isna
from pandas.util import testing as tm
from pandas._libs.tslib import iNaT


@pytest.mark.parametrize('klass', [Timestamp, Timedelta, Period])
def test_identity(klass):
    assert klass(None) is NaT

    result = klass(np.nan)
    assert result is NaT

    result = klass(None)
    assert result is NaT

    result = klass(iNaT)
    assert result is NaT

    result = klass(np.nan)
    assert result is NaT

    result = klass(float('nan'))
    assert result is NaT

    result = klass(NaT)
    assert result is NaT

    result = klass('NaT')
    assert result is NaT

    assert isna(klass('nat'))


@pytest.mark.parametrize('klass', [Timestamp, Timedelta, Period])
def test_equality(klass):

    # nat
    if klass is not Period:
        klass('').value == iNaT
    klass('nat').value == iNaT
    klass('NAT').value == iNaT
    klass(None).value == iNaT
    klass(np.nan).value == iNaT
    assert isna(klass('nat'))


@pytest.mark.parametrize('klass', [Timestamp, Timedelta])
def test_round_nat(klass):
    # GH#14940
    ts = klass('nat')
    for method in ["round", "floor", "ceil"]:
        round_method = getattr(ts, method)
        for freq in ["s", "5s", "min", "5min", "h", "5h"]:
            assert round_method(freq) is ts


def test_NaT_methods():
    # GH#9513
    # GH#17329 for `timestamp`
    raise_methods = ['astimezone', 'combine', 'ctime', 'dst',
                     'fromordinal', 'fromtimestamp', 'isocalendar',
                     'strftime', 'strptime', 'time', 'timestamp',
                     'timetuple', 'timetz', 'toordinal', 'tzname',
                     'utcfromtimestamp', 'utcnow', 'utcoffset',
                     'utctimetuple', 'timestamp']
    nat_methods = ['date', 'now', 'replace', 'to_datetime', 'today',
                   'tz_convert', 'tz_localize']
    nan_methods = ['weekday', 'isoweekday']

    for method in raise_methods:
        if hasattr(NaT, method):
            with pytest.raises(ValueError):
                getattr(NaT, method)()

    for method in nan_methods:
        if hasattr(NaT, method):
            assert np.isnan(getattr(NaT, method)())

    for method in nat_methods:
        if hasattr(NaT, method):
            # see GH#8254
            exp_warning = None
            if method == 'to_datetime':
                exp_warning = FutureWarning
            with tm.assert_produces_warning(
                    exp_warning, check_stacklevel=False):
                assert getattr(NaT, method)() is NaT

    # GH#12300
    assert NaT.isoformat() == 'NaT'


def test_NaT_docstrings():
    # GH#17327
    nat_names = dir(NaT)

    # NaT should have *most* of the Timestamp methods, with matching
    # docstrings.  The attributes that are not expected to be present in NaT
    # are private methods plus `ts_expected` below.
    ts_names = dir(Timestamp)
    ts_missing = [x for x in ts_names if x not in nat_names and
                  not x.startswith('_')]
    ts_missing.sort()
    ts_expected = ['freqstr', 'normalize',
                   'to_julian_date',
                   'to_period', 'tz']
    assert ts_missing == ts_expected

    ts_overlap = [x for x in nat_names if x in ts_names and
                  not x.startswith('_') and
                  callable(getattr(Timestamp, x))]
    for name in ts_overlap:
        tsdoc = getattr(Timestamp, name).__doc__
        natdoc = getattr(NaT, name).__doc__
        assert tsdoc == natdoc

    # NaT should have *most* of the Timedelta methods, with matching
    # docstrings.  The attributes that are not expected to be present in NaT
    # are private methods plus `td_expected` below.
    # For methods that are both Timestamp and Timedelta methods, the
    # Timestamp docstring takes priority.
    td_names = dir(Timedelta)
    td_missing = [x for x in td_names if x not in nat_names and
                  not x.startswith('_')]
    td_missing.sort()
    td_expected = ['components', 'delta', 'is_populated',
                   'to_pytimedelta', 'to_timedelta64', 'view']
    assert td_missing == td_expected

    td_overlap = [x for x in nat_names if x in td_names and
                  x not in ts_names and  # Timestamp __doc__ takes priority
                  not x.startswith('_') and
                  callable(getattr(Timedelta, x))]
    assert td_overlap == ['total_seconds']
    for name in td_overlap:
        tddoc = getattr(Timedelta, name).__doc__
        natdoc = getattr(NaT, name).__doc__
        assert tddoc == natdoc


@pytest.mark.parametrize('klass', [Timestamp, Timedelta])
def test_isoformat(klass):

    result = klass('NaT').isoformat()
    expected = 'NaT'
    assert result == expected


# TODO: split this test up?
def test_nat_arithmetic():
    # GH#6873
    i = 2
    f = 1.5

    for (left, right) in [(NaT, i), (NaT, f), (NaT, np.nan)]:
        assert left / right is NaT
        assert left * right is NaT
        assert right * left is NaT
        with pytest.raises(TypeError):
            right / left

    # Timestamp / datetime
    t = Timestamp('2014-01-01')
    dt = datetime(2014, 1, 1)
    for (left, right) in [(NaT, NaT), (NaT, t), (NaT, dt)]:
        # NaT __add__ or __sub__ Timestamp-like (or inverse) returns NaT
        assert right + left is NaT
        assert left + right is NaT
        assert left - right is NaT
        assert right - left is NaT

    # timedelta-like
    # offsets are tested in test_offsets.py

    delta = timedelta(3600)
    td = Timedelta('5s')

    for (left, right) in [(NaT, delta), (NaT, td)]:
        # NaT + timedelta-like returns NaT
        assert right + left is NaT
        assert left + right is NaT
        assert right - left is NaT
        assert left - right is NaT
        assert np.isnan(left / right)
        assert np.isnan(right / left)

    # GH#11718
    t_utc = Timestamp('2014-01-01', tz='UTC')
    t_tz = Timestamp('2014-01-01', tz='US/Eastern')
    dt_tz = pytz.timezone('Asia/Tokyo').localize(dt)

    for (left, right) in [(NaT, t_utc), (NaT, t_tz),
                          (NaT, dt_tz)]:
        # NaT __add__ or __sub__ Timestamp-like (or inverse) returns NaT
        assert right + left is NaT
        assert left + right is NaT
        assert left - right is NaT
        assert right - left is NaT

    # int addition / subtraction
    for (left, right) in [(NaT, 2), (NaT, 0), (NaT, -3)]:
        assert right + left is NaT
        assert left + right is NaT
        assert left - right is NaT
        assert right - left is NaT


def test_nat_rfloordiv_timedelta():
    # GH#18846
    # See also test_timedelta.TestTimedeltaArithmetic.test_floordiv
    td = Timedelta(hours=3, minutes=4)

    assert td // np.nan is NaT
    assert np.isnan(td // NaT)
    assert np.isnan(td // np.timedelta64('NaT'))


def test_nat_pinned_docstrings():
    # GH#17327
    assert NaT.ctime.__doc__ == datetime.ctime.__doc__
