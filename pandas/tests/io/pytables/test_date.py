import pytest
import datetime

import numpy as np
import pandas as pd
from pandas import Series, DataFrame, date_range, timedelta_range, bdate_range
from pandas.tests.io.pytables.common import (ensure_clean_store,
                                             _check_roundtrip)
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal, assert_series_equal


def test_calendar_roundtrip_issue():
    # 8591
    # doc example from tseries holiday section
    weekmask_egypt = 'Sun Mon Tue Wed Thu'
    holidays = ['2012-05-01',
                datetime.datetime(2013, 5, 1), np.datetime64('2014-05-01')]
    bday_egypt = pd.offsets.CustomBusinessDay(
        holidays=holidays, weekmask=weekmask_egypt)
    dt = datetime.datetime(2013, 4, 30)
    dts = date_range(dt, periods=5, freq=bday_egypt)

    s = (Series(dts.weekday, dts).map(
        Series('Mon Tue Wed Thu Fri Sat Sun'.split())))

    with ensure_clean_store() as store:
        store.put('fixed', s)
        result = store.select('fixed')
        assert_series_equal(result, s)

        store.append('table', s)
        result = store.select('table')
        assert_series_equal(result, s)


def test_roundtrip_tz_aware_index():
    # GH 17618
    time = pd.Timestamp('2000-01-01 01:00:00', tz='US/Eastern')
    df = pd.DataFrame(data=[0], index=[time])

    with ensure_clean_store() as store:
        store.put('frame', df, format='fixed')
        recons = store['frame']
        tm.assert_frame_equal(recons, df)
        assert recons.index[0].value == 946706400000000000


def test_timeseries_preepoch():
    dr = bdate_range('1/1/1940', '1/1/1960')
    ts = Series(np.random.randn(len(dr)), index=dr)
    try:
        _check_roundtrip(ts, tm.assert_series_equal)
    except OverflowError:
        pytest.skip('known failer on some windows platforms')


def test_can_serialize_dates():
    rng = [x.date() for x in bdate_range('1/1/2000', '1/30/2000')]
    frame = DataFrame(np.random.randn(len(rng), 4), index=rng)
    _check_roundtrip(frame, tm.assert_frame_equal)


def test_store_datetime_fractional_secs():
    with ensure_clean_store() as store:
        dt = datetime.datetime(2012, 1, 2, 3, 4, 5, 123456)
        series = Series([0], [dt])
        store['a'] = series
        assert store['a'].index[0] == dt


def test_tseries_indices_series():
    with ensure_clean_store() as store:
        idx = tm.makeDateIndex(10)
        ser = Series(np.random.randn(len(idx)), idx)
        store['a'] = ser
        result = store['a']

        tm.assert_series_equal(result, ser)
        assert result.index.freq == ser.index.freq
        tm.assert_class_equal(result.index, ser.index, obj="series index")

        idx = tm.makePeriodIndex(10)
        ser = Series(np.random.randn(len(idx)), idx)
        store['a'] = ser
        result = store['a']

        tm.assert_series_equal(result, ser)
        assert result.index.freq == ser.index.freq
        tm.assert_class_equal(result.index, ser.index, obj="series index")


def test_tseries_indices_frame():
    with ensure_clean_store() as store:
        idx = tm.makeDateIndex(10)
        df = DataFrame(np.random.randn(len(idx), 3), index=idx)
        store['a'] = df
        result = store['a']

        assert_frame_equal(result, df)
        assert result.index.freq == df.index.freq
        tm.assert_class_equal(result.index, df.index,
                              obj="dataframe index")

        idx = tm.makePeriodIndex(10)
        df = DataFrame(np.random.randn(len(idx), 3), idx)
        store['a'] = df
        result = store['a']

        assert_frame_equal(result, df)
        assert result.index.freq == df.index.freq
        tm.assert_class_equal(result.index, df.index,
                              obj="dataframe index")


def test_store_datetime_mixed():
    df = DataFrame(
        {'a': [1, 2, 3], 'b': [1., 2., 3.], 'c': ['a', 'b', 'c']})
    ts = tm.makeTimeSeries()
    df['d'] = ts.index[:3]
    _check_roundtrip(df, tm.assert_frame_equal)


def test_preserve_timedeltaindex_type():
    # GH9635
    # Storing TimedeltaIndexed DataFrames in fixed stores did not preserve
    # the type of the index.
    df = DataFrame(np.random.normal(size=(10, 5)))
    df.index = timedelta_range(
        start='0s', periods=10, freq='1s', name='example')

    with ensure_clean_store() as store:
        store['df'] = df
        assert_frame_equal(store['df'], df)
