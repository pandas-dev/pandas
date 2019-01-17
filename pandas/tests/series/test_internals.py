# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import NaT, Series, Timestamp
from pandas.core.internals.blocks import IntBlock
import pandas.util.testing as tm
from pandas.util.testing import assert_series_equal


class TestSeriesInternals(object):

    def test_convert_objects(self):

        s = Series([1., 2, 3], index=['a', 'b', 'c'])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates=False,
                                       convert_numeric=True)
        assert_series_equal(result, s)

        # force numeric conversion
        r = s.copy().astype('O')
        r['a'] = '1'
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False,
                                       convert_numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = '1.'
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False,
                                       convert_numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = 'garbled'
        expected = s.copy()
        expected['a'] = np.nan
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False,
                                       convert_numeric=True)
        assert_series_equal(result, expected)

        # GH 4119, not converting a mixed type (e.g.floats and object)
        s = Series([1, 'na', 3, 4])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_numeric=True)
        expected = Series([1, np.nan, 3, 4])
        assert_series_equal(result, expected)

        s = Series([1, '', 3, 4])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_numeric=True)
        expected = Series([1, np.nan, 3, 4])
        assert_series_equal(result, expected)

        # dates
        s = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0),
                    datetime(2001, 1, 3, 0, 0)])
        s2 = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0),
                     datetime(2001, 1, 3, 0, 0), 'foo', 1.0, 1,
                     Timestamp('20010104'), '20010105'],
                    dtype='O')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates=True,
                                       convert_numeric=False)
        expected = Series([Timestamp('20010101'), Timestamp('20010102'),
                           Timestamp('20010103')], dtype='M8[ns]')
        assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=True)
        assert_series_equal(result, expected)

        expected = Series([Timestamp('20010101'), Timestamp('20010102'),
                           Timestamp('20010103'),
                           NaT, NaT, NaT, Timestamp('20010104'),
                           Timestamp('20010105')], dtype='M8[ns]')
        with tm.assert_produces_warning(FutureWarning):
            result = s2.convert_objects(convert_dates='coerce',
                                        convert_numeric=False)
        assert_series_equal(result, expected)
        with tm.assert_produces_warning(FutureWarning):
            result = s2.convert_objects(convert_dates='coerce',
                                        convert_numeric=True)
        assert_series_equal(result, expected)

        # preserver all-nans (if convert_dates='coerce')
        s = Series(['foo', 'bar', 1, 1.0], dtype='O')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        expected = Series([NaT] * 2 + [Timestamp(1)] * 2)
        assert_series_equal(result, expected)

        # preserver if non-object
        s = Series([1], dtype='float32')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        assert_series_equal(result, s)

        # r = s.copy()
        # r[0] = np.nan
        # result = r.convert_objects(convert_dates=True,convert_numeric=False)
        # assert result.dtype == 'M8[ns]'

        # dateutil parses some single letters into today's value as a date
        for x in 'abcdefghijklmnopqrstuvwxyz':
            s = Series([x])
            with tm.assert_produces_warning(FutureWarning):
                result = s.convert_objects(convert_dates='coerce')
            assert_series_equal(result, s)
            s = Series([x.upper()])
            with tm.assert_produces_warning(FutureWarning):
                result = s.convert_objects(convert_dates='coerce')
            assert_series_equal(result, s)

    def test_convert_objects_preserve_bool(self):
        s = Series([1, True, 3, 5], dtype=object)
        with tm.assert_produces_warning(FutureWarning):
            r = s.convert_objects(convert_numeric=True)
        e = Series([1, 1, 3, 5], dtype='i8')
        tm.assert_series_equal(r, e)

    def test_convert_objects_preserve_all_bool(self):
        s = Series([False, True, False, False], dtype=object)
        with tm.assert_produces_warning(FutureWarning):
            r = s.convert_objects(convert_numeric=True)
        e = Series([False, True, False, False], dtype=bool)
        tm.assert_series_equal(r, e)

    # GH 10265
    def test_convert(self):
        # Tests: All to nans, coerce, true
        # Test coercion returns correct type
        s = Series(['a', 'b', 'c'])
        results = s._convert(datetime=True, coerce=True)
        expected = Series([NaT] * 3)
        assert_series_equal(results, expected)

        results = s._convert(numeric=True, coerce=True)
        expected = Series([np.nan] * 3)
        assert_series_equal(results, expected)

        expected = Series([NaT] * 3, dtype=np.dtype('m8[ns]'))
        results = s._convert(timedelta=True, coerce=True)
        assert_series_equal(results, expected)

        dt = datetime(2001, 1, 1, 0, 0)
        td = dt - datetime(2000, 1, 1, 0, 0)

        # Test coercion with mixed types
        s = Series(['a', '3.1415', dt, td])
        results = s._convert(datetime=True, coerce=True)
        expected = Series([NaT, NaT, dt, NaT])
        assert_series_equal(results, expected)

        results = s._convert(numeric=True, coerce=True)
        expected = Series([np.nan, 3.1415, np.nan, np.nan])
        assert_series_equal(results, expected)

        results = s._convert(timedelta=True, coerce=True)
        expected = Series([NaT, NaT, NaT, td],
                          dtype=np.dtype('m8[ns]'))
        assert_series_equal(results, expected)

        # Test standard conversion returns original
        results = s._convert(datetime=True)
        assert_series_equal(results, s)
        results = s._convert(numeric=True)
        expected = Series([np.nan, 3.1415, np.nan, np.nan])
        assert_series_equal(results, expected)
        results = s._convert(timedelta=True)
        assert_series_equal(results, s)

        # test pass-through and non-conversion when other types selected
        s = Series(['1.0', '2.0', '3.0'])
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([1.0, 2.0, 3.0])
        assert_series_equal(results, expected)
        results = s._convert(True, False, True)
        assert_series_equal(results, s)

        s = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 1, 0, 0)],
                   dtype='O')
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 1, 0,
                                                                0)])
        assert_series_equal(results, expected)
        results = s._convert(datetime=False, numeric=True, timedelta=True)
        assert_series_equal(results, s)

        td = datetime(2001, 1, 1, 0, 0) - datetime(2000, 1, 1, 0, 0)
        s = Series([td, td], dtype='O')
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([td, td])
        assert_series_equal(results, expected)
        results = s._convert(True, True, False)
        assert_series_equal(results, s)

        s = Series([1., 2, 3], index=['a', 'b', 'c'])
        result = s._convert(numeric=True)
        assert_series_equal(result, s)

        # force numeric conversion
        r = s.copy().astype('O')
        r['a'] = '1'
        result = r._convert(numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = '1.'
        result = r._convert(numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = 'garbled'
        result = r._convert(numeric=True)
        expected = s.copy()
        expected['a'] = np.nan
        assert_series_equal(result, expected)

        # GH 4119, not converting a mixed type (e.g.floats and object)
        s = Series([1, 'na', 3, 4])
        result = s._convert(datetime=True, numeric=True)
        expected = Series([1, np.nan, 3, 4])
        assert_series_equal(result, expected)

        s = Series([1, '', 3, 4])
        result = s._convert(datetime=True, numeric=True)
        assert_series_equal(result, expected)

        # dates
        s = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0),
                    datetime(2001, 1, 3, 0, 0)])
        s2 = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0),
                     datetime(2001, 1, 3, 0, 0), 'foo', 1.0, 1,
                     Timestamp('20010104'), '20010105'], dtype='O')

        result = s._convert(datetime=True)
        expected = Series([Timestamp('20010101'), Timestamp('20010102'),
                           Timestamp('20010103')], dtype='M8[ns]')
        assert_series_equal(result, expected)

        result = s._convert(datetime=True, coerce=True)
        assert_series_equal(result, expected)

        expected = Series([Timestamp('20010101'), Timestamp('20010102'),
                           Timestamp('20010103'), NaT, NaT, NaT,
                           Timestamp('20010104'), Timestamp('20010105')],
                          dtype='M8[ns]')
        result = s2._convert(datetime=True, numeric=False, timedelta=False,
                             coerce=True)
        assert_series_equal(result, expected)
        result = s2._convert(datetime=True, coerce=True)
        assert_series_equal(result, expected)

        s = Series(['foo', 'bar', 1, 1.0], dtype='O')
        result = s._convert(datetime=True, coerce=True)
        expected = Series([NaT] * 2 + [Timestamp(1)] * 2)
        assert_series_equal(result, expected)

        # preserver if non-object
        s = Series([1], dtype='float32')
        result = s._convert(datetime=True, coerce=True)
        assert_series_equal(result, s)

        # r = s.copy()
        # r[0] = np.nan
        # result = r._convert(convert_dates=True,convert_numeric=False)
        # assert result.dtype == 'M8[ns]'

        # dateutil parses some single letters into today's value as a date
        expected = Series([NaT])
        for x in 'abcdefghijklmnopqrstuvwxyz':
            s = Series([x])
            result = s._convert(datetime=True, coerce=True)
            assert_series_equal(result, expected)
            s = Series([x.upper()])
            result = s._convert(datetime=True, coerce=True)
            assert_series_equal(result, expected)

    def test_convert_no_arg_error(self):
        s = Series(['1.0', '2'])
        msg = r"At least one of datetime, numeric or timedelta must be True\."
        with pytest.raises(ValueError, match=msg):
            s._convert()

    def test_convert_preserve_bool(self):
        s = Series([1, True, 3, 5], dtype=object)
        r = s._convert(datetime=True, numeric=True)
        e = Series([1, 1, 3, 5], dtype='i8')
        tm.assert_series_equal(r, e)

    def test_convert_preserve_all_bool(self):
        s = Series([False, True, False, False], dtype=object)
        r = s._convert(datetime=True, numeric=True)
        e = Series([False, True, False, False], dtype=bool)
        tm.assert_series_equal(r, e)

    def test_constructor_no_pandas_array(self):
        ser = pd.Series([1, 2, 3])
        result = pd.Series(ser.array)
        tm.assert_series_equal(ser, result)
        assert isinstance(result._data.blocks[0], IntBlock)

    def test_from_array(self):
        result = pd.Series(pd.array(['1H', '2H'], dtype='timedelta64[ns]'))
        assert result._data.blocks[0].is_extension is False

        result = pd.Series(pd.array(['2015'], dtype='datetime64[ns]'))
        assert result._data.blocks[0].is_extension is False

    def test_from_list_dtype(self):
        result = pd.Series(['1H', '2H'], dtype='timedelta64[ns]')
        assert result._data.blocks[0].is_extension is False

        result = pd.Series(['2015'], dtype='datetime64[ns]')
        assert result._data.blocks[0].is_extension is False


def test_hasnans_unchached_for_series():
    # GH#19700
    idx = pd.Index([0, 1])
    assert idx.hasnans is False
    assert 'hasnans' in idx._cache
    ser = idx.to_series()
    assert ser.hasnans is False
    assert not hasattr(ser, '_cache')
    ser.iloc[-1] = np.nan
    assert ser.hasnans is True
    assert Series.hasnans.__doc__ == pd.Index.hasnans.__doc__
