import os
import locale
import codecs
import nose

import numpy as np
from numpy import iinfo

import pandas as pd
from pandas import (date_range, Index, _np_version_under1p9)
import pandas.util.testing as tm
from pandas.tools.util import cartesian_product, to_numeric

CURRENT_LOCALE = locale.getlocale()
LOCALE_OVERRIDE = os.environ.get('LOCALE_OVERRIDE', None)


class TestCartesianProduct(tm.TestCase):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result1, result2 = cartesian_product([x, y])
        expected1 = np.array(['A', 'A', 'B', 'B', 'C', 'C'])
        expected2 = np.array([1, 22, 1, 22, 1, 22])
        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    def test_datetimeindex(self):
        # regression test for GitHub issue #6439
        # make sure that the ordering on datetimeindex is consistent
        x = date_range('2000-01-01', periods=2)
        result1, result2 = [Index(y).day for y in cartesian_product([x, x])]
        expected1 = np.array([1, 1, 2, 2], dtype=np.int32)
        expected2 = np.array([1, 2, 1, 2], dtype=np.int32)
        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    def test_empty(self):
        # product of empty factors
        X = [[], [0, 1], []]
        Y = [[], [], ['a', 'b', 'c']]
        for x, y in zip(X, Y):
            expected1 = np.array([], dtype=np.asarray(x).dtype)
            expected2 = np.array([], dtype=np.asarray(y).dtype)
            result1, result2 = cartesian_product([x, y])
            tm.assert_numpy_array_equal(result1, expected1)
            tm.assert_numpy_array_equal(result2, expected2)

        # empty product (empty input):
        result = cartesian_product([])
        expected = []
        tm.assert_equal(result, expected)

    def test_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        msg = "Input must be a list-like of list-likes"
        for X in invalid_inputs:
            tm.assertRaisesRegexp(TypeError, msg, cartesian_product, X=X)


class TestLocaleUtils(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestLocaleUtils, cls).setUpClass()
        cls.locales = tm.get_locales()

        if not cls.locales:
            raise nose.SkipTest("No locales found")

        tm._skip_if_windows()

    @classmethod
    def tearDownClass(cls):
        super(TestLocaleUtils, cls).tearDownClass()
        del cls.locales

    def test_get_locales(self):
        # all systems should have at least a single locale
        assert len(tm.get_locales()) > 0

    def test_get_locales_prefix(self):
        if len(self.locales) == 1:
            raise nose.SkipTest("Only a single locale found, no point in "
                                "trying to test filtering locale prefixes")
        first_locale = self.locales[0]
        assert len(tm.get_locales(prefix=first_locale[:2])) > 0

    def test_set_locale(self):
        if len(self.locales) == 1:
            raise nose.SkipTest("Only a single locale found, no point in "
                                "trying to test setting another locale")

        if LOCALE_OVERRIDE is None:
            lang, enc = 'it_CH', 'UTF-8'
        elif LOCALE_OVERRIDE == 'C':
            lang, enc = 'en_US', 'ascii'
        else:
            lang, enc = LOCALE_OVERRIDE.split('.')

        enc = codecs.lookup(enc).name
        new_locale = lang, enc

        if not tm._can_set_locale(new_locale):
            with tm.assertRaises(locale.Error):
                with tm.set_locale(new_locale):
                    pass
        else:
            with tm.set_locale(new_locale) as normalized_locale:
                new_lang, new_enc = normalized_locale.split('.')
                new_enc = codecs.lookup(enc).name
                normalized_locale = new_lang, new_enc
                self.assertEqual(normalized_locale, new_locale)

        current_locale = locale.getlocale()
        self.assertEqual(current_locale, CURRENT_LOCALE)


class TestToNumeric(tm.TestCase):

    def test_series(self):
        s = pd.Series(['1', '-3.14', '7'])
        res = to_numeric(s)
        expected = pd.Series([1, -3.14, 7])
        tm.assert_series_equal(res, expected)

        s = pd.Series(['1', '-3.14', 7])
        res = to_numeric(s)
        tm.assert_series_equal(res, expected)

    def test_series_numeric(self):
        s = pd.Series([1, 3, 4, 5], index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

        s = pd.Series([1., 3., 4., 5.], index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

        # bool is regarded as numeric
        s = pd.Series([True, False, True, True],
                      index=list('ABCD'), name='XXX')
        res = to_numeric(s)
        tm.assert_series_equal(res, s)

    def test_error(self):
        s = pd.Series([1, -3.14, 'apple'])
        msg = 'Unable to parse string "apple" at position 2'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([1, -3.14, 'apple'])
        tm.assert_series_equal(res, expected)

        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1, -3.14, np.nan])
        tm.assert_series_equal(res, expected)

        s = pd.Series(['orange', 1, -3.14, 'apple'])
        msg = 'Unable to parse string "orange" at position 0'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

    def test_error_seen_bool(self):
        s = pd.Series([True, False, 'apple'])
        msg = 'Unable to parse string "apple" at position 2'
        with tm.assertRaisesRegexp(ValueError, msg):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([True, False, 'apple'])
        tm.assert_series_equal(res, expected)

        # coerces to float
        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1., 0., np.nan])
        tm.assert_series_equal(res, expected)

    def test_list(self):
        s = ['1', '-3.14', '7']
        res = to_numeric(s)
        expected = np.array([1, -3.14, 7])
        tm.assert_numpy_array_equal(res, expected)

    def test_list_numeric(self):
        s = [1, 3, 4, 5]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s, dtype=np.int64))

        s = [1., 3., 4., 5.]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s))

        # bool is regarded as numeric
        s = [True, False, True, True]
        res = to_numeric(s)
        tm.assert_numpy_array_equal(res, np.array(s))

    def test_numeric(self):
        s = pd.Series([1, -3.14, 7], dtype='O')
        res = to_numeric(s)
        expected = pd.Series([1, -3.14, 7])
        tm.assert_series_equal(res, expected)

        s = pd.Series([1, -3.14, 7])
        res = to_numeric(s)
        tm.assert_series_equal(res, expected)

    def test_all_nan(self):
        s = pd.Series(['a', 'b', 'c'])
        res = to_numeric(s, errors='coerce')
        expected = pd.Series([np.nan, np.nan, np.nan])
        tm.assert_series_equal(res, expected)

    def test_type_check(self):
        # GH 11776
        df = pd.DataFrame({'a': [1, -3.14, 7], 'b': ['4', '5', '6']})
        with tm.assertRaisesRegexp(TypeError, "1-d array"):
            to_numeric(df)
        for errors in ['ignore', 'raise', 'coerce']:
            with tm.assertRaisesRegexp(TypeError, "1-d array"):
                to_numeric(df, errors=errors)

    def test_scalar(self):
        self.assertEqual(pd.to_numeric(1), 1)
        self.assertEqual(pd.to_numeric(1.1), 1.1)

        self.assertEqual(pd.to_numeric('1'), 1)
        self.assertEqual(pd.to_numeric('1.1'), 1.1)

        with tm.assertRaises(ValueError):
            to_numeric('XX', errors='raise')

        self.assertEqual(to_numeric('XX', errors='ignore'), 'XX')
        self.assertTrue(np.isnan(to_numeric('XX', errors='coerce')))

    def test_numeric_dtypes(self):
        idx = pd.Index([1, 2, 3], name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, idx)

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.values)

        idx = pd.Index([1., np.nan, 3., np.nan], name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, idx)

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.values)

    def test_str(self):
        idx = pd.Index(['1', '2', '3'], name='xxx')
        exp = np.array([1, 2, 3], dtype='int64')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(exp, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(exp, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, exp)

        idx = pd.Index(['1.5', '2.7', '3.4'], name='xxx')
        exp = np.array([1.5, 2.7, 3.4])
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(exp, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(exp, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, exp)

    def test_datetimelike(self):
        for tz in [None, 'US/Eastern', 'Asia/Tokyo']:
            idx = pd.date_range('20130101', periods=3, tz=tz, name='xxx')
            res = pd.to_numeric(idx)
            tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

            res = pd.to_numeric(pd.Series(idx, name='xxx'))
            tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

            res = pd.to_numeric(idx.values)
            tm.assert_numpy_array_equal(res, idx.asi8)

    def test_timedelta(self):
        idx = pd.timedelta_range('1 days', periods=3, freq='D', name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

        res = pd.to_numeric(pd.Series(idx, name='xxx'))
        tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

        res = pd.to_numeric(idx.values)
        tm.assert_numpy_array_equal(res, idx.asi8)

    def test_period(self):
        idx = pd.period_range('2011-01', periods=3, freq='M', name='xxx')
        res = pd.to_numeric(idx)
        tm.assert_index_equal(res, pd.Index(idx.asi8, name='xxx'))

        # ToDo: enable when we can support native PeriodDtype
        # res = pd.to_numeric(pd.Series(idx, name='xxx'))
        # tm.assert_series_equal(res, pd.Series(idx.asi8, name='xxx'))

    def test_non_hashable(self):
        # Test for Bug #13324
        s = pd.Series([[10.0, 2], 1.0, 'apple'])
        res = pd.to_numeric(s, errors='coerce')
        tm.assert_series_equal(res, pd.Series([np.nan, 1.0, np.nan]))

        res = pd.to_numeric(s, errors='ignore')
        tm.assert_series_equal(res, pd.Series([[10.0, 2], 1.0, 'apple']))

        with self.assertRaisesRegexp(TypeError, "Invalid object type"):
            pd.to_numeric(s)

    def test_downcast(self):
        # see gh-13352
        mixed_data = ['1', 2, 3]
        int_data = [1, 2, 3]
        date_data = np.array(['1970-01-02', '1970-01-03',
                              '1970-01-04'], dtype='datetime64[D]')

        invalid_downcast = 'unsigned-integer'
        msg = 'invalid downcasting method provided'

        smallest_int_dtype = np.dtype(np.typecodes['Integer'][0])
        smallest_uint_dtype = np.dtype(np.typecodes['UnsignedInteger'][0])

        # support below np.float32 is rare and far between
        float_32_char = np.dtype(np.float32).char
        smallest_float_dtype = float_32_char

        for data in (mixed_data, int_data, date_data):
            with self.assertRaisesRegexp(ValueError, msg):
                pd.to_numeric(data, downcast=invalid_downcast)

            expected = np.array([1, 2, 3], dtype=np.int64)

            res = pd.to_numeric(data)
            tm.assert_numpy_array_equal(res, expected)

            res = pd.to_numeric(data, downcast=None)
            tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_int_dtype)

            for signed_downcast in ('integer', 'signed'):
                res = pd.to_numeric(data, downcast=signed_downcast)
                tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_uint_dtype)
            res = pd.to_numeric(data, downcast='unsigned')
            tm.assert_numpy_array_equal(res, expected)

            expected = np.array([1, 2, 3], dtype=smallest_float_dtype)
            res = pd.to_numeric(data, downcast='float')
            tm.assert_numpy_array_equal(res, expected)

        # if we can't successfully cast the given
        # data to a numeric dtype, do not bother
        # with the downcast parameter
        data = ['foo', 2, 3]
        expected = np.array(data, dtype=object)
        res = pd.to_numeric(data, errors='ignore',
                            downcast='unsigned')
        tm.assert_numpy_array_equal(res, expected)

        # cannot cast to an unsigned integer because
        # we have a negative number
        data = ['-1', 2, 3]
        expected = np.array([-1, 2, 3], dtype=np.int64)
        res = pd.to_numeric(data, downcast='unsigned')
        tm.assert_numpy_array_equal(res, expected)

        # cannot cast to an integer (signed or unsigned)
        # because we have a float number
        data = ['1.1', 2, 3]
        expected = np.array([1.1, 2, 3], dtype=np.float64)

        for downcast in ('integer', 'signed', 'unsigned'):
            res = pd.to_numeric(data, downcast=downcast)
            tm.assert_numpy_array_equal(res, expected)

        # the smallest integer dtype need not be np.(u)int8
        data = ['256', 257, 258]

        for downcast, expected_dtype in zip(
                ['integer', 'signed', 'unsigned'],
                [np.int16, np.int16, np.uint16]):
            expected = np.array([256, 257, 258], dtype=expected_dtype)
            res = pd.to_numeric(data, downcast=downcast)
            tm.assert_numpy_array_equal(res, expected)

    def test_downcast_limits(self):
        # Test the limits of each downcast. Bug: #14401.
        # Check to make sure numpy is new enough to run this test.
        if _np_version_under1p9:
            raise nose.SkipTest("Numpy version is under 1.9")

        i = 'integer'
        u = 'unsigned'
        dtype_downcast_min_max = [
            ('int8', i, [iinfo(np.int8).min, iinfo(np.int8).max]),
            ('int16', i, [iinfo(np.int16).min, iinfo(np.int16).max]),
            ('int32', i, [iinfo(np.int32).min, iinfo(np.int32).max]),
            ('int64', i, [iinfo(np.int64).min, iinfo(np.int64).max]),
            ('uint8', u, [iinfo(np.uint8).min, iinfo(np.uint8).max]),
            ('uint16', u, [iinfo(np.uint16).min, iinfo(np.uint16).max]),
            ('uint32', u, [iinfo(np.uint32).min, iinfo(np.uint32).max]),
            # Test will be skipped until there is more uint64 support.
            # ('uint64', u, [iinfo(uint64).min, iinfo(uint64).max]),
            ('int16', i, [iinfo(np.int8).min, iinfo(np.int8).max + 1]),
            ('int32', i, [iinfo(np.int16).min, iinfo(np.int16).max + 1]),
            ('int64', i, [iinfo(np.int32).min, iinfo(np.int32).max + 1]),
            ('int16', i, [iinfo(np.int8).min - 1, iinfo(np.int16).max]),
            ('int32', i, [iinfo(np.int16).min - 1, iinfo(np.int32).max]),
            ('int64', i, [iinfo(np.int32).min - 1, iinfo(np.int64).max]),
            ('uint16', u, [iinfo(np.uint8).min, iinfo(np.uint8).max + 1]),
            ('uint32', u, [iinfo(np.uint16).min, iinfo(np.uint16).max + 1]),
            # Test will be skipped until there is more uint64 support.
            # ('uint64', u, [iinfo(np.uint32).min, iinfo(np.uint32).max + 1]),
        ]

        for dtype, downcast, min_max in dtype_downcast_min_max:
            series = pd.to_numeric(pd.Series(min_max), downcast=downcast)
            tm.assert_equal(series.dtype, dtype)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
