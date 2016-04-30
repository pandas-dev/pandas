import os
import locale
import codecs
import nose

import numpy as np
from numpy.testing import assert_equal

import pandas as pd
from pandas import date_range, Index
import pandas.util.testing as tm
from pandas.tools.util import cartesian_product, to_numeric

CURRENT_LOCALE = locale.getlocale()
LOCALE_OVERRIDE = os.environ.get('LOCALE_OVERRIDE', None)


class TestCartesianProduct(tm.TestCase):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result = cartesian_product([x, y])
        expected = [np.array(['A', 'A', 'B', 'B', 'C', 'C']),
                    np.array([1, 22, 1, 22, 1, 22])]
        assert_equal(result, expected)

    def test_datetimeindex(self):
        # regression test for GitHub issue #6439
        # make sure that the ordering on datetimeindex is consistent
        x = date_range('2000-01-01', periods=2)
        result = [Index(y).day for y in cartesian_product([x, x])]
        expected = [np.array([1, 1, 2, 2]), np.array([1, 2, 1, 2])]
        assert_equal(result, expected)


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

        if LOCALE_OVERRIDE is not None:
            lang, enc = LOCALE_OVERRIDE.split('.')
        else:
            lang, enc = 'it_CH', 'UTF-8'

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
        with tm.assertRaises(ValueError):
            to_numeric(s, errors='raise')

        res = to_numeric(s, errors='ignore')
        expected = pd.Series([1, -3.14, 'apple'])
        tm.assert_series_equal(res, expected)

        res = to_numeric(s, errors='coerce')
        expected = pd.Series([1, -3.14, np.nan])
        tm.assert_series_equal(res, expected)

    def test_error_seen_bool(self):
        s = pd.Series([True, False, 'apple'])
        with tm.assertRaises(ValueError):
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
        tm.assert_numpy_array_equal(res, np.array(s))

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


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
