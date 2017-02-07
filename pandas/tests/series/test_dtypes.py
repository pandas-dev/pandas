# coding=utf-8
# pylint: disable-msg=E1101,W0612

import sys
from datetime import datetime
import string

from numpy import nan
import numpy as np

from pandas import Series, Timestamp, Timedelta, DataFrame, date_range

from pandas.compat import lrange, range, u
from pandas import compat
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesDtypes(TestData, tm.TestCase):

    def test_astype(self):
        s = Series(np.random.randn(5), name='foo')

        for dtype in ['float32', 'float64', 'int64', 'int32']:
            astyped = s.astype(dtype)
            self.assertEqual(astyped.dtype, dtype)
            self.assertEqual(astyped.name, s.name)

    def test_dtype(self):

        self.assertEqual(self.ts.dtype, np.dtype('float64'))
        self.assertEqual(self.ts.dtypes, np.dtype('float64'))
        self.assertEqual(self.ts.ftype, 'float64:dense')
        self.assertEqual(self.ts.ftypes, 'float64:dense')
        assert_series_equal(self.ts.get_dtype_counts(), Series(1, ['float64']))
        assert_series_equal(self.ts.get_ftype_counts(), Series(
            1, ['float64:dense']))

    def test_astype_cast_nan_inf_int(self):
        # GH14265, check nan and inf raise error when converting to int
        types = [np.int32, np.int64]
        values = [np.nan, np.inf]
        msg = 'Cannot convert non-finite values \\(NA or inf\\) to integer'

        for this_type in types:
            for this_val in values:
                s = Series([this_val])
                with self.assertRaisesRegexp(ValueError, msg):
                    s.astype(this_type)

    def test_astype_cast_object_int(self):
        arr = Series(["car", "house", "tree", "1"])

        self.assertRaises(ValueError, arr.astype, int)
        self.assertRaises(ValueError, arr.astype, np.int64)
        self.assertRaises(ValueError, arr.astype, np.int8)

        arr = Series(['1', '2', '3', '4'], dtype=object)
        result = arr.astype(int)
        self.assert_series_equal(result, Series(np.arange(1, 5)))

    def test_astype_datetimes(self):
        import pandas.tslib as tslib

        s = Series(tslib.iNaT, dtype='M8[ns]', index=lrange(5))
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

        s = Series([datetime(2001, 1, 2, 0, 0)])
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

        s = Series([datetime(2001, 1, 2, 0, 0) for i in range(3)])
        s[1] = np.nan
        self.assertEqual(s.dtype, 'M8[ns]')
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

    def test_astype_str(self):
        # GH4405
        digits = string.digits
        s1 = Series([digits * 10, tm.rands(63), tm.rands(64), tm.rands(1000)])
        s2 = Series([digits * 10, tm.rands(63), tm.rands(64), nan, 1.0])
        types = (compat.text_type, np.str_)
        for typ in types:
            for s in (s1, s2):
                res = s.astype(typ)
                expec = s.map(compat.text_type)
                assert_series_equal(res, expec)

        # GH9757
        # Test str and unicode on python 2.x and just str on python 3.x
        for tt in set([str, compat.text_type]):
            ts = Series([Timestamp('2010-01-04 00:00:00')])
            s = ts.astype(tt)
            expected = Series([tt('2010-01-04')])
            assert_series_equal(s, expected)

            ts = Series([Timestamp('2010-01-04 00:00:00', tz='US/Eastern')])
            s = ts.astype(tt)
            expected = Series([tt('2010-01-04 00:00:00-05:00')])
            assert_series_equal(s, expected)

            td = Series([Timedelta(1, unit='d')])
            s = td.astype(tt)
            expected = Series([tt('1 days 00:00:00.000000000')])
            assert_series_equal(s, expected)

    def test_astype_unicode(self):

        # GH7758
        # a bit of magic is required to set default encoding encoding to utf-8
        digits = string.digits
        test_series = [
            Series([digits * 10, tm.rands(63), tm.rands(64), tm.rands(1000)]),
            Series([u('データーサイエンス、お前はもう死んでいる')]),

        ]

        former_encoding = None
        if not compat.PY3:
            # in python we can force the default encoding for this test
            former_encoding = sys.getdefaultencoding()
            reload(sys)  # noqa
            sys.setdefaultencoding("utf-8")
        if sys.getdefaultencoding() == "utf-8":
            test_series.append(Series([u('野菜食べないとやばい')
                                       .encode("utf-8")]))
        for s in test_series:
            res = s.astype("unicode")
            expec = s.map(compat.text_type)
            assert_series_equal(res, expec)
        # restore the former encoding
        if former_encoding is not None and former_encoding != "utf-8":
            reload(sys)  # noqa
            sys.setdefaultencoding(former_encoding)

    def test_astype_dict(self):
        # GH7271
        s = Series(range(0, 10, 2), name='abc')

        result = s.astype({'abc': str})
        expected = Series(['0', '2', '4', '6', '8'], name='abc')
        assert_series_equal(result, expected)

        result = s.astype({'abc': 'float64'})
        expected = Series([0.0, 2.0, 4.0, 6.0, 8.0], dtype='float64',
                          name='abc')
        assert_series_equal(result, expected)

        self.assertRaises(KeyError, s.astype, {'abc': str, 'def': str})
        self.assertRaises(KeyError, s.astype, {0: str})

    def test_complexx(self):
        # GH4819
        # complex access for ndarray compat
        a = np.arange(5, dtype=np.float64)
        b = Series(a + 4j * a)
        tm.assert_numpy_array_equal(a, b.real)
        tm.assert_numpy_array_equal(4 * a, b.imag)

        b.real = np.arange(5) + 5
        tm.assert_numpy_array_equal(a + 5, b.real)
        tm.assert_numpy_array_equal(4 * a, b.imag)

    def test_arg_for_errors_in_astype(self):
        # issue #14878

        sr = Series([1, 2, 3])

        with self.assertRaises(ValueError):
            sr.astype(np.float64, errors=False)

        with tm.assert_produces_warning(FutureWarning):
            sr.astype(np.int8, raise_on_error=True)

        sr.astype(np.int8, errors='raise')

    def test_intercept_astype_object(self):
        series = Series(date_range('1/1/2000', periods=10))

        # this test no longer makes sense as series is by default already
        # M8[ns]
        expected = series.astype('object')

        df = DataFrame({'a': series,
                        'b': np.random.randn(len(series))})
        exp_dtypes = Series([np.dtype('datetime64[ns]'),
                             np.dtype('float64')], index=['a', 'b'])
        tm.assert_series_equal(df.dtypes, exp_dtypes)

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())

        df = DataFrame({'a': series, 'b': ['foo'] * len(series)})

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())
