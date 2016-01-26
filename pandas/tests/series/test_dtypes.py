# coding=utf-8
# pylint: disable-msg=E1101,W0612

import sys
from datetime import datetime
import string

from numpy import nan
import numpy as np

from pandas import Series
from pandas.tseries.index import Timestamp
from pandas.tseries.tdi import Timedelta

from pandas.compat import lrange, range, u
from pandas import compat
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesDtypes(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

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

    def test_astype_cast_nan_int(self):
        df = Series([1.0, 2.0, 3.0, np.nan])
        self.assertRaises(ValueError, df.astype, np.int64)

    def test_astype_cast_object_int(self):
        arr = Series(["car", "house", "tree", "1"])

        self.assertRaises(ValueError, arr.astype, int)
        self.assertRaises(ValueError, arr.astype, np.int64)
        self.assertRaises(ValueError, arr.astype, np.int8)

        arr = Series(['1', '2', '3', '4'], dtype=object)
        result = arr.astype(int)
        self.assert_numpy_array_equal(result, np.arange(1, 5))

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

    def test_complexx(self):

        # GH4819
        # complex access for ndarray compat
        a = np.arange(5)
        b = Series(a + 4j * a)
        tm.assert_almost_equal(a, b.real)
        tm.assert_almost_equal(4 * a, b.imag)

        b.real = np.arange(5) + 5
        tm.assert_almost_equal(a + 5, b.real)
        tm.assert_almost_equal(4 * a, b.imag)
