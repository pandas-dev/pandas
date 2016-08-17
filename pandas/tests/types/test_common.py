# -*- coding: utf-8 -*-

import nose
import numpy as np

from pandas.types.dtypes import DatetimeTZDtype, PeriodDtype, CategoricalDtype
from pandas.types.common import pandas_dtype, is_dtype_equal

import pandas.util.testing as tm

_multiprocess_can_split_ = True


class TestPandasDtype(tm.TestCase):

    def test_numpy_dtype(self):
        for dtype in ['M8[ns]', 'm8[ns]', 'object', 'float64', 'int64']:
            self.assertEqual(pandas_dtype(dtype), np.dtype(dtype))

    def test_numpy_string_dtype(self):
        # do not parse freq-like string as period dtype
        self.assertEqual(pandas_dtype('U'), np.dtype('U'))
        self.assertEqual(pandas_dtype('S'), np.dtype('S'))

    def test_datetimetz_dtype(self):
        for dtype in ['datetime64[ns, US/Eastern]',
                      'datetime64[ns, Asia/Tokyo]',
                      'datetime64[ns, UTC]']:
            self.assertIs(pandas_dtype(dtype), DatetimeTZDtype(dtype))
            self.assertEqual(pandas_dtype(dtype), DatetimeTZDtype(dtype))
            self.assertEqual(pandas_dtype(dtype), dtype)

    def test_categorical_dtype(self):
        self.assertEqual(pandas_dtype('category'), CategoricalDtype())

    def test_period_dtype(self):
        for dtype in ['period[D]', 'period[3M]', 'period[U]',
                      'Period[D]', 'Period[3M]', 'Period[U]']:
            self.assertIs(pandas_dtype(dtype), PeriodDtype(dtype))
            self.assertEqual(pandas_dtype(dtype), PeriodDtype(dtype))
            self.assertEqual(pandas_dtype(dtype), dtype)


def test_dtype_equal():
    assert is_dtype_equal(np.int64, np.int64)
    assert not is_dtype_equal(np.int64, np.float64)

    p1 = PeriodDtype('D')
    p2 = PeriodDtype('D')
    assert is_dtype_equal(p1, p2)
    assert not is_dtype_equal(np.int64, p1)

    p3 = PeriodDtype('2D')
    assert not is_dtype_equal(p1, p3)

    assert not DatetimeTZDtype.is_dtype(np.int64)
    assert not PeriodDtype.is_dtype(np.int64)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
