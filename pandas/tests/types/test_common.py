# -*- coding: utf-8 -*-

import nose
import numpy as np

from pandas.types.dtypes import DatetimeTZDtype, PeriodDtype, CategoricalDtype
from pandas.types.common import pandas_dtype, is_dtype_equal

_multiprocess_can_split_ = True


def test_pandas_dtype():

    assert pandas_dtype('datetime64[ns, US/Eastern]') == DatetimeTZDtype(
        'datetime64[ns, US/Eastern]')
    assert pandas_dtype('category') == CategoricalDtype()
    for dtype in ['M8[ns]', 'm8[ns]', 'object', 'float64', 'int64']:
        assert pandas_dtype(dtype) == np.dtype(dtype)


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
