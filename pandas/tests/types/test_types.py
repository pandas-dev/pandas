# -*- coding: utf-8 -*-
import nose
import numpy as np

from pandas import NaT
from pandas.types.api import (DatetimeTZDtype, CategoricalDtype,
                              na_value_for_dtype, pandas_dtype)


def test_pandas_dtype():

    assert pandas_dtype('datetime64[ns, US/Eastern]') == DatetimeTZDtype(
        'datetime64[ns, US/Eastern]')
    assert pandas_dtype('category') == CategoricalDtype()
    for dtype in ['M8[ns]', 'm8[ns]', 'object', 'float64', 'int64']:
        assert pandas_dtype(dtype) == np.dtype(dtype)


def test_na_value_for_dtype():
    for dtype in [np.dtype('M8[ns]'), np.dtype('m8[ns]'),
                  DatetimeTZDtype('datetime64[ns, US/Eastern]')]:
        assert na_value_for_dtype(dtype) is NaT

    for dtype in ['u1', 'u2', 'u4', 'u8',
                  'i1', 'i2', 'i4', 'i8']:
        assert na_value_for_dtype(np.dtype(dtype)) == 0

    for dtype in ['bool']:
        assert na_value_for_dtype(np.dtype(dtype)) is False

    for dtype in ['f2', 'f4', 'f8']:
        assert np.isnan(na_value_for_dtype(np.dtype(dtype)))

    for dtype in ['O']:
        assert np.isnan(na_value_for_dtype(np.dtype(dtype)))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
