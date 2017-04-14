# -*- coding: utf-8 -*-

import pytest
import numpy as np

from pandas.core.typed.dtypes import (
    DatetimeTZDtype, PeriodDtype, CategoricalDtype)
from pandas.core.typed.common import (
    pandas_dtype, is_dtype_equal)

import pandas.util.testing as tm


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


dtypes = dict(datetime_tz=pandas_dtype('datetime64[ns, US/Eastern]'),
              datetime=pandas_dtype('datetime64[ns]'),
              timedelta=pandas_dtype('timedelta64[ns]'),
              period=PeriodDtype('D'),
              integer=np.dtype(np.int64),
              float=np.dtype(np.float64),
              object=np.dtype(np.object),
              category=pandas_dtype('category'))


@pytest.mark.parametrize('name1,dtype1',
                         list(dtypes.items()),
                         ids=lambda x: str(x))
@pytest.mark.parametrize('name2,dtype2',
                         list(dtypes.items()),
                         ids=lambda x: str(x))
def test_dtype_equal(name1, dtype1, name2, dtype2):

    # match equal to self, but not equal to other
    assert is_dtype_equal(dtype1, dtype1)
    if name1 != name2:
        assert not is_dtype_equal(dtype1, dtype2)


def test_dtype_equal_strict():

    # we are strict on kind equality
    for dtype in [np.int8, np.int16, np.int32]:
        assert not is_dtype_equal(np.int64, dtype)

    for dtype in [np.float32]:
        assert not is_dtype_equal(np.float64, dtype)

    # strict w.r.t. PeriodDtype
    assert not is_dtype_equal(PeriodDtype('D'),
                              PeriodDtype('2D'))

    # strict w.r.t. datetime64
    assert not is_dtype_equal(
        pandas_dtype('datetime64[ns, US/Eastern]'),
        pandas_dtype('datetime64[ns, CET]'))

    # see gh-15941: no exception should be raised
    assert not is_dtype_equal(None, None)


def get_is_dtype_funcs():
    """
    Get all functions in pandas.core.typed.common that
    begin with 'is_' and end with 'dtype'

    """
    import pandas.core.typed.common as com

    fnames = [f for f in dir(com) if (f.startswith('is_') and
                                      f.endswith('dtype'))]
    return [getattr(com, fname) for fname in fnames]


@pytest.mark.parametrize('func',
                         get_is_dtype_funcs(),
                         ids=lambda x: x.__name__)
def test_get_dtype_error_catch(func):
    # see gh-15941
    #
    # No exception should be raised.

    assert not func(None)
