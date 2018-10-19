# -*- coding: utf-8 -*-
import operator

import numpy as np
import pytest

from pandas import Series

import pandas as pd
import pandas.util.testing as tm


# ------------------------------------------------------------------
# Comparisons

class TestSeriesComparison(object):

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes(self, opname):
        # GH#15115
        ser = Series([1, 3, 2], index=range(3))
        const = 2

        result = getattr(ser, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes_empty(self, opname):
        # GH#15115 empty Series case
        ser = Series([1, 3, 2], index=range(3))
        empty = ser.iloc[:0]
        const = 2

        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.le, operator.lt,
                                    operator.ge, operator.gt])
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('baz', 'baz', 'baz')])
    def test_ser_cmp_result_names(self, names, op):
        # datetime64 dtype
        dti = pd.date_range('1949-06-07 03:00:00',
                            freq='H', periods=5, name=names[0])
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # datetime64tz dtype
        dti = dti.tz_localize('US/Central')
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # timedelta64 dtype
        tdi = dti - dti.shift(1)
        ser = Series(tdi).rename(names[1])
        result = op(ser, tdi)
        assert result.name == names[2]

        # categorical
        if op in [operator.eq, operator.ne]:
            # categorical dtype comparisons raise for inequalities
            cidx = tdi.astype('category')
            ser = Series(cidx).rename(names[1])
            result = op(ser, cidx)
            assert result.name == names[2]


def test_pow_ops_object():
    # 22922
    # pow is weird with masking & 1, so testing here
    a = Series([1, np.nan, 1, np.nan], dtype=object)
    b = Series([1, np.nan, np.nan, 1], dtype=object)
    result = a ** b
    expected = Series(a.values ** b.values, dtype=object)
    tm.assert_series_equal(result, expected)

    result = b ** a
    expected = Series(b.values ** a.values, dtype=object)

    tm.assert_series_equal(result, expected)
