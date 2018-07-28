# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
from datetime import timedelta

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from pandas import Timedelta


# ------------------------------------------------------------------
# Numeric dtypes Arithmetic with Timedelta Scalar

class TestNumericArraylikeArithmeticWithTimedeltaScalar(object):

    @pytest.mark.parametrize('box', [
        pd.Index,
        pd.Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="block.eval incorrect",
                                             strict=True))
    ])
    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 11)),
        pd.UInt64Index(range(1, 11)),
        pd.Float64Index(range(1, 11)),
        pd.RangeIndex(1, 11)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_index_mul_timedelta(self, scalar_td, index, box):
        # GH#19333

        if (box is pd.Series and
                type(scalar_td) is timedelta and index.dtype == 'f8'):
            raise pytest.xfail(reason="Cannot multiply timedelta by float")

        expected = pd.timedelta_range('1 days', '10 days')

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = index * scalar_td
        tm.assert_equal(result, expected)

        commute = scalar_td * index
        tm.assert_equal(commute, expected)

    @pytest.mark.parametrize('box', [pd.Index, pd.Series, pd.DataFrame])
    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 3)),
        pd.UInt64Index(range(1, 3)),
        pd.Float64Index(range(1, 3)),
        pd.RangeIndex(1, 3)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_index_rdiv_timedelta(self, scalar_td, index, box):

        if box is pd.Series and type(scalar_td) is timedelta:
            raise pytest.xfail(reason="TODO: Figure out why this case fails")
        if box is pd.DataFrame and isinstance(scalar_td, timedelta):
            raise pytest.xfail(reason="TODO: Figure out why this case fails")

        expected = pd.TimedeltaIndex(['1 Day', '12 Hours'])

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = scalar_td / index
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            index / scalar_td


# ------------------------------------------------------------------
# Timedelta64[ns] dtype Arithmetic Operations


class TestTimedeltaArraylikeInvalidArithmeticOps(object):

    @pytest.mark.parametrize('box', [
        pd.Index,
        pd.Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="raises ValueError "
                                                    "instead of TypeError",
                                             strict=True))
    ])
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_td64series_pow_invalid(self, scalar_td, box):
        td1 = pd.Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        td1 = tm.box_expected(td1, box)

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot|not supported'
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td ** td1

        with tm.assert_raises_regex(TypeError, pattern):
            td1 ** scalar_td
