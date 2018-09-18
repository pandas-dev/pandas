# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

from datetime import datetime, timedelta
import operator

from numpy import nan
import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, isna, bdate_range,
                    NaT, date_range, timedelta_range, Categorical)
from pandas.core.indexes.datetimes import Timestamp
import pandas.core.nanops as nanops
from pandas.core import ops

from pandas.compat import range
from pandas import compat
from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal)
import pandas.util.testing as tm

from .common import TestData


class TestSeriesComparisons(object):
    def test_comparisons(self):
        left = np.random.randn(10)
        right = np.random.randn(10)
        left[:3] = np.nan

        result = nanops.nangt(left, right)
        with np.errstate(invalid='ignore'):
            expected = (left > right).astype('O')
        expected[:3] = np.nan

        assert_almost_equal(result, expected)

        s = Series(['a', 'b', 'c'])
        s2 = Series([False, True, False])

        # it works!
        exp = Series([False, False, False])
        assert_series_equal(s == s2, exp)
        assert_series_equal(s2 == s, exp)

    def test_categorical_comparisons(self):
        # GH 8938
        # allow equality comparisons
        a = Series(list('abc'), dtype="category")
        b = Series(list('abc'), dtype="object")
        c = Series(['a', 'b', 'cc'], dtype="object")
        d = Series(list('acb'), dtype="object")
        e = Categorical(list('abc'))
        f = Categorical(list('acb'))

        # vs scalar
        assert not (a == 'a').all()
        assert ((a != 'a') == ~(a == 'a')).all()

        assert not ('a' == a).all()
        assert (a == 'a')[0]
        assert ('a' == a)[0]
        assert not ('a' != a)[0]

        # vs list-like
        assert (a == a).all()
        assert not (a != a).all()

        assert (a == list(a)).all()
        assert (a == b).all()
        assert (b == a).all()
        assert ((~(a == b)) == (a != b)).all()
        assert ((~(b == a)) == (b != a)).all()

        assert not (a == c).all()
        assert not (c == a).all()
        assert not (a == d).all()
        assert not (d == a).all()

        # vs a cat-like
        assert (a == e).all()
        assert (e == a).all()
        assert not (a == f).all()
        assert not (f == a).all()

        assert ((~(a == e) == (a != e)).all())
        assert ((~(e == a) == (e != a)).all())
        assert ((~(a == f) == (a != f)).all())
        assert ((~(f == a) == (f != a)).all())

        # non-equality is not comparable
        with pytest.raises(TypeError):
            a < b
        with pytest.raises(TypeError):
            b < a
        with pytest.raises(TypeError):
            a > b
        with pytest.raises(TypeError):
            b > a

    def test_comparison_tuples(self):
        # GH11339
        # comparisons vs tuple
        s = Series([(1, 1), (1, 2)])

        result = s == (1, 2)
        expected = Series([False, True])
        assert_series_equal(result, expected)

        result = s != (1, 2)
        expected = Series([True, False])
        assert_series_equal(result, expected)

        result = s == (0, 0)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        result = s != (0, 0)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        s = Series([(1, 1), (1, 1)])

        result = s == (1, 1)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        result = s != (1, 1)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        s = Series([frozenset([1]), frozenset([1, 2])])

        result = s == frozenset([1])
        expected = Series([True, False])
        assert_series_equal(result, expected)

    def test_comparison_operators_with_nas(self):
        ser = Series(bdate_range('1/1/2000', periods=10), dtype=object)
        ser[::2] = np.nan

        # test that comparisons work
        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            val = ser[5]

            f = getattr(operator, op)
            result = f(ser, val)

            expected = f(ser.dropna(), val).reindex(ser.index)

            if op == 'ne':
                expected = expected.fillna(True).astype(bool)
            else:
                expected = expected.fillna(False).astype(bool)

            assert_series_equal(result, expected)

            # fffffffuuuuuuuuuuuu
            # result = f(val, s)
            # expected = f(val, s.dropna()).reindex(s.index)
            # assert_series_equal(result, expected)

    @pytest.mark.parametrize('bool_op', [operator.and_,
                                         operator.or_, operator.xor])
    def test_bool_operators_with_nas(self, bool_op):
        # boolean &, |, ^ should work with object arrays and propagate NAs
        ser = Series(bdate_range('1/1/2000', periods=10), dtype=object)
        ser[::2] = np.nan

        mask = ser.isna()
        filled = ser.fillna(ser[0])

        result = bool_op(ser < ser[9], ser > ser[3])

        expected = bool_op(filled < filled[9], filled > filled[3])
        expected[mask] = False
        assert_series_equal(result, expected)

    def test_unequal_categorical_comparison_raises_type_error(self):
        # unequal comparison should raise for unordered cats
        cat = Series(Categorical(list("abc")))
        with pytest.raises(TypeError):
            cat > "b"

        cat = Series(Categorical(list("abc"), ordered=False))
        with pytest.raises(TypeError):
            cat > "b"

        # https://github.com/pandas-dev/pandas/issues/9836#issuecomment-92123057
        # and following comparisons with scalars not in categories should raise
        # for unequal comps, but not for equal/not equal
        cat = Series(Categorical(list("abc"), ordered=True))

        with pytest.raises(TypeError):
            cat < "d"
        with pytest.raises(TypeError):
            cat > "d"
        with pytest.raises(TypeError):
            "d" < cat
        with pytest.raises(TypeError):
            "d" > cat

        tm.assert_series_equal(cat == "d", Series([False, False, False]))
        tm.assert_series_equal(cat != "d", Series([True, True, True]))

    @pytest.mark.parametrize('pair', [
        ([pd.Timestamp('2011-01-01'), NaT, pd.Timestamp('2011-01-03')],
         [NaT, NaT, pd.Timestamp('2011-01-03')]),

        ([pd.Timedelta('1 days'), NaT, pd.Timedelta('3 days')],
         [NaT, NaT, pd.Timedelta('3 days')]),

        ([pd.Period('2011-01', freq='M'), NaT, pd.Period('2011-03', freq='M')],
         [NaT, NaT, pd.Period('2011-03', freq='M')])])
    @pytest.mark.parametrize('reverse', [True, False])
    @pytest.mark.parametrize('box', [Series, Index])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_nat_comparisons(self, dtype, box, reverse, pair):
        l, r = pair
        if reverse:
            # add lhs / rhs switched data
            l, r = r, l

        left = Series(l, dtype=dtype)
        right = box(r, dtype=dtype)
        # Series, Index

        expected = Series([False, False, True])
        assert_series_equal(left == right, expected)

        expected = Series([True, True, False])
        assert_series_equal(left != right, expected)

        expected = Series([False, False, False])
        assert_series_equal(left < right, expected)

        expected = Series([False, False, False])
        assert_series_equal(left > right, expected)

        expected = Series([False, False, True])
        assert_series_equal(left >= right, expected)

        expected = Series([False, False, True])
        assert_series_equal(left <= right, expected)

    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        with pytest.raises(ValueError):
            a < b

        a = Series([1, 2])
        b = Series([2, 3, 4])
        with pytest.raises(ValueError):
            a == b

    def test_comparison_label_based(self):

        # GH 4947
        # comparisons should be label based

        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False], list('abc'))

        expected = Series([False, True, False], list('abc'))
        result = a & b
        assert_series_equal(result, expected)

        expected = Series([True, True, False], list('abc'))
        result = a | b
        assert_series_equal(result, expected)

        expected = Series([True, False, False], list('abc'))
        result = a ^ b
        assert_series_equal(result, expected)

        # rhs is bigger
        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False, True], list('abcd'))

        expected = Series([False, True, False, False], list('abcd'))
        result = a & b
        assert_series_equal(result, expected)

        expected = Series([True, True, False, False], list('abcd'))
        result = a | b
        assert_series_equal(result, expected)

        # filling

        # vs empty
        result = a & Series([])
        expected = Series([False, False, False], list('bca'))
        assert_series_equal(result, expected)

        result = a | Series([])
        expected = Series([True, False, True], list('bca'))
        assert_series_equal(result, expected)

        # vs non-matching
        result = a & Series([1], ['z'])
        expected = Series([False, False, False, False], list('abcz'))
        assert_series_equal(result, expected)

        result = a | Series([1], ['z'])
        expected = Series([True, True, False, False], list('abcz'))
        assert_series_equal(result, expected)

        # identity
        # we would like s[s|e] == s to hold for any e, whether empty or not
        for e in [Series([]), Series([1], ['z']),
                  Series(np.nan, b.index), Series(np.nan, a.index)]:
            result = a[a | e]
            assert_series_equal(result, a[a])

        for e in [Series(['z'])]:
            if compat.PY3:
                with tm.assert_produces_warning(RuntimeWarning):
                    result = a[a | e]
            else:
                result = a[a | e]
            assert_series_equal(result, a[a])

        # vs scalars
        index = list('bca')
        t = Series([True, False, True])

        for v in [True, 1, 2]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, True, True], index=index)
            assert_series_equal(result, expected)

        for v in [np.nan, 'foo']:
            with pytest.raises(TypeError):
                t | v

        for v in [False, 0]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, False, True], index=index)
            assert_series_equal(result, expected)

        for v in [True, 1]:
            result = Series([True, False, True], index=index) & v
            expected = Series([True, False, True], index=index)
            assert_series_equal(result, expected)

        for v in [False, 0]:
            result = Series([True, False, True], index=index) & v
            expected = Series([False, False, False], index=index)
            assert_series_equal(result, expected)
        for v in [np.nan]:
            with pytest.raises(TypeError):
                t & v

    def test_comparison_flex_basic(self):
        left = pd.Series(np.random.randn(10))
        right = pd.Series(np.random.randn(10))

        assert_series_equal(left.eq(right), left == right)
        assert_series_equal(left.ne(right), left != right)
        assert_series_equal(left.le(right), left < right)
        assert_series_equal(left.lt(right), left <= right)
        assert_series_equal(left.gt(right), left > right)
        assert_series_equal(left.ge(right), left >= right)

        # axis
        for axis in [0, None, 'index']:
            assert_series_equal(left.eq(right, axis=axis), left == right)
            assert_series_equal(left.ne(right, axis=axis), left != right)
            assert_series_equal(left.le(right, axis=axis), left < right)
            assert_series_equal(left.lt(right, axis=axis), left <= right)
            assert_series_equal(left.gt(right, axis=axis), left > right)
            assert_series_equal(left.ge(right, axis=axis), left >= right)

        #
        msg = 'No axis named 1 for object type'
        for op in ['eq', 'ne', 'le', 'le', 'gt', 'ge']:
            with tm.assert_raises_regex(ValueError, msg):
                getattr(left, op)(right, axis=1)

    def test_comparison_flex_alignment(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        assert_series_equal(left.eq(right), exp)

        exp = pd.Series([True, True, False, True], index=list('abcd'))
        assert_series_equal(left.ne(right), exp)

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        assert_series_equal(left.le(right), exp)

        exp = pd.Series([False, False, False, False], index=list('abcd'))
        assert_series_equal(left.lt(right), exp)

        exp = pd.Series([False, True, True, False], index=list('abcd'))
        assert_series_equal(left.ge(right), exp)

        exp = pd.Series([False, True, False, False], index=list('abcd'))
        assert_series_equal(left.gt(right), exp)

    def test_comparison_flex_alignment_fill(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        assert_series_equal(left.eq(right, fill_value=2), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        assert_series_equal(left.ne(right, fill_value=2), exp)

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        assert_series_equal(left.le(right, fill_value=0), exp)

        exp = pd.Series([False, False, False, True], index=list('abcd'))
        assert_series_equal(left.lt(right, fill_value=0), exp)

        exp = pd.Series([True, True, True, False], index=list('abcd'))
        assert_series_equal(left.ge(right, fill_value=0), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        assert_series_equal(left.gt(right, fill_value=0), exp)

    def test_logical_ops_with_index(self):
        # GH22092
        ser = Series([True, True, False, False])
        idx1 = Index([True, False, True, False])
        idx2 = Index([1, 0, 1, 0])

        expected = Series([True, False, False, False])
        result1 = ser & idx1
        assert_series_equal(result1, expected)
        result2 = ser & idx2
        assert_series_equal(result2, expected)

        expected = Series([True, True, True, False])
        result1 = ser | idx1
        assert_series_equal(result1, expected)
        result2 = ser | idx2
        assert_series_equal(result2, expected)

        expected = Series([False, True, True, False])
        result1 = ser ^ idx1
        assert_series_equal(result1, expected)
        result2 = ser ^ idx2
        assert_series_equal(result2, expected)

    def test_ne(self):
        ts = Series([3, 4, 5, 6, 7], [3, 4, 5, 6, 7], dtype=float)
        expected = [True, True, False, True, True]
        assert tm.equalContents(ts.index != 5, expected)
        assert tm.equalContents(~(ts.index == 5), expected)

    def test_comp_ops_df_compat(self):
        # GH 1134
        s1 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s2 = pd.Series([2, 2, 2], index=list('ABD'), name='x')

        s3 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s4 = pd.Series([2, 2, 2, 2], index=list('ABCD'), name='x')

        for left, right in [(s1, s2), (s2, s1), (s3, s4), (s4, s3)]:

            msg = "Can only compare identically-labeled Series objects"
            with tm.assert_raises_regex(ValueError, msg):
                left == right

            with tm.assert_raises_regex(ValueError, msg):
                left != right

            with tm.assert_raises_regex(ValueError, msg):
                left < right

            msg = "Can only compare identically-labeled DataFrame objects"
            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() == right.to_frame()

            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() != right.to_frame()

            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() < right.to_frame()


class TestDatetimeSeriesArithmetic(object):

    def test_operators_datetimelike_invalid(self, all_arithmetic_operators):
        # these are all TypeEror ops
        op_str = all_arithmetic_operators

        def check(get_ser, test_ser):

            # check that we are getting a TypeError
            # with 'operate' (from core/ops.py) for the ops that are not
            # defined
            op = getattr(get_ser, op_str, None)
            with tm.assert_raises_regex(TypeError, 'operate|cannot'):
                op(test_ser)

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        if op_str not in ['__sub__', '__rsub__']:
            check(dt1, dt2)

        # ## datetime64 with timetimedelta ###
        # TODO(jreback) __rsub__ should raise?
        if op_str not in ['__add__', '__radd__', '__sub__']:
            check(dt1, td1)

        # 8260, 10763
        # datetime64 with tz
        tz = 'US/Eastern'
        dt1 = Series(date_range('2000-01-01 09:00:00', periods=5,
                                tz=tz), name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan
        td1 = Series(timedelta_range('1 days 1 min', periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan

        if op_str not in ['__add__', '__radd__', '__sub__', '__rsub__']:
            check(dt2, td2)

    def test_operators_datetimelike(self):

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        dt1 - dt2
        dt2 - dt1

        # ## datetime64 with timetimedelta ###
        dt1 + td1
        td1 + dt1
        dt1 - td1
        # TODO: Decide if this ought to work.
        # td1 - dt1

        # ## timetimedelta with datetime64 ###
        td1 + dt1
        dt1 + td1


class TestSeriesOperators(TestData):
    @pytest.mark.parametrize(
        'ts',
        [
            (lambda x: x, lambda x: x * 2, False),
            (lambda x: x, lambda x: x[::2], False),
            (lambda x: x, lambda x: 5, True),
            (lambda x: tm.makeFloatSeries(),
             lambda x: tm.makeFloatSeries(),
             True)
        ])
    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul', 'floordiv',
                                        'truediv', 'div', 'pow'])
    def test_op_method(self, opname, ts):
        # check that Series.{opname} behaves like Series.__{opname}__,
        tser = tm.makeTimeSeries().rename('ts')

        series = ts[0](tser)
        other = ts[1](tser)
        check_reverse = ts[2]

        if opname == 'div' and compat.PY3:
            pytest.skip('div test only for Py3')

        op = getattr(Series, opname)

        if op == 'div':
            alt = operator.truediv
        else:
            alt = getattr(operator, opname)

        result = op(series, other)
        expected = alt(series, other)
        assert_almost_equal(result, expected)
        if check_reverse:
            rop = getattr(Series, "r" + opname)
            result = rop(series, other)
            expected = alt(other, series)
            assert_almost_equal(result, expected)

    def test_neg(self):
        assert_series_equal(-self.series, -1 * self.series)

    def test_invert(self):
        assert_series_equal(-(self.series < 0), ~(self.series < 0))

    def test_operators_empty_int_corner(self):
        s1 = Series([], [], dtype=np.int32)
        s2 = Series({'x': 0.})
        assert_series_equal(s1 * s2, Series([np.nan], index=['x']))

    def test_ops_datetimelike_align(self):
        # GH 7500
        # datetimelike ops need to align
        dt = Series(date_range('2012-1-1', periods=3, freq='D'))
        dt.iloc[2] = np.nan
        dt2 = dt[::-1]

        expected = Series([timedelta(0), timedelta(0), pd.NaT])
        # name is reset
        result = dt2 - dt
        assert_series_equal(result, expected)

        expected = Series(expected, name=0)
        result = (dt2.to_frame() - dt.to_frame())[0]
        assert_series_equal(result, expected)

    @pytest.mark.parametrize('op', [
        operator.and_,
        operator.or_,
        pytest.param(ops.rand_,
                     marks=pytest.mark.xfail(reason="GH#22092 Index "
                                                    "implementation returns "
                                                    "Index",
                                             raises=AssertionError,
                                             strict=True)),
        pytest.param(ops.ror_,
                     marks=pytest.mark.xfail(reason="GH#22092 Index "
                                                    "implementation raises",
                                             raises=ValueError, strict=True))
    ])
    def test_bool_ops_with_index(self, op):
        # GH#19792
        ser = Series([True, False, True])
        idx = pd.Index([False, True, True])

        expected = Series([op(ser[n], idx[n]) for n in range(len(ser))])

        result = op(ser, idx)
        assert_series_equal(result, expected)

    def test_operators_bitwise(self):
        # GH 9016: support bitwise op for integer types
        index = list('bca')

        s_tft = Series([True, False, True], index=index)
        s_fff = Series([False, False, False], index=index)
        s_tff = Series([True, False, False], index=index)
        s_empty = Series([])

        # TODO: unused
        # s_0101 = Series([0, 1, 0, 1])

        s_0123 = Series(range(4), dtype='int64')
        s_3333 = Series([3] * 4)
        s_4444 = Series([4] * 4)

        res = s_tft & s_empty
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_tft | s_empty
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & s_3333
        expected = Series(range(4), dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123 | s_4444
        expected = Series(range(4, 8), dtype='int64')
        assert_series_equal(res, expected)

        s_a0b1c0 = Series([1], list('b'))

        res = s_tft & s_a0b1c0
        expected = s_tff.reindex(list('abc'))
        assert_series_equal(res, expected)

        res = s_tft | s_a0b1c0
        expected = s_tft.reindex(list('abc'))
        assert_series_equal(res, expected)

        n0 = 0
        res = s_tft & n0
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_0123 & n0
        expected = Series([0] * 4)
        assert_series_equal(res, expected)

        n1 = 1
        res = s_tft & n1
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & n1
        expected = Series([0, 1, 0, 1])
        assert_series_equal(res, expected)

        s_1111 = Series([1] * 4, dtype='int8')
        res = s_0123 & s_1111
        expected = Series([0, 1, 0, 1], dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123.astype(np.int16) | s_1111.astype(np.int32)
        expected = Series([1, 1, 3, 3], dtype='int32')
        assert_series_equal(res, expected)

        with pytest.raises(TypeError):
            s_1111 & 'a'
        with pytest.raises(TypeError):
            s_1111 & ['a', 'b', 'c', 'd']
        with pytest.raises(TypeError):
            s_0123 & np.NaN
        with pytest.raises(TypeError):
            s_0123 & 3.14
        with pytest.raises(TypeError):
            s_0123 & [0.1, 4, 3.14, 2]

        # s_0123 will be all false now because of reindexing like s_tft
        if compat.PY3:
            # unable to sort incompatible object via .union.
            exp = Series([False] * 7, index=['b', 'c', 'a', 0, 1, 2, 3])
            with tm.assert_produces_warning(RuntimeWarning):
                assert_series_equal(s_tft & s_0123, exp)
        else:
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'a', 'b', 'c'])
            assert_series_equal(s_tft & s_0123, exp)

        # s_tft will be all false now because of reindexing like s_0123
        if compat.PY3:
            # unable to sort incompatible object via .union.
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'b', 'c', 'a'])
            with tm.assert_produces_warning(RuntimeWarning):
                assert_series_equal(s_0123 & s_tft, exp)
        else:
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'a', 'b', 'c'])
            assert_series_equal(s_0123 & s_tft, exp)

        assert_series_equal(s_0123 & False, Series([False] * 4))
        assert_series_equal(s_0123 ^ False, Series([False, True, True, True]))
        assert_series_equal(s_0123 & [False], Series([False] * 4))
        assert_series_equal(s_0123 & (False), Series([False] * 4))
        assert_series_equal(s_0123 & Series([False, np.NaN, False, False]),
                            Series([False] * 4))

        s_ftft = Series([False, True, False, True])
        assert_series_equal(s_0123 & Series([0.1, 4, -3.14, 2]), s_ftft)

        s_abNd = Series(['a', 'b', np.NaN, 'd'])
        res = s_0123 & s_abNd
        expected = s_ftft
        assert_series_equal(res, expected)

    def test_scalar_na_cmp_corners(self):
        s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])

        def tester(a, b):
            return a & b

        with pytest.raises(TypeError):
            s & datetime(2005, 1, 1)

        s = Series([2, 3, 4, 5, 6, 7, 8, 9, datetime(2005, 1, 1)])
        s[::2] = np.nan

        expected = Series(True, index=s.index)
        expected[::2] = False
        result = s & list(s)
        assert_series_equal(result, expected)

        d = DataFrame({'A': s})
        # TODO: Fix this exception - needs to be fixed! (see GH5035)
        # (previously this was a TypeError because series returned
        # NotImplemented

        # this is an alignment issue; these are equivalent
        # https://github.com/pandas-dev/pandas/issues/5284

        pytest.raises(ValueError, lambda: d.__and__(s, axis='columns'))
        pytest.raises(ValueError, tester, s, d)

        # this is wrong as its not a boolean result
        # result = d.__and__(s,axis='index')

    def test_operators_corner(self):
        series = self.ts

        empty = Series([], index=Index([]))

        result = series + empty
        assert np.isnan(result).all()

        result = empty + Series([], index=Index([]))
        assert len(result) == 0

        # TODO: this returned NotImplemented earlier, what to do?
        # deltas = Series([timedelta(1)] * 5, index=np.arange(5))
        # sub_deltas = deltas[::2]
        # deltas5 = deltas * 5
        # deltas = deltas + sub_deltas

        # float + int
        int_ts = self.ts.astype(int)[:-5]
        added = self.ts + int_ts
        expected = Series(self.ts.values[:-5] + int_ts.values,
                          index=self.ts.index[:-5], name='ts')
        tm.assert_series_equal(added[:-5], expected)

    pairings = []
    for op in ['add', 'sub', 'mul', 'pow', 'truediv', 'floordiv']:
        fv = 0
        lop = getattr(Series, op)
        lequiv = getattr(operator, op)
        rop = getattr(Series, 'r' + op)
        # bind op at definition time...
        requiv = lambda x, y, op=op: getattr(operator, op)(y, x)
        pairings.append((lop, lequiv, fv))
        pairings.append((rop, requiv, fv))
    if compat.PY3:
        pairings.append((Series.div, operator.truediv, 1))
        pairings.append((Series.rdiv, lambda x, y: operator.truediv(y, x), 1))
    else:
        pairings.append((Series.div, operator.div, 1))
        pairings.append((Series.rdiv, lambda x, y: operator.div(y, x), 1))

    @pytest.mark.parametrize('op, equiv_op, fv', pairings)
    def test_operators_combine(self, op, equiv_op, fv):
        def _check_fill(meth, op, a, b, fill_value=0):
            exp_index = a.index.union(b.index)
            a = a.reindex(exp_index)
            b = b.reindex(exp_index)

            amask = isna(a)
            bmask = isna(b)

            exp_values = []
            for i in range(len(exp_index)):
                with np.errstate(all='ignore'):
                    if amask[i]:
                        if bmask[i]:
                            exp_values.append(nan)
                            continue
                        exp_values.append(op(fill_value, b[i]))
                    elif bmask[i]:
                        if amask[i]:
                            exp_values.append(nan)
                            continue
                        exp_values.append(op(a[i], fill_value))
                    else:
                        exp_values.append(op(a[i], b[i]))

            result = meth(a, b, fill_value=fill_value)
            expected = Series(exp_values, exp_index)
            assert_series_equal(result, expected)

        a = Series([nan, 1., 2., 3., nan], index=np.arange(5))
        b = Series([nan, 1, nan, 3, nan, 4.], index=np.arange(6))

        result = op(a, b)
        exp = equiv_op(a, b)
        assert_series_equal(result, exp)
        _check_fill(op, equiv_op, a, b, fill_value=fv)
        # should accept axis=0 or axis='rows'
        op(a, b, axis=0)

    def test_operators_na_handling(self):
        from decimal import Decimal
        from datetime import date
        s = Series([Decimal('1.3'), Decimal('2.3')],
                   index=[date(2012, 1, 1), date(2012, 1, 2)])

        result = s + s.shift(1)
        result2 = s.shift(1) + s
        assert isna(result[0])
        assert isna(result2[0])

    def test_op_duplicate_index(self):
        # GH14227
        s1 = Series([1, 2], index=[1, 1])
        s2 = Series([10, 10], index=[1, 2])
        result = s1 + s2
        expected = pd.Series([11, 12, np.nan], index=[1, 1, 2])
        assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "test_input,error_type",
        [
            (pd.Series([]), ValueError),

            # For strings, or any Series with dtype 'O'
            (pd.Series(['foo', 'bar', 'baz']), TypeError),
            (pd.Series([(1,), (2,)]), TypeError),

            # For mixed data types
            (
                pd.Series(['foo', 'foo', 'bar', 'bar', None, np.nan, 'baz']),
                TypeError
            ),
        ]
    )
    def test_assert_idxminmax_raises(self, test_input, error_type):
        """
        Cases where ``Series.argmax`` and related should raise an exception
        """
        with pytest.raises(error_type):
            test_input.idxmin()
        with pytest.raises(error_type):
            test_input.idxmin(skipna=False)
        with pytest.raises(error_type):
            test_input.idxmax()
        with pytest.raises(error_type):
            test_input.idxmax(skipna=False)

    def test_idxminmax_with_inf(self):
        # For numeric data with NA and Inf (GH #13595)
        s = pd.Series([0, -np.inf, np.inf, np.nan])

        assert s.idxmin() == 1
        assert np.isnan(s.idxmin(skipna=False))

        assert s.idxmax() == 2
        assert np.isnan(s.idxmax(skipna=False))

        # Using old-style behavior that treats floating point nan, -inf, and
        # +inf as missing
        with pd.option_context('mode.use_inf_as_na', True):
            assert s.idxmin() == 0
            assert np.isnan(s.idxmin(skipna=False))
            assert s.idxmax() == 0
            np.isnan(s.idxmax(skipna=False))


class TestSeriesOperationsDataFrameCompat(object):

    def test_bool_ops_df_compat(self):
        # GH 1134
        s1 = pd.Series([True, False, True], index=list('ABC'), name='x')
        s2 = pd.Series([True, True, False], index=list('ABD'), name='x')

        exp = pd.Series([True, False, False, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s1 & s2, exp)
        assert_series_equal(s2 & s1, exp)

        # True | np.nan => True
        exp = pd.Series([True, True, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s1 | s2, exp)
        # np.nan | True => np.nan, filled with False
        exp = pd.Series([True, True, False, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s2 | s1, exp)

        # DataFrame doesn't fill nan with False
        exp = pd.DataFrame({'x': [True, False, np.nan, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s1.to_frame() & s2.to_frame(), exp)
        assert_frame_equal(s2.to_frame() & s1.to_frame(), exp)

        exp = pd.DataFrame({'x': [True, True, np.nan, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s1.to_frame() | s2.to_frame(), exp)
        assert_frame_equal(s2.to_frame() | s1.to_frame(), exp)

        # different length
        s3 = pd.Series([True, False, True], index=list('ABC'), name='x')
        s4 = pd.Series([True, True, True, True], index=list('ABCD'), name='x')

        exp = pd.Series([True, False, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s3 & s4, exp)
        assert_series_equal(s4 & s3, exp)

        # np.nan | True => np.nan, filled with False
        exp = pd.Series([True, True, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s3 | s4, exp)
        # True | np.nan => True
        exp = pd.Series([True, True, True, True],
                        index=list('ABCD'), name='x')
        assert_series_equal(s4 | s3, exp)

        exp = pd.DataFrame({'x': [True, False, True, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s3.to_frame() & s4.to_frame(), exp)
        assert_frame_equal(s4.to_frame() & s3.to_frame(), exp)

        exp = pd.DataFrame({'x': [True, True, True, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s3.to_frame() | s4.to_frame(), exp)
        assert_frame_equal(s4.to_frame() | s3.to_frame(), exp)
