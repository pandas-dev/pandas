# coding=utf-8
from datetime import datetime
import operator

import pytest

import numpy as np
import pandas as pd

from pandas import Index, Series, NaT, Categorical, date_range, bdate_range
from pandas.core import nanops
from pandas import compat

import pandas.util.testing as tm


class TestSeriesComparisons(object):
    def test_series_comparison_scalars(self):
        series = Series(date_range('1/1/2000', periods=10))

        val = datetime(2000, 1, 4)
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

        val = series[5]
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

    def test_comparisons(self):
        left = np.random.randn(10)
        right = np.random.randn(10)
        left[:3] = np.nan

        result = nanops.nangt(left, right)
        with np.errstate(invalid='ignore'):
            expected = (left > right).astype('O')
        expected[:3] = np.nan

        tm.assert_almost_equal(result, expected)

        s = Series(['a', 'b', 'c'])
        s2 = Series([False, True, False])

        # it works!
        exp = Series([False, False, False])
        tm.assert_series_equal(s == s2, exp)
        tm.assert_series_equal(s2 == s, exp)

    def test_operator_series_comparison_zerorank(self):
        # GH 13006
        result = np.float64(0) > pd.Series([1, 2, 3])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)
        result = pd.Series([1, 2, 3]) < np.float64(0)
        expected = pd.Series([1, 2, 3]) < 0.0
        tm.assert_series_equal(result, expected)
        result = np.array([0, 1, 2])[0] > pd.Series([0, 1, 2])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)

    def test_object_comparisons(self):
        s = Series(['a', 'b', np.nan, 'c', 'a'])

        result = s == 'a'
        expected = Series([True, False, False, False, True])
        tm.assert_series_equal(result, expected)

        result = s < 'a'
        expected = Series([False, False, False, False, False])
        tm.assert_series_equal(result, expected)

        result = s != 'a'
        expected = -(s == 'a')
        tm.assert_series_equal(result, expected)

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
        pytest.raises(TypeError, lambda: a < b)
        pytest.raises(TypeError, lambda: b < a)
        pytest.raises(TypeError, lambda: a > b)
        pytest.raises(TypeError, lambda: b > a)

    def test_comparison_tuples(self):
        # GH11339
        # comparisons vs tuple
        s = Series([(1, 1), (1, 2)])

        result = s == (1, 2)
        expected = Series([False, True])
        tm.assert_series_equal(result, expected)

        result = s != (1, 2)
        expected = Series([True, False])
        tm.assert_series_equal(result, expected)

        result = s == (0, 0)
        expected = Series([False, False])
        tm.assert_series_equal(result, expected)

        result = s != (0, 0)
        expected = Series([True, True])
        tm.assert_series_equal(result, expected)

        s = Series([(1, 1), (1, 1)])

        result = s == (1, 1)
        expected = Series([True, True])
        tm.assert_series_equal(result, expected)

        result = s != (1, 1)
        expected = Series([False, False])
        tm.assert_series_equal(result, expected)

        s = Series([frozenset([1]), frozenset([1, 2])])

        result = s == frozenset([1])
        expected = Series([True, False])
        tm.assert_series_equal(result, expected)

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

            tm.assert_series_equal(result, expected)

            # fffffffuuuuuuuuuuuu
            # result = f(val, s)
            # expected = f(val, s.dropna()).reindex(s.index)
            # tm.assert_series_equal(result, expected)

            # boolean &, |, ^ should work with object arrays and propagate NAs

        ops = ['and_', 'or_', 'xor']
        mask = ser.isna()
        for bool_op in ops:
            func = getattr(operator, bool_op)

            filled = ser.fillna(ser[0])

            result = func(ser < ser[9], ser > ser[3])

            expected = func(filled < filled[9], filled > filled[3])
            expected[mask] = False
            tm.assert_series_equal(result, expected)

    def test_comparison_object_numeric_nas(self):
        ser = Series(np.random.randn(10), dtype=object)
        shifted = ser.shift(2)

        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            func = getattr(operator, op)

            result = func(ser, shifted)
            expected = func(ser.astype(float), shifted.astype(float))
            tm.assert_series_equal(result, expected)

    def test_comparison_invalid(self):
        # GH4968
        # invalid date/int comparisons
        s = Series(range(5))
        s2 = Series(date_range('20010101', periods=5))

        for (x, y) in [(s, s2), (s2, s)]:
            pytest.raises(TypeError, lambda: x == y)
            pytest.raises(TypeError, lambda: x != y)
            pytest.raises(TypeError, lambda: x >= y)
            pytest.raises(TypeError, lambda: x > y)
            pytest.raises(TypeError, lambda: x < y)
            pytest.raises(TypeError, lambda: x <= y)

    def test_unequal_categorical_comparison_raises_type_error(self):
        # unequal comparison should raise for unordered cats
        cat = Series(Categorical(list("abc")))

        def f():
            cat > "b"

        pytest.raises(TypeError, f)
        cat = Series(Categorical(list("abc"), ordered=False))

        def f():
            cat > "b"

        pytest.raises(TypeError, f)

        # https://github.com/pandas-dev/pandas/issues/9836#issuecomment-92123057
        # and following comparisons with scalars not in categories should raise
        # for unequal comps, but not for equal/not equal
        cat = Series(Categorical(list("abc"), ordered=True))

        pytest.raises(TypeError, lambda: cat < "d")
        pytest.raises(TypeError, lambda: cat > "d")
        pytest.raises(TypeError, lambda: "d" < cat)
        pytest.raises(TypeError, lambda: "d" > cat)

        tm.assert_series_equal(cat == "d", Series([False, False, False]))
        tm.assert_series_equal(cat != "d", Series([True, True, True]))

    @pytest.mark.parametrize('dtype', [None, object])
    def test_more_na_comparisons(self, dtype):
        left = Series(['a', np.nan, 'c'], dtype=dtype)
        right = Series(['a', np.nan, 'd'], dtype=dtype)

        result = left == right
        expected = Series([True, False, False])
        tm.assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        tm.assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        tm.assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        tm.assert_series_equal(result, expected)

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
        tm.assert_series_equal(left == right, expected)

        expected = Series([True, True, False])
        tm.assert_series_equal(left != right, expected)

        expected = Series([False, False, False])
        tm.assert_series_equal(left < right, expected)

        expected = Series([False, False, False])
        tm.assert_series_equal(left > right, expected)

        expected = Series([False, False, True])
        tm.assert_series_equal(left >= right, expected)

        expected = Series([False, False, True])
        tm.assert_series_equal(left <= right, expected)

    @pytest.mark.parametrize('data', [
        [pd.Timestamp('2011-01-01'), NaT, pd.Timestamp('2011-01-03')],
        [pd.Timedelta('1 days'), NaT, pd.Timedelta('3 days')],
        [pd.Period('2011-01', freq='M'), NaT, pd.Period('2011-03', freq='M')]
    ])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_nat_comparisons_scalar(self, dtype, data):
        left = Series(data, dtype=dtype)

        expected = Series([False, False, False])
        tm.assert_series_equal(left == pd.NaT, expected)
        tm.assert_series_equal(pd.NaT == left, expected)

        expected = Series([True, True, True])
        tm.assert_series_equal(left != pd.NaT, expected)
        tm.assert_series_equal(pd.NaT != left, expected)

        expected = Series([False, False, False])
        tm.assert_series_equal(left < pd.NaT, expected)
        tm.assert_series_equal(pd.NaT > left, expected)
        tm.assert_series_equal(left <= pd.NaT, expected)
        tm.assert_series_equal(pd.NaT >= left, expected)

        tm.assert_series_equal(left > pd.NaT, expected)
        tm.assert_series_equal(pd.NaT < left, expected)
        tm.assert_series_equal(left >= pd.NaT, expected)
        tm.assert_series_equal(pd.NaT <= left, expected)

    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        pytest.raises(ValueError, a.__lt__, b)

        a = Series([1, 2])
        b = Series([2, 3, 4])
        pytest.raises(ValueError, a.__eq__, b)

    def test_comparison_label_based(self):

        # GH 4947
        # comparisons should be label based

        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False], list('abc'))

        expected = Series([False, True, False], list('abc'))
        result = a & b
        tm.assert_series_equal(result, expected)

        expected = Series([True, True, False], list('abc'))
        result = a | b
        tm.assert_series_equal(result, expected)

        expected = Series([True, False, False], list('abc'))
        result = a ^ b
        tm.assert_series_equal(result, expected)

        # rhs is bigger
        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False, True], list('abcd'))

        expected = Series([False, True, False, False], list('abcd'))
        result = a & b
        tm.assert_series_equal(result, expected)

        expected = Series([True, True, False, False], list('abcd'))
        result = a | b
        tm.assert_series_equal(result, expected)

        # filling

        # vs empty
        result = a & Series([])
        expected = Series([False, False, False], list('bca'))
        tm.assert_series_equal(result, expected)

        result = a | Series([])
        expected = Series([True, False, True], list('bca'))
        tm.assert_series_equal(result, expected)

        # vs non-matching
        result = a & Series([1], ['z'])
        expected = Series([False, False, False, False], list('abcz'))
        tm.assert_series_equal(result, expected)

        result = a | Series([1], ['z'])
        expected = Series([True, True, False, False], list('abcz'))
        tm.assert_series_equal(result, expected)

        # identity
        # we would like s[s|e] == s to hold for any e, whether empty or not
        for e in [Series([]), Series([1], ['z']),
                  Series(np.nan, b.index), Series(np.nan, a.index)]:
            result = a[a | e]
            tm.assert_series_equal(result, a[a])

        for e in [Series(['z'])]:
            if compat.PY3:
                with tm.assert_produces_warning(RuntimeWarning):
                    result = a[a | e]
            else:
                result = a[a | e]
            tm.assert_series_equal(result, a[a])

        # vs scalars
        index = list('bca')
        t = Series([True, False, True])

        for v in [True, 1, 2]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, True, True], index=index)
            tm.assert_series_equal(result, expected)

        for v in [np.nan, 'foo']:
            pytest.raises(TypeError, lambda: t | v)

        for v in [False, 0]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, False, True], index=index)
            tm.assert_series_equal(result, expected)

        for v in [True, 1]:
            result = Series([True, False, True], index=index) & v
            expected = Series([True, False, True], index=index)
            tm.assert_series_equal(result, expected)

        for v in [False, 0]:
            result = Series([True, False, True], index=index) & v
            expected = Series([False, False, False], index=index)
            tm.assert_series_equal(result, expected)
        for v in [np.nan]:
            pytest.raises(TypeError, lambda: t & v)

    def test_comparison_flex_basic(self):
        left = pd.Series(np.random.randn(10))
        right = pd.Series(np.random.randn(10))

        tm.assert_series_equal(left.eq(right), left == right)
        tm.assert_series_equal(left.ne(right), left != right)
        tm.assert_series_equal(left.le(right), left < right)
        tm.assert_series_equal(left.lt(right), left <= right)
        tm.assert_series_equal(left.gt(right), left > right)
        tm.assert_series_equal(left.ge(right), left >= right)

        # axis
        for axis in [0, None, 'index']:
            tm.assert_series_equal(left.eq(right, axis=axis), left == right)
            tm.assert_series_equal(left.ne(right, axis=axis), left != right)
            tm.assert_series_equal(left.le(right, axis=axis), left < right)
            tm.assert_series_equal(left.lt(right, axis=axis), left <= right)
            tm.assert_series_equal(left.gt(right, axis=axis), left > right)
            tm.assert_series_equal(left.ge(right, axis=axis), left >= right)

        #
        msg = 'No axis named 1 for object type'
        for op in ['eq', 'ne', 'le', 'le', 'gt', 'ge']:
            with tm.assert_raises_regex(ValueError, msg):
                getattr(left, op)(right, axis=1)

    def test_comparison_flex_alignment(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        tm.assert_series_equal(left.eq(right), exp)

        exp = pd.Series([True, True, False, True], index=list('abcd'))
        tm.assert_series_equal(left.ne(right), exp)

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        tm.assert_series_equal(left.le(right), exp)

        exp = pd.Series([False, False, False, False], index=list('abcd'))
        tm.assert_series_equal(left.lt(right), exp)

        exp = pd.Series([False, True, True, False], index=list('abcd'))
        tm.assert_series_equal(left.ge(right), exp)

        exp = pd.Series([False, True, False, False], index=list('abcd'))
        tm.assert_series_equal(left.gt(right), exp)

    def test_comparison_flex_alignment_fill(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        tm.assert_series_equal(left.eq(right, fill_value=2), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        tm.assert_series_equal(left.ne(right, fill_value=2), exp)

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        tm.assert_series_equal(left.le(right, fill_value=0), exp)

        exp = pd.Series([False, False, False, True], index=list('abcd'))
        tm.assert_series_equal(left.lt(right, fill_value=0), exp)

        exp = pd.Series([True, True, True, False], index=list('abcd'))
        tm.assert_series_equal(left.ge(right, fill_value=0), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        tm.assert_series_equal(left.gt(right, fill_value=0), exp)

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
