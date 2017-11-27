# -*- coding: utf-8 -*-
import pytest
from datetime import datetime

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, DataFrame, Timestamp,
                    CategoricalIndex, date_range, DatetimeIndex,
                    period_range, timedelta_range, NaT,
                    Interval, IntervalIndex)
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.dtypes.common import is_float_dtype, is_integer_dtype
from pandas.tests.categorical.common import TestCategoricalBlock


class TestCategoricalConstructors(object):

    def test_validate_ordered(self):
        # see gh-14058
        exp_msg = "'ordered' must either be 'True' or 'False'"
        exp_err = TypeError

        # This should be a boolean.
        ordered = np.array([0, 1, 2])

        with tm.assert_raises_regex(exp_err, exp_msg):
            Categorical([1, 2, 3], ordered=ordered)

        with tm.assert_raises_regex(exp_err, exp_msg):
            Categorical.from_codes([0, 0, 1], categories=['a', 'b', 'c'],
                                   ordered=ordered)

    def test_constructor_empty(self):
        # GH 17248
        c = Categorical([])
        expected = Index([])
        tm.assert_index_equal(c.categories, expected)

        c = Categorical([], categories=[1, 2, 3])
        expected = pd.Int64Index([1, 2, 3])
        tm.assert_index_equal(c.categories, expected)

    def test_constructor_tuples(self):
        values = np.array([(1,), (1, 2), (1,), (1, 2)], dtype=object)
        result = Categorical(values)
        expected = Index([(1,), (1, 2)], tupleize_cols=False)
        tm.assert_index_equal(result.categories, expected)
        assert result.ordered is False

    def test_constructor_tuples_datetimes(self):
        # numpy will auto reshape when all of the tuples are the
        # same len, so add an extra one with 2 items and slice it off
        values = np.array([(Timestamp('2010-01-01'),),
                           (Timestamp('2010-01-02'),),
                           (Timestamp('2010-01-01'),),
                           (Timestamp('2010-01-02'),),
                           ('a', 'b')], dtype=object)[:-1]
        result = Categorical(values)
        expected = Index([(Timestamp('2010-01-01'),),
                          (Timestamp('2010-01-02'),)], tupleize_cols=False)
        tm.assert_index_equal(result.categories, expected)

    def test_constructor_unsortable(self):

        # it works!
        arr = np.array([1, 2, 3, datetime.now()], dtype='O')
        factor = Categorical(arr, ordered=False)
        assert not factor.ordered

        # this however will raise as cannot be sorted
        pytest.raises(
            TypeError, lambda: Categorical(arr, ordered=True))

    def test_constructor_interval(self):
        result = Categorical([Interval(1, 2), Interval(2, 3), Interval(3, 6)],
                             ordered=True)
        ii = IntervalIndex.from_intervals([Interval(1, 2),
                                           Interval(2, 3),
                                           Interval(3, 6)])
        exp = Categorical(ii, ordered=True)
        tm.assert_categorical_equal(result, exp)
        tm.assert_index_equal(result.categories, ii)

    def test_constructor(self):

        exp_arr = np.array(["a", "b", "c", "a", "b", "c"], dtype=np.object_)
        c1 = Categorical(exp_arr)
        tm.assert_numpy_array_equal(c1.__array__(), exp_arr)
        c2 = Categorical(exp_arr, categories=["a", "b", "c"])
        tm.assert_numpy_array_equal(c2.__array__(), exp_arr)
        c2 = Categorical(exp_arr, categories=["c", "b", "a"])
        tm.assert_numpy_array_equal(c2.__array__(), exp_arr)

        # categories must be unique
        def f():
            Categorical([1, 2], [1, 2, 2])

        pytest.raises(ValueError, f)

        def f():
            Categorical(["a", "b"], ["a", "b", "b"])

        pytest.raises(ValueError, f)

        # The default should be unordered
        c1 = Categorical(["a", "b", "c", "a"])
        assert not c1.ordered

        # Categorical as input
        c1 = Categorical(["a", "b", "c", "a"])
        c2 = Categorical(c1)
        tm.assert_categorical_equal(c1, c2)

        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        c2 = Categorical(c1)
        tm.assert_categorical_equal(c1, c2)

        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "c", "b"])
        c2 = Categorical(c1)
        tm.assert_categorical_equal(c1, c2)

        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "c", "b"])
        c2 = Categorical(c1, categories=["a", "b", "c"])
        tm.assert_numpy_array_equal(c1.__array__(), c2.__array__())
        tm.assert_index_equal(c2.categories, Index(["a", "b", "c"]))

        # Series of dtype category
        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        c2 = Categorical(Series(c1))
        tm.assert_categorical_equal(c1, c2)

        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "c", "b"])
        c2 = Categorical(Series(c1))
        tm.assert_categorical_equal(c1, c2)

        # Series
        c1 = Categorical(["a", "b", "c", "a"])
        c2 = Categorical(Series(["a", "b", "c", "a"]))
        tm.assert_categorical_equal(c1, c2)

        c1 = Categorical(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        c2 = Categorical(Series(["a", "b", "c", "a"]),
                         categories=["a", "b", "c", "d"])
        tm.assert_categorical_equal(c1, c2)

        # This should result in integer categories, not float!
        cat = Categorical([1, 2, 3, np.nan], categories=[1, 2, 3])
        assert is_integer_dtype(cat.categories)

        # https://github.com/pandas-dev/pandas/issues/3678
        cat = Categorical([np.nan, 1, 2, 3])
        assert is_integer_dtype(cat.categories)

        # this should result in floats
        cat = Categorical([np.nan, 1, 2., 3])
        assert is_float_dtype(cat.categories)

        cat = Categorical([np.nan, 1., 2., 3.])
        assert is_float_dtype(cat.categories)

        # This doesn't work -> this would probably need some kind of "remember
        # the original type" feature to try to cast the array interface result
        # to...

        # vals = np.asarray(cat[cat.notna()])
        # assert is_integer_dtype(vals)

        # corner cases
        cat = Categorical([1])
        assert len(cat.categories) == 1
        assert cat.categories[0] == 1
        assert len(cat.codes) == 1
        assert cat.codes[0] == 0

        cat = Categorical(["a"])
        assert len(cat.categories) == 1
        assert cat.categories[0] == "a"
        assert len(cat.codes) == 1
        assert cat.codes[0] == 0

        # Scalars should be converted to lists
        cat = Categorical(1)
        assert len(cat.categories) == 1
        assert cat.categories[0] == 1
        assert len(cat.codes) == 1
        assert cat.codes[0] == 0

        # two arrays
        #  - when the first is an integer dtype and the second is not
        #  - when the resulting codes are all -1/NaN
        with tm.assert_produces_warning(None):
            c_old = Categorical([0, 1, 2, 0, 1, 2],
                                categories=["a", "b", "c"])  # noqa

        with tm.assert_produces_warning(None):
            c_old = Categorical([0, 1, 2, 0, 1, 2],  # noqa
                                categories=[3, 4, 5])

        # the next one are from the old docs
        with tm.assert_produces_warning(None):
            c_old2 = Categorical([0, 1, 2, 0, 1, 2], [1, 2, 3])  # noqa
            cat = Categorical([1, 2], categories=[1, 2, 3])

        # this is a legitimate constructor
        with tm.assert_produces_warning(None):
            c = Categorical(np.array([], dtype='int64'),  # noqa
                            categories=[3, 2, 1], ordered=True)

    def test_constructor_not_sequence(self):
        # https://github.com/pandas-dev/pandas/issues/16022
        with pytest.raises(TypeError):
            Categorical(['a', 'b'], categories='a')

    def test_constructor_with_null(self):

        # Cannot have NaN in categories
        with pytest.raises(ValueError):
            Categorical([np.nan, "a", "b", "c"],
                        categories=[np.nan, "a", "b", "c"])

        with pytest.raises(ValueError):
            Categorical([None, "a", "b", "c"],
                        categories=[None, "a", "b", "c"])

        with pytest.raises(ValueError):
            Categorical(DatetimeIndex(['nat', '20160101']),
                        categories=[NaT, Timestamp('20160101')])

    def test_constructor_with_index(self):
        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        tm.assert_categorical_equal(ci.values, Categorical(ci))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        tm.assert_categorical_equal(ci.values,
                                    Categorical(ci.astype(object),
                                                categories=ci.categories))

    def test_constructor_with_generator(self):
        # This was raising an Error in isna(single_val).any() because isna
        # returned a scalar for a generator
        xrange = range

        exp = Categorical([0, 1, 2])
        cat = Categorical((x for x in [0, 1, 2]))
        tm.assert_categorical_equal(cat, exp)
        cat = Categorical(xrange(3))
        tm.assert_categorical_equal(cat, exp)

        # This uses xrange internally
        from pandas.core.index import MultiIndex
        MultiIndex.from_product([range(5), ['a', 'b', 'c']])

        # check that categories accept generators and sequences
        cat = Categorical([0, 1, 2], categories=(x for x in [0, 1, 2]))
        tm.assert_categorical_equal(cat, exp)
        cat = Categorical([0, 1, 2], categories=xrange(3))
        tm.assert_categorical_equal(cat, exp)

    def test_constructor_with_datetimelike(self):

        # 12077
        # constructor wwth a datetimelike and NaT

        for dtl in [date_range('1995-01-01 00:00:00', periods=5, freq='s'),
                    date_range('1995-01-01 00:00:00', periods=5,
                               freq='s', tz='US/Eastern'),
                    timedelta_range('1 day', periods=5, freq='s')]:

            s = Series(dtl)
            c = Categorical(s)
            expected = type(dtl)(s)
            expected.freq = None
            tm.assert_index_equal(c.categories, expected)
            tm.assert_numpy_array_equal(c.codes, np.arange(5, dtype='int8'))

            # with NaT
            s2 = s.copy()
            s2.iloc[-1] = NaT
            c = Categorical(s2)
            expected = type(dtl)(s2.dropna())
            expected.freq = None
            tm.assert_index_equal(c.categories, expected)

            exp = np.array([0, 1, 2, 3, -1], dtype=np.int8)
            tm.assert_numpy_array_equal(c.codes, exp)

            result = repr(c)
            assert 'NaT' in result

    def test_constructor_from_index_series_datetimetz(self):
        idx = date_range('2015-01-01 10:00', freq='D', periods=3,
                         tz='US/Eastern')
        result = Categorical(idx)
        tm.assert_index_equal(result.categories, idx)

        result = Categorical(Series(idx))
        tm.assert_index_equal(result.categories, idx)

    def test_constructor_from_index_series_timedelta(self):
        idx = timedelta_range('1 days', freq='D', periods=3)
        result = Categorical(idx)
        tm.assert_index_equal(result.categories, idx)

        result = Categorical(Series(idx))
        tm.assert_index_equal(result.categories, idx)

    def test_constructor_from_index_series_period(self):
        idx = period_range('2015-01-01', freq='D', periods=3)
        result = Categorical(idx)
        tm.assert_index_equal(result.categories, idx)

        result = Categorical(Series(idx))
        tm.assert_index_equal(result.categories, idx)

    def test_constructor_invariant(self):
        # GH 14190
        vals = [
            np.array([1., 1.2, 1.8, np.nan]),
            np.array([1, 2, 3], dtype='int64'),
            ['a', 'b', 'c', np.nan],
            [pd.Period('2014-01'), pd.Period('2014-02'), NaT],
            [Timestamp('2014-01-01'), Timestamp('2014-01-02'), NaT],
            [Timestamp('2014-01-01', tz='US/Eastern'),
             Timestamp('2014-01-02', tz='US/Eastern'), NaT],
        ]
        for val in vals:
            c = Categorical(val)
            c2 = Categorical(c)
            tm.assert_categorical_equal(c, c2)

    @pytest.mark.parametrize('ordered', [True, False])
    def test_constructor_with_dtype(self, ordered):
        categories = ['b', 'a', 'c']
        dtype = CategoricalDtype(categories, ordered=ordered)
        result = Categorical(['a', 'b', 'a', 'c'], dtype=dtype)
        expected = Categorical(['a', 'b', 'a', 'c'], categories=categories,
                               ordered=ordered)
        tm.assert_categorical_equal(result, expected)
        assert result.ordered is ordered

    def test_constructor_dtype_and_others_raises(self):
        dtype = CategoricalDtype(['a', 'b'], ordered=True)
        with tm.assert_raises_regex(ValueError, "Cannot"):
            Categorical(['a', 'b'], categories=['a', 'b'], dtype=dtype)

        with tm.assert_raises_regex(ValueError, "Cannot"):
            Categorical(['a', 'b'], ordered=True, dtype=dtype)

        with tm.assert_raises_regex(ValueError, "Cannot"):
            Categorical(['a', 'b'], ordered=False, dtype=dtype)

    @pytest.mark.parametrize('categories', [
        None, ['a', 'b'], ['a', 'c'],
    ])
    @pytest.mark.parametrize('ordered', [True, False])
    def test_constructor_str_category(self, categories, ordered):
        result = Categorical(['a', 'b'], categories=categories,
                             ordered=ordered, dtype='category')
        expected = Categorical(['a', 'b'], categories=categories,
                               ordered=ordered)
        tm.assert_categorical_equal(result, expected)

    def test_constructor_str_unknown(self):
        with tm.assert_raises_regex(ValueError, "Unknown `dtype`"):
            Categorical([1, 2], dtype="foo")

    def test_constructor_from_categorical_with_dtype(self):
        dtype = CategoricalDtype(['a', 'b', 'c'], ordered=True)
        values = Categorical(['a', 'b', 'd'])
        result = Categorical(values, dtype=dtype)
        # We use dtype.categories, not values.categories
        expected = Categorical(['a', 'b', 'd'], categories=['a', 'b', 'c'],
                               ordered=True)
        tm.assert_categorical_equal(result, expected)

    def test_constructor_from_categorical_with_unknown_dtype(self):
        dtype = CategoricalDtype(None, ordered=True)
        values = Categorical(['a', 'b', 'd'])
        result = Categorical(values, dtype=dtype)
        # We use values.categories, not dtype.categories
        expected = Categorical(['a', 'b', 'd'], categories=['a', 'b', 'd'],
                               ordered=True)
        tm.assert_categorical_equal(result, expected)

    def test_contructor_from_categorical_string(self):
        values = Categorical(['a', 'b', 'd'])
        # use categories, ordered
        result = Categorical(values, categories=['a', 'b', 'c'], ordered=True,
                             dtype='category')
        expected = Categorical(['a', 'b', 'd'], categories=['a', 'b', 'c'],
                               ordered=True)
        tm.assert_categorical_equal(result, expected)

        # No string
        result = Categorical(values, categories=['a', 'b', 'c'], ordered=True)
        tm.assert_categorical_equal(result, expected)

    def test_constructor_with_categorical_categories(self):
        # GH17884
        expected = Categorical(['a', 'b'], categories=['a', 'b', 'c'])

        result = Categorical(
            ['a', 'b'], categories=Categorical(['a', 'b', 'c']))
        tm.assert_categorical_equal(result, expected)

        result = Categorical(
            ['a', 'b'], categories=CategoricalIndex(['a', 'b', 'c']))
        tm.assert_categorical_equal(result, expected)

    def test_from_codes(self):

        # too few categories
        def f():
            Categorical.from_codes([1, 2], [1, 2])

        pytest.raises(ValueError, f)

        # no int codes
        def f():
            Categorical.from_codes(["a"], [1, 2])

        pytest.raises(ValueError, f)

        # no unique categories
        def f():
            Categorical.from_codes([0, 1, 2], ["a", "a", "b"])

        pytest.raises(ValueError, f)

        # NaN categories included
        def f():
            Categorical.from_codes([0, 1, 2], ["a", "b", np.nan])

        pytest.raises(ValueError, f)

        # too negative
        def f():
            Categorical.from_codes([-2, 1, 2], ["a", "b", "c"])

        pytest.raises(ValueError, f)

        exp = Categorical(["a", "b", "c"], ordered=False)
        res = Categorical.from_codes([0, 1, 2], ["a", "b", "c"])
        tm.assert_categorical_equal(exp, res)

        # Not available in earlier numpy versions
        if hasattr(np.random, "choice"):
            codes = np.random.choice([0, 1], 5, p=[0.9, 0.1])
            Categorical.from_codes(codes, categories=["train", "test"])

    def test_from_codes_with_categorical_categories(self):
        # GH17884
        expected = Categorical(['a', 'b'], categories=['a', 'b', 'c'])

        result = Categorical.from_codes(
            [0, 1], categories=Categorical(['a', 'b', 'c']))
        tm.assert_categorical_equal(result, expected)

        result = Categorical.from_codes(
            [0, 1], categories=CategoricalIndex(['a', 'b', 'c']))
        tm.assert_categorical_equal(result, expected)

        # non-unique Categorical still raises
        with pytest.raises(ValueError):
            Categorical.from_codes([0, 1], Categorical(['a', 'b', 'a']))

    @pytest.mark.parametrize('dtype', [None, 'category'])
    def test_from_inferred_categories(self, dtype):
        cats = ['a', 'b']
        codes = np.array([0, 0, 1, 1], dtype='i8')
        result = Categorical._from_inferred_categories(cats, codes, dtype)
        expected = Categorical.from_codes(codes, cats)
        tm.assert_categorical_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, 'category'])
    def test_from_inferred_categories_sorts(self, dtype):
        cats = ['b', 'a']
        codes = np.array([0, 1, 1, 1], dtype='i8')
        result = Categorical._from_inferred_categories(cats, codes, dtype)
        expected = Categorical.from_codes([1, 0, 0, 0], ['a', 'b'])
        tm.assert_categorical_equal(result, expected)

    def test_from_inferred_categories_dtype(self):
        cats = ['a', 'b', 'd']
        codes = np.array([0, 1, 0, 2], dtype='i8')
        dtype = CategoricalDtype(['c', 'b', 'a'], ordered=True)
        result = Categorical._from_inferred_categories(cats, codes, dtype)
        expected = Categorical(['a', 'b', 'a', 'd'],
                               categories=['c', 'b', 'a'],
                               ordered=True)
        tm.assert_categorical_equal(result, expected)

    def test_from_inferred_categories_coerces(self):
        cats = ['1', '2', 'bad']
        codes = np.array([0, 0, 1, 2], dtype='i8')
        dtype = CategoricalDtype([1, 2])
        result = Categorical._from_inferred_categories(cats, codes, dtype)
        expected = Categorical([1, 1, 2, np.nan])
        tm.assert_categorical_equal(result, expected)

    def test_construction_with_ordered(self):
        # GH 9347, 9190
        cat = Categorical([0, 1, 2])
        assert not cat.ordered
        cat = Categorical([0, 1, 2], ordered=False)
        assert not cat.ordered
        cat = Categorical([0, 1, 2], ordered=True)
        assert cat.ordered


class TestCategoricalBlockConstructorsWithFactor(TestCategoricalBlock):

    def test_basic(self):

        # test basic creation / coercion of categoricals
        s = Series(self.factor, name='A')
        assert s.dtype == 'category'
        assert len(s) == len(self.factor)
        str(s.values)
        str(s)

        # in a frame
        df = DataFrame({'A': self.factor})
        result = df['A']
        tm.assert_series_equal(result, s)
        result = df.iloc[:, 0]
        tm.assert_series_equal(result, s)
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        df = DataFrame({'A': s})
        result = df['A']
        tm.assert_series_equal(result, s)
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        # multiples
        df = DataFrame({'A': s, 'B': s, 'C': 1})
        result1 = df['A']
        result2 = df['B']
        tm.assert_series_equal(result1, s)
        tm.assert_series_equal(result2, s, check_names=False)
        assert result2.name == 'B'
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        # GH8623
        x = DataFrame([[1, 'John P. Doe'], [2, 'Jane Dove'],
                       [1, 'John P. Doe']],
                      columns=['person_id', 'person_name'])
        x['person_name'] = Categorical(x.person_name
                                       )  # doing this breaks transform

        expected = x.iloc[0].person_name
        result = x.person_name.iloc[0]
        assert result == expected

        result = x.person_name[0]
        assert result == expected

        result = x.person_name.loc[0]
        assert result == expected


class TestCategoricalBlockConstructors(object):

    def test_construction_series(self):

        l = [1, 2, 3, 1]
        exp = Series(l).astype('category')
        res = Series(l, dtype='category')
        tm.assert_series_equal(res, exp)

        l = ["a", "b", "c", "a"]
        exp = Series(l).astype('category')
        res = Series(l, dtype='category')
        tm.assert_series_equal(res, exp)

        # insert into frame with different index
        # GH 8076
        index = date_range('20000101', periods=3)
        expected = Series(Categorical(values=[np.nan, np.nan, np.nan],
                                      categories=['a', 'b', 'c']))
        expected.index = index

        expected = DataFrame({'x': expected})
        df = DataFrame(
            {'x': Series(['a', 'b', 'c'], dtype='category')}, index=index)
        tm.assert_frame_equal(df, expected)

    def test_construction_frame(self):

        # GH8626

        # dict creation
        df = DataFrame({'A': list('abc')}, dtype='category')
        expected = Series(list('abc'), dtype='category', name='A')
        tm.assert_series_equal(df['A'], expected)

        # to_frame
        s = Series(list('abc'), dtype='category')
        result = s.to_frame()
        expected = Series(list('abc'), dtype='category', name=0)
        tm.assert_series_equal(result[0], expected)
        result = s.to_frame(name='foo')
        expected = Series(list('abc'), dtype='category', name='foo')
        tm.assert_series_equal(result['foo'], expected)

        # list-like creation
        df = DataFrame(list('abc'), dtype='category')
        expected = Series(list('abc'), dtype='category', name=0)
        tm.assert_series_equal(df[0], expected)

        # ndim != 1
        df = DataFrame([Categorical(list('abc'))])
        expected = DataFrame({0: Series(list('abc'), dtype='category')})
        tm.assert_frame_equal(df, expected)

        df = DataFrame([Categorical(list('abc')), Categorical(list('abd'))])
        expected = DataFrame({0: Series(list('abc'), dtype='category'),
                              1: Series(list('abd'), dtype='category')},
                             columns=[0, 1])
        tm.assert_frame_equal(df, expected)

        # mixed
        df = DataFrame([Categorical(list('abc')), list('def')])
        expected = DataFrame({0: Series(list('abc'), dtype='category'),
                              1: list('def')}, columns=[0, 1])
        tm.assert_frame_equal(df, expected)

        # invalid (shape)
        pytest.raises(ValueError,
                      lambda: DataFrame([Categorical(list('abc')),
                                         Categorical(list('abdefg'))]))

        # ndim > 1
        pytest.raises(NotImplementedError,
                      lambda: Categorical(np.array([list('abcd')])))

    def test_creation_astype(self):
        l = ["a", "b", "c", "a"]
        s = Series(l)
        exp = Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        l = [1, 2, 3, 1]
        s = Series(l)
        exp = Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        df = DataFrame({"cats": [1, 2, 3, 4, 5, 6],
                        "vals": [1, 2, 3, 4, 5, 6]})
        cats = Categorical([1, 2, 3, 4, 5, 6])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        df = DataFrame({"cats": ['a', 'b', 'b', 'a', 'a', 'd'],
                        "vals": [1, 2, 3, 4, 5, 6]})
        cats = Categorical(['a', 'b', 'b', 'a', 'a', 'd'])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        # with keywords
        l = ["a", "b", "c", "a"]
        s = Series(l)
        exp = Series(Categorical(l, ordered=True))
        res = s.astype(CategoricalDtype(None, ordered=True))
        tm.assert_series_equal(res, exp)

        exp = Series(Categorical(l, categories=list('abcdef'), ordered=True))
        res = s.astype(CategoricalDtype(list('abcdef'), ordered=True))
        tm.assert_series_equal(res, exp)
