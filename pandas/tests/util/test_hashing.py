import pytest
import datetime

from warnings import catch_warnings
import numpy as np
import pandas as pd

from pandas import DataFrame, Series, Index, MultiIndex
from pandas.util import hash_array, hash_pandas_object
from pandas.core.util.hashing import hash_tuples, hash_tuple, _hash_scalar
import pandas.util.testing as tm


class TestHashing(object):

    @pytest.fixture(params=[
        Series([1, 2, 3] * 3, dtype='int32'),
        Series([None, 2.5, 3.5] * 3, dtype='float32'),
        Series(['a', 'b', 'c'] * 3, dtype='category'),
        Series(['d', 'e', 'f'] * 3),
        Series([True, False, True] * 3),
        Series(pd.date_range('20130101', periods=9)),
        Series(pd.date_range('20130101', periods=9, tz='US/Eastern')),
        Series(pd.timedelta_range('2000', periods=9))])
    def series(self, request):
        return request.param

    def test_consistency(self):
        # check that our hash doesn't change because of a mistake
        # in the actual code; this is the ground truth
        result = hash_pandas_object(Index(['foo', 'bar', 'baz']))
        expected = Series(np.array([3600424527151052760, 1374399572096150070,
                                    477881037637427054], dtype='uint64'),
                          index=['foo', 'bar', 'baz'])
        tm.assert_series_equal(result, expected)

    def test_hash_array(self, series):
        a = series.values
        tm.assert_numpy_array_equal(hash_array(a), hash_array(a))

    def test_hash_array_mixed(self):
        result1 = hash_array(np.array([3, 4, 'All']))
        result2 = hash_array(np.array(['3', '4', 'All']))
        result3 = hash_array(np.array([3, 4, 'All'], dtype=object))
        tm.assert_numpy_array_equal(result1, result2)
        tm.assert_numpy_array_equal(result1, result3)

    @pytest.mark.parametrize('val', [5, 'foo', pd.Timestamp('20130101')])
    def test_hash_array_errors(self, val):
        msg = 'must pass a ndarray-like'
        with tm.assert_raises_regex(TypeError, msg):
            hash_array(val)

    def check_equal(self, obj, **kwargs):
        a = hash_pandas_object(obj, **kwargs)
        b = hash_pandas_object(obj, **kwargs)
        tm.assert_series_equal(a, b)

        kwargs.pop('index', None)
        a = hash_pandas_object(obj, **kwargs)
        b = hash_pandas_object(obj, **kwargs)
        tm.assert_series_equal(a, b)

    def check_not_equal_with_index(self, obj):

        # check that we are not hashing the same if
        # we include the index
        if not isinstance(obj, Index):
            a = hash_pandas_object(obj, index=True)
            b = hash_pandas_object(obj, index=False)
            if len(obj):
                assert not (a == b).all()

    def test_hash_tuples(self):
        tups = [(1, 'one'), (1, 'two'), (2, 'one')]
        result = hash_tuples(tups)
        expected = hash_pandas_object(MultiIndex.from_tuples(tups)).values
        tm.assert_numpy_array_equal(result, expected)

        result = hash_tuples(tups[0])
        assert result == expected[0]

    @pytest.mark.parametrize('tup', [
        (1, 'one'), (1, np.nan), (1.0, pd.NaT, 'A'),
        ('A', pd.Timestamp("2012-01-01"))])
    def test_hash_tuple(self, tup):
        # test equivalence between hash_tuples and hash_tuple
        result = hash_tuple(tup)
        expected = hash_tuples([tup])[0]
        assert result == expected

    @pytest.mark.parametrize('val', [
        1, 1.4, 'A', b'A', u'A', pd.Timestamp("2012-01-01"),
        pd.Timestamp("2012-01-01", tz='Europe/Brussels'),
        datetime.datetime(2012, 1, 1),
        pd.Timestamp("2012-01-01", tz='EST').to_pydatetime(),
        pd.Timedelta('1 days'), datetime.timedelta(1),
        pd.Period('2012-01-01', freq='D'), pd.Interval(0, 1),
        np.nan, pd.NaT, None])
    def test_hash_scalar(self, val):
        result = _hash_scalar(val)
        expected = hash_array(np.array([val], dtype=object), categorize=True)
        assert result[0] == expected[0]

    @pytest.mark.parametrize('val', [5, 'foo', pd.Timestamp('20130101')])
    def test_hash_tuples_err(self, val):
        msg = 'must be convertible to a list-of-tuples'
        with tm.assert_raises_regex(TypeError, msg):
            hash_tuples(val)

    def test_multiindex_unique(self):
        mi = MultiIndex.from_tuples([(118, 472), (236, 118),
                                     (51, 204), (102, 51)])
        assert mi.is_unique
        result = hash_pandas_object(mi)
        assert result.is_unique

    def test_multiindex_objects(self):
        mi = MultiIndex(levels=[['b', 'd', 'a'], [1, 2, 3]],
                        labels=[[0, 1, 0, 2], [2, 0, 0, 1]],
                        names=['col1', 'col2'])
        recons = mi._sort_levels_monotonic()

        # these are equal
        assert mi.equals(recons)
        assert Index(mi.values).equals(Index(recons.values))

        # _hashed_values and hash_pandas_object(..., index=False)
        # equivalency
        expected = hash_pandas_object(
            mi, index=False).values
        result = mi._hashed_values
        tm.assert_numpy_array_equal(result, expected)

        expected = hash_pandas_object(
            recons, index=False).values
        result = recons._hashed_values
        tm.assert_numpy_array_equal(result, expected)

        expected = mi._hashed_values
        result = recons._hashed_values

        # values should match, but in different order
        tm.assert_numpy_array_equal(np.sort(result),
                                    np.sort(expected))

    @pytest.mark.parametrize('obj', [
        Series([1, 2, 3]),
        Series([1.0, 1.5, 3.2]),
        Series([1.0, 1.5, np.nan]),
        Series([1.0, 1.5, 3.2], index=[1.5, 1.1, 3.3]),
        Series(['a', 'b', 'c']),
        Series(['a', np.nan, 'c']),
        Series(['a', None, 'c']),
        Series([True, False, True]),
        Series(),
        Index([1, 2, 3]),
        Index([True, False, True]),
        DataFrame({'x': ['a', 'b', 'c'], 'y': [1, 2, 3]}),
        DataFrame(),
        tm.makeMissingDataframe(),
        tm.makeMixedDataFrame(),
        tm.makeTimeDataFrame(),
        tm.makeTimeSeries(),
        tm.makeTimedeltaIndex(),
        tm.makePeriodIndex(),
        Series(tm.makePeriodIndex()),
        Series(pd.date_range('20130101', periods=3, tz='US/Eastern')),
        MultiIndex.from_product([range(5), ['foo', 'bar', 'baz'],
                                 pd.date_range('20130101', periods=2)]),
        MultiIndex.from_product([pd.CategoricalIndex(list('aabc')), range(3)])
    ])
    def test_hash_pandas_object(self, obj):
        self.check_equal(obj)
        self.check_not_equal_with_index(obj)

    def test_hash_pandas_object2(self, series):
        self.check_equal(series)
        self.check_not_equal_with_index(series)

    @pytest.mark.parametrize('obj', [
        Series([], dtype='float64'), Series([], dtype='object'), Index([])])
    def test_hash_pandas_empty_object(self, obj):
        # these are by-definition the same with
        # or w/o the index as the data is empty
        self.check_equal(obj)

    @pytest.mark.parametrize('s1', [
        Series(['a', 'b', 'c', 'd']),
        Series([1000, 2000, 3000, 4000]),
        Series(pd.date_range(0, periods=4))])
    @pytest.mark.parametrize('categorize', [True, False])
    def test_categorical_consistency(self, s1, categorize):
        # GH15143
        # Check that categoricals hash consistent with their values, not codes
        # This should work for categoricals of any dtype
        s2 = s1.astype('category').cat.set_categories(s1)
        s3 = s2.cat.set_categories(list(reversed(s1)))

        # These should all hash identically
        h1 = hash_pandas_object(s1, categorize=categorize)
        h2 = hash_pandas_object(s2, categorize=categorize)
        h3 = hash_pandas_object(s3, categorize=categorize)
        tm.assert_series_equal(h1, h2)
        tm.assert_series_equal(h1, h3)

    def test_categorical_with_nan_consistency(self):
        c = pd.Categorical.from_codes(
            [-1, 0, 1, 2, 3, 4],
            categories=pd.date_range('2012-01-01', periods=5, name='B'))
        expected = hash_array(c, categorize=False)
        c = pd.Categorical.from_codes(
            [-1, 0],
            categories=[pd.Timestamp('2012-01-01')])
        result = hash_array(c, categorize=False)
        assert result[0] in expected
        assert result[1] in expected

    def test_pandas_errors(self):
        with pytest.raises(TypeError):
            hash_pandas_object(pd.Timestamp('20130101'))

        with catch_warnings(record=True):
            obj = tm.makePanel()

        with pytest.raises(TypeError):
            hash_pandas_object(obj)

    def test_hash_keys(self):
        # using different hash keys, should have different hashes
        # for the same data

        # this only matters for object dtypes
        obj = Series(list('abc'))
        a = hash_pandas_object(obj, hash_key='9876543210123456')
        b = hash_pandas_object(obj, hash_key='9876543210123465')
        assert (a != b).all()

    def test_invalid_key(self):
        # this only matters for object dtypes
        msg = 'key should be a 16-byte string encoded'
        with tm.assert_raises_regex(ValueError, msg):
            hash_pandas_object(Series(list('abc')), hash_key='foo')

    def test_alread_encoded(self):
        # if already encoded then ok

        obj = Series(list('abc')).str.encode('utf8')
        self.check_equal(obj)

    def test_alternate_encoding(self):

        obj = Series(list('abc'))
        self.check_equal(obj, encoding='ascii')

    @pytest.mark.parametrize('l_exp', range(8))
    @pytest.mark.parametrize('l_add', [0, 1])
    def test_same_len_hash_collisions(self, l_exp, l_add):
        length = 2**(l_exp + 8) + l_add
        s = tm.rands_array(length, 2)
        result = hash_array(s, 'utf8')
        assert not result[0] == result[1]

    def test_hash_collisions(self):

        # hash collisions are bad
        # https://github.com/pandas-dev/pandas/issues/14711#issuecomment-264885726
        L = ['Ingrid-9Z9fKIZmkO7i7Cn51Li34pJm44fgX6DYGBNj3VPlOH50m7HnBlPxfIwFMrcNJNMP6PSgLmwWnInciMWrCSAlLEvt7JkJl4IxiMrVbXSa8ZQoVaq5xoQPjltuJEfwdNlO6jo8qRRHvD8sBEBMQASrRa6TsdaPTPCBo3nwIBpE7YzzmyH0vMBhjQZLx1aCT7faSEx7PgFxQhHdKFWROcysamgy9iVj8DO2Fmwg1NNl93rIAqC3mdqfrCxrzfvIY8aJdzin2cHVzy3QUJxZgHvtUtOLxoqnUHsYbNTeq0xcLXpTZEZCxD4PGubIuCNf32c33M7HFsnjWSEjE2yVdWKhmSVodyF8hFYVmhYnMCztQnJrt3O8ZvVRXd5IKwlLexiSp4h888w7SzAIcKgc3g5XQJf6MlSMftDXm9lIsE1mJNiJEv6uY6pgvC3fUPhatlR5JPpVAHNSbSEE73MBzJrhCAbOLXQumyOXigZuPoME7QgJcBalliQol7YZ9',  # noqa
             'Tim-b9MddTxOWW2AT1Py6vtVbZwGAmYCjbp89p8mxsiFoVX4FyDOF3wFiAkyQTUgwg9sVqVYOZo09Dh1AzhFHbgij52ylF0SEwgzjzHH8TGY8Lypart4p4onnDoDvVMBa0kdthVGKl6K0BDVGzyOXPXKpmnMF1H6rJzqHJ0HywfwS4XYpVwlAkoeNsiicHkJUFdUAhG229INzvIAiJuAHeJDUoyO4DCBqtoZ5TDend6TK7Y914yHlfH3g1WZu5LksKv68VQHJriWFYusW5e6ZZ6dKaMjTwEGuRgdT66iU5nqWTHRH8WSzpXoCFwGcTOwyuqPSe0fTe21DVtJn1FKj9F9nEnR9xOvJUO7E0piCIF4Ad9yAIDY4DBimpsTfKXCu1vdHpKYerzbndfuFe5AhfMduLYZJi5iAw8qKSwR5h86ttXV0Mc0QmXz8dsRvDgxjXSmupPxBggdlqUlC828hXiTPD7am0yETBV0F3bEtvPiNJfremszcV8NcqAoARMe']  # noqa

        # these should be different!
        result1 = hash_array(np.asarray(L[0:1], dtype=object), 'utf8')
        expected1 = np.array([14963968704024874985], dtype=np.uint64)
        tm.assert_numpy_array_equal(result1, expected1)

        result2 = hash_array(np.asarray(L[1:2], dtype=object), 'utf8')
        expected2 = np.array([16428432627716348016], dtype=np.uint64)
        tm.assert_numpy_array_equal(result2, expected2)

        result = hash_array(np.asarray(L, dtype=object), 'utf8')
        tm.assert_numpy_array_equal(
            result, np.concatenate([expected1, expected2], axis=0))
