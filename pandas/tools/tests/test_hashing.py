import numpy as np
import pandas as pd

from pandas import DataFrame, Series, Index
from pandas.tools.hashing import hash_array, hash_pandas_object
import pandas.util.testing as tm


class TestHashing(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.df = DataFrame(
            {'i32': np.array([1, 2, 3] * 3, dtype='int32'),
             'f32': np.array([None, 2.5, 3.5] * 3, dtype='float32'),
             'cat': Series(['a', 'b', 'c'] * 3).astype('category'),
             'obj': Series(['d', 'e', 'f'] * 3),
             'bool': np.array([True, False, True] * 3),
             'dt': Series(pd.date_range('20130101', periods=9)),
             'dt_tz': Series(pd.date_range('20130101', periods=9,
                                           tz='US/Eastern')),
             'td': Series(pd.timedelta_range('2000', periods=9))})

    def test_consistency(self):
        # check that our hash doesn't change because of a mistake
        # in the actual code; this is the ground truth
        result = hash_pandas_object(Index(['foo', 'bar', 'baz']))
        expected = Series(np.array([3600424527151052760, 1374399572096150070,
                                    477881037637427054], dtype='uint64'),
                          index=['foo', 'bar', 'baz'])
        tm.assert_series_equal(result, expected)

    def test_hash_array(self):
        for name, s in self.df.iteritems():
            a = s.values
            tm.assert_numpy_array_equal(hash_array(a), hash_array(a))

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
            self.assertFalse((a == b).all())

    def test_hash_pandas_object(self):

        for obj in [Series([1, 2, 3]),
                    Series([1.0, 1.5, 3.2]),
                    Series([1.0, 1.5, np.nan]),
                    Series([1.0, 1.5, 3.2], index=[1.5, 1.1, 3.3]),
                    Series(['a', 'b', 'c']),
                    Series(['a', np.nan, 'c']),
                    Series(['a', None, 'c']),
                    Series([True, False, True]),
                    Index([1, 2, 3]),
                    Index([True, False, True]),
                    DataFrame({'x': ['a', 'b', 'c'], 'y': [1, 2, 3]}),
                    tm.makeMissingDataframe(),
                    tm.makeMixedDataFrame(),
                    tm.makeTimeDataFrame(),
                    tm.makeTimeSeries(),
                    tm.makeTimedeltaIndex()]:
            self.check_equal(obj)
            self.check_not_equal_with_index(obj)

    def test_hash_pandas_object2(self):
        for name, s in self.df.iteritems():
            self.check_equal(s)
            self.check_not_equal_with_index(s)

    def test_hash_pandas_empty_object(self):
        for obj in [Series([], dtype='float64'),
                    Series([], dtype='object'),
                    Index([])]:
            self.check_equal(obj)

            # these are by-definition the same with
            # or w/o the index as the data is empty

    def test_errors(self):

        for obj in [pd.Timestamp('20130101'), tm.makePanel()]:
            def f():
                hash_pandas_object(f)

            self.assertRaises(TypeError, f)

    def test_hash_keys(self):
        # using different hash keys, should have different hashes
        # for the same data

        # this only matters for object dtypes
        obj = Series(list('abc'))
        a = hash_pandas_object(obj, hash_key='9876543210123456')
        b = hash_pandas_object(obj, hash_key='9876543210123465')
        self.assertTrue((a != b).all())

    def test_invalid_key(self):
        # this only matters for object dtypes
        def f():
            hash_pandas_object(Series(list('abc')), hash_key='foo')
        self.assertRaises(ValueError, f)

    def test_unsupported_objects(self):

        # mixed objects are not supported
        obj = Series(['1', 2, 3])

        def f():
            hash_pandas_object(obj)
        self.assertRaises(TypeError, f)

        # MultiIndex are represented as tuples
        obj = Series([1, 2, 3], index=pd.MultiIndex.from_tuples(
            [('a', 1), ('a', 2), ('b', 1)]))

        def f():
            hash_pandas_object(obj)
        self.assertRaises(TypeError, f)

    def test_alread_encoded(self):
        # if already encoded then ok

        obj = Series(list('abc')).str.encode('utf8')
        self.check_equal(obj)

    def test_alternate_encoding(self):

        obj = Series(list('abc'))
        self.check_equal(obj, encoding='ascii')

    def test_same_len_hash_collisions(self):

        for l in range(8):
            length = 2**(l + 8) + 1
            s = tm.rands_array(length, 2)
            result = hash_array(s, 'utf8')
            self.assertFalse(result[0] == result[1])

        for l in range(8):
            length = 2**(l + 8)
            s = tm.rands_array(length, 2)
            result = hash_array(s, 'utf8')
            self.assertFalse(result[0] == result[1])

    def test_hash_collisions(self):

        # hash collisions are bad
        # https://github.com/pandas-dev/pandas/issues/14711#issuecomment-264885726
        L = ['Ingrid-9Z9fKIZmkO7i7Cn51Li34pJm44fgX6DYGBNj3VPlOH50m7HnBlPxfIwFMrcNJNMP6PSgLmwWnInciMWrCSAlLEvt7JkJl4IxiMrVbXSa8ZQoVaq5xoQPjltuJEfwdNlO6jo8qRRHvD8sBEBMQASrRa6TsdaPTPCBo3nwIBpE7YzzmyH0vMBhjQZLx1aCT7faSEx7PgFxQhHdKFWROcysamgy9iVj8DO2Fmwg1NNl93rIAqC3mdqfrCxrzfvIY8aJdzin2cHVzy3QUJxZgHvtUtOLxoqnUHsYbNTeq0xcLXpTZEZCxD4PGubIuCNf32c33M7HFsnjWSEjE2yVdWKhmSVodyF8hFYVmhYnMCztQnJrt3O8ZvVRXd5IKwlLexiSp4h888w7SzAIcKgc3g5XQJf6MlSMftDXm9lIsE1mJNiJEv6uY6pgvC3fUPhatlR5JPpVAHNSbSEE73MBzJrhCAbOLXQumyOXigZuPoME7QgJcBalliQol7YZ9',  # noqa
             'Tim-b9MddTxOWW2AT1Py6vtVbZwGAmYCjbp89p8mxsiFoVX4FyDOF3wFiAkyQTUgwg9sVqVYOZo09Dh1AzhFHbgij52ylF0SEwgzjzHH8TGY8Lypart4p4onnDoDvVMBa0kdthVGKl6K0BDVGzyOXPXKpmnMF1H6rJzqHJ0HywfwS4XYpVwlAkoeNsiicHkJUFdUAhG229INzvIAiJuAHeJDUoyO4DCBqtoZ5TDend6TK7Y914yHlfH3g1WZu5LksKv68VQHJriWFYusW5e6ZZ6dKaMjTwEGuRgdT66iU5nqWTHRH8WSzpXoCFwGcTOwyuqPSe0fTe21DVtJn1FKj9F9nEnR9xOvJUO7E0piCIF4Ad9yAIDY4DBimpsTfKXCu1vdHpKYerzbndfuFe5AhfMduLYZJi5iAw8qKSwR5h86ttXV0Mc0QmXz8dsRvDgxjXSmupPxBggdlqUlC828hXiTPD7am0yETBV0F3bEtvPiNJfremszcV8NcqAoARMe']  # noqa

        # these should be different!
        result1 = hash_array(np.asarray(L[0:1], dtype=object), 'utf8')
        expected1 = np.array([14963968704024874985], dtype=np.uint64)
        self.assert_numpy_array_equal(result1, expected1)

        result2 = hash_array(np.asarray(L[1:2], dtype=object), 'utf8')
        expected2 = np.array([16428432627716348016], dtype=np.uint64)
        self.assert_numpy_array_equal(result2, expected2)

        result = hash_array(np.asarray(L, dtype=object), 'utf8')
        self.assert_numpy_array_equal(
            result, np.concatenate([expected1, expected2], axis=0))
