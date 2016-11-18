# -*- coding: utf-8 -*-

# TODO(wesm): fix long line flake8 issues
# flake8: noqa

import pandas.util.testing as tm
from pandas.indexes.api import Index, CategoricalIndex
from .common import Base

from pandas.compat import range, PY3

import numpy as np

from pandas import Categorical, compat, notnull
from pandas.util.testing import assert_almost_equal
import pandas.core.config as cf
import pandas as pd

if PY3:
    unicode = lambda x: x


class TestCategoricalIndex(Base, tm.TestCase):
    _holder = CategoricalIndex

    def setUp(self):
        self.indices = dict(catIndex=tm.makeCategoricalIndex(100))
        self.setup_indices()

    def create_index(self, categories=None, ordered=False):
        if categories is None:
            categories = list('cab')
        return CategoricalIndex(
            list('aabbca'), categories=categories, ordered=ordered)

    def test_construction(self):

        ci = self.create_index(categories=list('abcd'))
        categories = ci.categories

        result = Index(ci)
        tm.assert_index_equal(result, ci, exact=True)
        self.assertFalse(result.ordered)

        result = Index(ci.values)
        tm.assert_index_equal(result, ci, exact=True)
        self.assertFalse(result.ordered)

        # empty
        result = CategoricalIndex(categories=categories)
        self.assert_index_equal(result.categories, Index(categories))
        tm.assert_numpy_array_equal(result.codes, np.array([], dtype='int8'))
        self.assertFalse(result.ordered)

        # passing categories
        result = CategoricalIndex(list('aabbca'), categories=categories)
        self.assert_index_equal(result.categories, Index(categories))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, 2, 0], dtype='int8'))

        c = pd.Categorical(list('aabbca'))
        result = CategoricalIndex(c)
        self.assert_index_equal(result.categories, Index(list('abc')))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, 2, 0], dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(c, categories=categories)
        self.assert_index_equal(result.categories, Index(categories))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, 2, 0], dtype='int8'))
        self.assertFalse(result.ordered)

        ci = CategoricalIndex(c, categories=list('abcd'))
        result = CategoricalIndex(ci)
        self.assert_index_equal(result.categories, Index(categories))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, 2, 0], dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(ci, categories=list('ab'))
        self.assert_index_equal(result.categories, Index(list('ab')))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, -1, 0],
                                             dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(ci, categories=list('ab'), ordered=True)
        self.assert_index_equal(result.categories, Index(list('ab')))
        tm.assert_numpy_array_equal(result.codes,
                                    np.array([0, 0, 1, 1, -1, 0],
                                             dtype='int8'))
        self.assertTrue(result.ordered)

        # turn me to an Index
        result = Index(np.array(ci))
        self.assertIsInstance(result, Index)
        self.assertNotIsInstance(result, CategoricalIndex)

    def test_construction_with_dtype(self):

        # specify dtype
        ci = self.create_index(categories=list('abc'))

        result = Index(np.array(ci), dtype='category')
        tm.assert_index_equal(result, ci, exact=True)

        result = Index(np.array(ci).tolist(), dtype='category')
        tm.assert_index_equal(result, ci, exact=True)

        # these are generally only equal when the categories are reordered
        ci = self.create_index()

        result = Index(
            np.array(ci), dtype='category').reorder_categories(ci.categories)
        tm.assert_index_equal(result, ci, exact=True)

        # make sure indexes are handled
        expected = CategoricalIndex([0, 1, 2], categories=[0, 1, 2],
                                    ordered=True)
        idx = Index(range(3))
        result = CategoricalIndex(idx, categories=idx, ordered=True)
        tm.assert_index_equal(result, expected, exact=True)

    def test_disallow_set_ops(self):

        # GH 10039
        # set ops (+/-) raise TypeError
        idx = pd.Index(pd.Categorical(['a', 'b']))

        self.assertRaises(TypeError, lambda: idx - idx)
        self.assertRaises(TypeError, lambda: idx + idx)
        self.assertRaises(TypeError, lambda: idx - ['a', 'b'])
        self.assertRaises(TypeError, lambda: idx + ['a', 'b'])
        self.assertRaises(TypeError, lambda: ['a', 'b'] - idx)
        self.assertRaises(TypeError, lambda: ['a', 'b'] + idx)

    def test_method_delegation(self):

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.set_categories(list('cab'))
        tm.assert_index_equal(result, CategoricalIndex(
            list('aabbca'), categories=list('cab')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.rename_categories(list('efg'))
        tm.assert_index_equal(result, CategoricalIndex(
            list('ffggef'), categories=list('efg')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.add_categories(['d'])
        tm.assert_index_equal(result, CategoricalIndex(
            list('aabbca'), categories=list('cabd')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.remove_categories(['c'])
        tm.assert_index_equal(result, CategoricalIndex(
            list('aabb') + [np.nan] + ['a'], categories=list('ab')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.as_unordered()
        tm.assert_index_equal(result, ci)

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.as_ordered()
        tm.assert_index_equal(result, CategoricalIndex(
            list('aabbca'), categories=list('cabdef'), ordered=True))

        # invalid
        self.assertRaises(ValueError, lambda: ci.set_categories(
            list('cab'), inplace=True))

    def test_contains(self):

        ci = self.create_index(categories=list('cabdef'))

        self.assertTrue('a' in ci)
        self.assertTrue('z' not in ci)
        self.assertTrue('e' not in ci)
        self.assertTrue(np.nan not in ci)

        # assert codes NOT in index
        self.assertFalse(0 in ci)
        self.assertFalse(1 in ci)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            ci = CategoricalIndex(
                list('aabbca'), categories=list('cabdef') + [np.nan])
        self.assertFalse(np.nan in ci)

        ci = CategoricalIndex(
            list('aabbca') + [np.nan], categories=list('cabdef'))
        self.assertTrue(np.nan in ci)

    def test_min_max(self):

        ci = self.create_index(ordered=False)
        self.assertRaises(TypeError, lambda: ci.min())
        self.assertRaises(TypeError, lambda: ci.max())

        ci = self.create_index(ordered=True)

        self.assertEqual(ci.min(), 'c')
        self.assertEqual(ci.max(), 'b')

    def test_map(self):
        ci = pd.CategoricalIndex(list('ABABC'), categories=list('CBA'),
                                 ordered=True)
        result = ci.map(lambda x: x.lower())
        exp = pd.Categorical(list('ababc'), categories=list('cba'),
                             ordered=True)
        tm.assert_categorical_equal(result, exp)

        ci = pd.CategoricalIndex(list('ABABC'), categories=list('BAC'),
                                 ordered=False, name='XXX')
        result = ci.map(lambda x: x.lower())
        exp = pd.Categorical(list('ababc'), categories=list('bac'),
                             ordered=False)
        tm.assert_categorical_equal(result, exp)

        tm.assert_numpy_array_equal(ci.map(lambda x: 1),
                                    np.array([1] * 5, dtype=np.int64))

        # change categories dtype
        ci = pd.CategoricalIndex(list('ABABC'), categories=list('BAC'),
                                 ordered=False)
        def f(x):
            return {'A': 10, 'B': 20, 'C': 30}.get(x)

        result = ci.map(f)
        exp = pd.Categorical([10, 20, 10, 20, 30], categories=[20, 10, 30],
                             ordered=False)
        tm.assert_categorical_equal(result, exp)

    def test_where(self):
        i = self.create_index()
        result = i.where(notnull(i))
        expected = i
        tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = pd.CategoricalIndex([np.nan, np.nan] + i[2:].tolist(),
                                 categories=i.categories)
        result = i.where(notnull(i2))
        expected = i2
        tm.assert_index_equal(result, expected)

    def test_append(self):

        ci = self.create_index()
        categories = ci.categories

        # append cats with the same categories
        result = ci[:3].append(ci[3:])
        tm.assert_index_equal(result, ci, exact=True)

        foos = [ci[:1], ci[1:3], ci[3:]]
        result = foos[0].append(foos[1:])
        tm.assert_index_equal(result, ci, exact=True)

        # empty
        result = ci.append([])
        tm.assert_index_equal(result, ci, exact=True)

        # appending with different categories or reoreded is not ok
        self.assertRaises(
            TypeError,
            lambda: ci.append(ci.values.set_categories(list('abcd'))))
        self.assertRaises(
            TypeError,
            lambda: ci.append(ci.values.reorder_categories(list('abc'))))

        # with objects
        result = ci.append(Index(['c', 'a']))
        expected = CategoricalIndex(list('aabbcaca'), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # invalid objects
        self.assertRaises(TypeError, lambda: ci.append(Index(['a', 'd'])))

        # GH14298 - if base object is not categorical -> coerce to object
        result = Index(['c', 'a']).append(ci)
        expected = Index(list('caaabbca'))
        tm.assert_index_equal(result, expected, exact=True)

    def test_insert(self):

        ci = self.create_index()
        categories = ci.categories

        # test 0th element
        result = ci.insert(0, 'a')
        expected = CategoricalIndex(list('aaabbca'), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # test Nth element that follows Python list behavior
        result = ci.insert(-1, 'a')
        expected = CategoricalIndex(list('aabbcaa'), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # test empty
        result = CategoricalIndex(categories=categories).insert(0, 'a')
        expected = CategoricalIndex(['a'], categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # invalid
        self.assertRaises(TypeError, lambda: ci.insert(0, 'd'))

    def test_delete(self):

        ci = self.create_index()
        categories = ci.categories

        result = ci.delete(0)
        expected = CategoricalIndex(list('abbca'), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        result = ci.delete(-1)
        expected = CategoricalIndex(list('aabbc'), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = ci.delete(10)

    def test_astype(self):

        ci = self.create_index()
        result = ci.astype('category')
        tm.assert_index_equal(result, ci, exact=True)

        result = ci.astype(object)
        self.assert_index_equal(result, Index(np.array(ci)))

        # this IS equal, but not the same class
        self.assertTrue(result.equals(ci))
        self.assertIsInstance(result, Index)
        self.assertNotIsInstance(result, CategoricalIndex)

    def test_reindex_base(self):

        # determined by cat ordering
        idx = self.create_index()
        expected = np.array([4, 0, 1, 5, 2, 3], dtype=np.intp)

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with tm.assertRaisesRegexp(ValueError, 'Invalid fill method'):
            idx.get_indexer(idx, method='invalid')

    def test_reindexing(self):

        ci = self.create_index()
        oidx = Index(np.array(ci))

        for n in [1, 2, 5, len(ci)]:
            finder = oidx[np.random.randint(0, len(ci), size=n)]
            expected = oidx.get_indexer_non_unique(finder)[0]

            actual = ci.get_indexer(finder)
            tm.assert_numpy_array_equal(expected.values, actual, check_dtype=False)

    def test_reindex_dtype(self):
        c = CategoricalIndex(['a', 'b', 'c', 'a'])
        res, indexer = c.reindex(['a', 'c'])
        tm.assert_index_equal(res, Index(['a', 'a', 'c']), exact=True)
        tm.assert_numpy_array_equal(indexer,
                                    np.array([0, 3, 2], dtype=np.int64))

        c = CategoricalIndex(['a', 'b', 'c', 'a'])
        res, indexer = c.reindex(Categorical(['a', 'c']))

        exp = CategoricalIndex(['a', 'a', 'c'], categories=['a', 'c'])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer,
                                    np.array([0, 3, 2], dtype=np.int64))

        c = CategoricalIndex(['a', 'b', 'c', 'a'],
                             categories=['a', 'b', 'c', 'd'])
        res, indexer = c.reindex(['a', 'c'])
        exp = Index(['a', 'a', 'c'], dtype='object')
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer,
                                    np.array([0, 3, 2], dtype=np.int64))

        c = CategoricalIndex(['a', 'b', 'c', 'a'],
                             categories=['a', 'b', 'c', 'd'])
        res, indexer = c.reindex(Categorical(['a', 'c']))
        exp = CategoricalIndex(['a', 'a', 'c'], categories=['a', 'c'])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer,
                                    np.array([0, 3, 2], dtype=np.int64))

    def test_duplicates(self):

        idx = CategoricalIndex([0, 0, 0], name='foo')
        self.assertFalse(idx.is_unique)
        self.assertTrue(idx.has_duplicates)

        expected = CategoricalIndex([0], name='foo')
        self.assert_index_equal(idx.drop_duplicates(), expected)
        self.assert_index_equal(idx.unique(), expected)

    def test_get_indexer(self):

        idx1 = CategoricalIndex(list('aabcde'), categories=list('edabc'))
        idx2 = CategoricalIndex(list('abf'))

        for indexer in [idx2, list('abf'), Index(list('abf'))]:
            r1 = idx1.get_indexer(idx2)
            assert_almost_equal(r1, np.array([0, 1, 2, -1], dtype=np.intp))

        self.assertRaises(NotImplementedError,
                          lambda: idx2.get_indexer(idx1, method='pad'))
        self.assertRaises(NotImplementedError,
                          lambda: idx2.get_indexer(idx1, method='backfill'))
        self.assertRaises(NotImplementedError,
                          lambda: idx2.get_indexer(idx1, method='nearest'))

    def test_get_loc(self):
        # GH 12531
        cidx1 = CategoricalIndex(list('abcde'), categories=list('edabc'))
        idx1 = Index(list('abcde'))
        self.assertEqual(cidx1.get_loc('a'), idx1.get_loc('a'))
        self.assertEqual(cidx1.get_loc('e'), idx1.get_loc('e'))

        for i in [cidx1, idx1]:
            with tm.assertRaises(KeyError):
                i.get_loc('NOT-EXIST')

        # non-unique
        cidx2 = CategoricalIndex(list('aacded'), categories=list('edabc'))
        idx2 = Index(list('aacded'))

        # results in bool array
        res = cidx2.get_loc('d')
        self.assert_numpy_array_equal(res, idx2.get_loc('d'))
        self.assert_numpy_array_equal(res, np.array([False, False, False,
                                                     True, False, True]))
        # unique element results in scalar
        res = cidx2.get_loc('e')
        self.assertEqual(res, idx2.get_loc('e'))
        self.assertEqual(res, 4)

        for i in [cidx2, idx2]:
            with tm.assertRaises(KeyError):
                i.get_loc('NOT-EXIST')

        # non-unique, slicable
        cidx3 = CategoricalIndex(list('aabbb'), categories=list('abc'))
        idx3 = Index(list('aabbb'))

        # results in slice
        res = cidx3.get_loc('a')
        self.assertEqual(res, idx3.get_loc('a'))
        self.assertEqual(res, slice(0, 2, None))

        res = cidx3.get_loc('b')
        self.assertEqual(res, idx3.get_loc('b'))
        self.assertEqual(res, slice(2, 5, None))

        for i in [cidx3, idx3]:
            with tm.assertRaises(KeyError):
                i.get_loc('c')

    def test_repr_roundtrip(self):

        ci = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        str(ci)
        tm.assert_index_equal(eval(repr(ci)), ci, exact=True)

        # formatting
        if PY3:
            str(ci)
        else:
            compat.text_type(ci)

        # long format
        # this is not reprable
        ci = CategoricalIndex(np.random.randint(0, 5, size=100))
        if PY3:
            str(ci)
        else:
            compat.text_type(ci)

    def test_isin(self):

        ci = CategoricalIndex(
            list('aabca') + [np.nan], categories=['c', 'a', 'b'])
        tm.assert_numpy_array_equal(
            ci.isin(['c']),
            np.array([False, False, False, True, False, False]))
        tm.assert_numpy_array_equal(
            ci.isin(['c', 'a', 'b']), np.array([True] * 5 + [False]))
        tm.assert_numpy_array_equal(
            ci.isin(['c', 'a', 'b', np.nan]), np.array([True] * 6))

        # mismatched categorical -> coerced to ndarray so doesn't matter
        tm.assert_numpy_array_equal(
            ci.isin(ci.set_categories(list('abcdefghi'))), np.array([True] *
                                                                    6))
        tm.assert_numpy_array_equal(
            ci.isin(ci.set_categories(list('defghi'))),
            np.array([False] * 5 + [True]))

    def test_identical(self):

        ci1 = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        ci2 = CategoricalIndex(['a', 'b'], categories=['a', 'b', 'c'],
                               ordered=True)
        self.assertTrue(ci1.identical(ci1))
        self.assertTrue(ci1.identical(ci1.copy()))
        self.assertFalse(ci1.identical(ci2))

    def test_ensure_copied_data(self):
        # Check the "copy" argument of each Index.__new__ is honoured
        # GH12309
        # Must be tested separately from other indexes because
        # self.value is not an ndarray
        _base = lambda ar : ar if ar.base is None else ar.base
        for index in self.indices.values():
            result = CategoricalIndex(index.values, copy=True)
            tm.assert_index_equal(index, result)
            self.assertIsNot(_base(index.values), _base(result.values))

            result = CategoricalIndex(index.values, copy=False)
            self.assertIs(_base(index.values), _base(result.values))

    def test_equals_categorical(self):

        ci1 = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        ci2 = CategoricalIndex(['a', 'b'], categories=['a', 'b', 'c'],
                               ordered=True)

        self.assertTrue(ci1.equals(ci1))
        self.assertFalse(ci1.equals(ci2))
        self.assertTrue(ci1.equals(ci1.astype(object)))
        self.assertTrue(ci1.astype(object).equals(ci1))

        self.assertTrue((ci1 == ci1).all())
        self.assertFalse((ci1 != ci1).all())
        self.assertFalse((ci1 > ci1).all())
        self.assertFalse((ci1 < ci1).all())
        self.assertTrue((ci1 <= ci1).all())
        self.assertTrue((ci1 >= ci1).all())

        self.assertFalse((ci1 == 1).all())
        self.assertTrue((ci1 == Index(['a', 'b'])).all())
        self.assertTrue((ci1 == ci1.values).all())

        # invalid comparisons
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            ci1 == Index(['a', 'b', 'c'])
        self.assertRaises(TypeError, lambda: ci1 == ci2)
        self.assertRaises(
            TypeError, lambda: ci1 == Categorical(ci1.values, ordered=False))
        self.assertRaises(
            TypeError,
            lambda: ci1 == Categorical(ci1.values, categories=list('abc')))

        # tests
        # make sure that we are testing for category inclusion properly
        ci = CategoricalIndex(list('aabca'), categories=['c', 'a', 'b'])
        self.assertFalse(ci.equals(list('aabca')))
        self.assertFalse(ci.equals(CategoricalIndex(list('aabca'))))
        self.assertTrue(ci.equals(ci.copy()))

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            ci = CategoricalIndex(list('aabca'),
                                  categories=['c', 'a', 'b', np.nan])
        self.assertFalse(ci.equals(list('aabca')))
        self.assertFalse(ci.equals(CategoricalIndex(list('aabca'))))
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertTrue(ci.equals(ci.copy()))

        ci = CategoricalIndex(list('aabca') + [np.nan],
                              categories=['c', 'a', 'b'])
        self.assertFalse(ci.equals(list('aabca')))
        self.assertFalse(ci.equals(CategoricalIndex(list('aabca'))))
        self.assertTrue(ci.equals(ci.copy()))

        ci = CategoricalIndex(list('aabca') + [np.nan],
                              categories=['c', 'a', 'b'])
        self.assertFalse(ci.equals(list('aabca') + [np.nan]))
        self.assertFalse(ci.equals(CategoricalIndex(list('aabca') + [np.nan])))
        self.assertTrue(ci.equals(ci.copy()))

    def test_string_categorical_index_repr(self):
        # short
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'])
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc'], categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc'], categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'] * 10)
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb',
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
                  u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
                  u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc',
                  u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
                 categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category')"""

            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'] * 100)
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  ...
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category', length=300)"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
                  u'ccc', u'a',
                  ...
                  u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
                  u'bb', u'ccc'],
                 categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category', length=300)"""

            self.assertEqual(unicode(idx), expected)

        # larger categories
        idx = pd.CategoricalIndex(list('abcdefghijklmmo'))
        if PY3:
            expected = u"""CategoricalIndex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                  'm', 'm', 'o'],
                 categories=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', ...], ordered=False, dtype='category')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'b', u'c', u'd', u'e', u'f', u'g', u'h', u'i', u'j',
                  u'k', u'l', u'm', u'm', u'o'],
                 categories=[u'a', u'b', u'c', u'd', u'e', u'f', u'g', u'h', ...], ordered=False, dtype='category')"""

            self.assertEqual(unicode(idx), expected)

        # short
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'])
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう'], categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 10)
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう',
                  u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""

            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 100)
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ',
                  ...
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category', length=300)"""

            self.assertEqual(unicode(idx), expected)

        # larger categories
        idx = pd.CategoricalIndex(list(u'あいうえおかきくけこさしすせそ'))
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ', 'さ', 'し',
                  'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', u'け', u'こ',
                  u'さ', u'し', u'す', u'せ', u'そ'],
                 categories=[u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', ...], ordered=False, dtype='category')"""

            self.assertEqual(unicode(idx), expected)

        # Emable Unicode option -----------------------------------------
        with cf.option_context('display.unicode.east_asian_width', True):

            # short
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'])
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう'], categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(unicode(idx), expected)

            # multiple lines
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 10)
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""

                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""

                self.assertEqual(unicode(idx), expected)

            # truncated
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 100)
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""

                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ',
                  ...
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category', length=300)"""

                self.assertEqual(unicode(idx), expected)

            # larger categories
            idx = pd.CategoricalIndex(list(u'あいうえおかきくけこさしすせそ'))
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ',
                  'さ', 'し', 'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""

                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く',
                  u'け', u'こ', u'さ', u'し', u'す', u'せ', u'そ'],
                 categories=[u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', ...], ordered=False, dtype='category')"""

                self.assertEqual(unicode(idx), expected)

    def test_fillna_categorical(self):
        # GH 11343
        idx = CategoricalIndex([1.0, np.nan, 3.0, 1.0], name='x')
        # fill by value in categories
        exp = CategoricalIndex([1.0, 1.0, 3.0, 1.0], name='x')
        self.assert_index_equal(idx.fillna(1.0), exp)

        # fill by value not in categories raises ValueError
        with tm.assertRaisesRegexp(ValueError,
                                   'fill value must be in categories'):
            idx.fillna(2.0)

    def test_take_fill_value(self):
        # GH 12631

        # numeric category
        idx = pd.CategoricalIndex([1, 2, 3], name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.CategoricalIndex([2, 1, 3], name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.CategoricalIndex([2, 1, np.nan], categories=[1, 2, 3],
                                       name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.CategoricalIndex([2, 1, 3], name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # object category
        idx = pd.CategoricalIndex(list('CBA'), categories=list('ABC'),
                                  ordered=True, name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.CategoricalIndex(list('BCA'), categories=list('ABC'),
                                       ordered=True, name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.CategoricalIndex(['B', 'C', np.nan],
                                       categories=list('ABC'), ordered=True,
                                       name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.CategoricalIndex(list('BCA'), categories=list('ABC'),
                                       ordered=True, name='xxx')
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_take_fill_value_datetime(self):

        # datetime category
        idx = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               name='xxx')
        idx = pd.CategoricalIndex(idx)
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        expected = pd.CategoricalIndex(expected)
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                    name='xxx')
        exp_cats = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'])
        expected = pd.CategoricalIndex(expected, categories=exp_cats)
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        expected = pd.CategoricalIndex(expected)
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_take_invalid_kwargs(self):
        idx = pd.CategoricalIndex([1, 2, 3], name='foo')
        indices = [1, 0, -1]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        tm.assertRaisesRegexp(TypeError, msg, idx.take,
                              indices, foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, mode='clip')
