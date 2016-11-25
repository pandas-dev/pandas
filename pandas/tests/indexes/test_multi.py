# -*- coding: utf-8 -*-

from datetime import timedelta
from itertools import product
import nose
import re
import warnings

from pandas import (DataFrame, date_range, period_range, MultiIndex, Index,
                    CategoricalIndex, compat)
from pandas.core.common import PerformanceWarning
from pandas.indexes.base import InvalidIndexError
from pandas.compat import range, lrange, u, PY3, long, lzip

import numpy as np

from pandas.util.testing import (assert_almost_equal, assertRaises,
                                 assertRaisesRegexp, assert_copy)

import pandas.util.testing as tm

import pandas as pd
from pandas.lib import Timestamp

from .common import Base


class TestMultiIndex(Base, tm.TestCase):
    _holder = MultiIndex
    _multiprocess_can_split_ = True
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setUp(self):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.indices = dict(index=MultiIndex(levels=[major_axis, minor_axis],
                                             labels=[major_labels, minor_labels
                                                     ], names=self.index_names,
                                             verify_integrity=False))
        self.setup_indices()

    def create_index(self):
        return self.index

    def test_boolean_context_compat2(self):

        # boolean context compat
        # GH7897
        i1 = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        i2 = MultiIndex.from_tuples([('A', 1), ('A', 3)])
        common = i1.intersection(i2)

        def f():
            if common:
                pass

        tm.assertRaisesRegexp(ValueError, 'The truth value of a', f)

    def test_labels_dtypes(self):

        # GH 8456
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        self.assertTrue(i.labels[0].dtype == 'int8')
        self.assertTrue(i.labels[1].dtype == 'int8')

        i = MultiIndex.from_product([['a'], range(40)])
        self.assertTrue(i.labels[1].dtype == 'int8')
        i = MultiIndex.from_product([['a'], range(400)])
        self.assertTrue(i.labels[1].dtype == 'int16')
        i = MultiIndex.from_product([['a'], range(40000)])
        self.assertTrue(i.labels[1].dtype == 'int32')

        i = pd.MultiIndex.from_product([['a'], range(1000)])
        self.assertTrue((i.labels[0] >= 0).all())
        self.assertTrue((i.labels[1] >= 0).all())

    def test_where(self):
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])

        def f():
            i.where(True)

        self.assertRaises(NotImplementedError, f)

    def test_repeat(self):
        reps = 2
        numbers = [1, 2, 3]
        names = np.array(['foo', 'bar'])

        m = MultiIndex.from_product([
            numbers, names], names=names)
        expected = MultiIndex.from_product([
            numbers, names.repeat(reps)], names=names)
        tm.assert_index_equal(m.repeat(reps), expected)

    def test_numpy_repeat(self):
        reps = 2
        numbers = [1, 2, 3]
        names = np.array(['foo', 'bar'])

        m = MultiIndex.from_product([
            numbers, names], names=names)
        expected = MultiIndex.from_product([
            numbers, names.repeat(reps)], names=names)
        tm.assert_index_equal(np.repeat(m, reps), expected)

        msg = "the 'axis' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.repeat, m, reps, axis=1)

    def test_set_name_methods(self):
        # so long as these are synonyms, we don't need to test set_names
        self.assertEqual(self.index.rename, self.index.set_names)
        new_names = [name + "SUFFIX" for name in self.index_names]
        ind = self.index.set_names(new_names)
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, new_names)
        with assertRaisesRegexp(ValueError, "^Length"):
            ind.set_names(new_names + new_names)
        new_names2 = [name + "SUFFIX2" for name in new_names]
        res = ind.set_names(new_names2, inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, new_names2)

        # set names for specific level (# GH7792)
        ind = self.index.set_names(new_names[0], level=0)
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, [new_names[0], self.index_names[1]])

        res = ind.set_names(new_names2[0], level=0, inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, [new_names2[0], self.index_names[1]])

        # set names for multiple levels
        ind = self.index.set_names(new_names, level=[0, 1])
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, new_names)

        res = ind.set_names(new_names2, level=[0, 1], inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, new_names2)

    def test_set_levels(self):
        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        levels = self.index.levels
        new_levels = [[lev + 'a' for lev in level] for level in levels]

        def assert_matching(actual, expected, check_dtype=False):
            # avoid specifying internal representation
            # as much as possible
            self.assertEqual(len(actual), len(expected))
            for act, exp in zip(actual, expected):
                act = np.asarray(act)
                exp = np.asarray(exp)
                tm.assert_numpy_array_equal(act, exp, check_dtype=check_dtype)

        # level changing [w/o mutation]
        ind2 = self.index.set_levels(new_levels)
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # level changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, new_levels)

        # level changing specific level [w/o mutation]
        ind2 = self.index.set_levels(new_levels[0], level=0)
        assert_matching(ind2.levels, [new_levels[0], levels[1]])
        assert_matching(self.index.levels, levels)

        ind2 = self.index.set_levels(new_levels[1], level=1)
        assert_matching(ind2.levels, [levels[0], new_levels[1]])
        assert_matching(self.index.levels, levels)

        # level changing multiple levels [w/o mutation]
        ind2 = self.index.set_levels(new_levels, level=[0, 1])
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # level changing specific level [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels[0], level=0, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, [new_levels[0], levels[1]])
        assert_matching(self.index.levels, levels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels[1], level=1, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, [levels[0], new_levels[1]])
        assert_matching(self.index.levels, levels)

        # level changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, level=[0, 1],
                                         inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # illegal level changing should not change levels
        # GH 13754
        original_index = self.index.copy()
        for inplace in [True, False]:
            with assertRaisesRegexp(ValueError, "^On"):
                self.index.set_levels(['c'], level=0, inplace=inplace)
            assert_matching(self.index.levels, original_index.levels,
                            check_dtype=True)

            with assertRaisesRegexp(ValueError, "^On"):
                self.index.set_labels([0, 1, 2, 3, 4, 5], level=0,
                                      inplace=inplace)
            assert_matching(self.index.labels, original_index.labels,
                            check_dtype=True)

            with assertRaisesRegexp(TypeError, "^Levels"):
                self.index.set_levels('c', level=0, inplace=inplace)
            assert_matching(self.index.levels, original_index.levels,
                            check_dtype=True)

            with assertRaisesRegexp(TypeError, "^Labels"):
                self.index.set_labels(1, level=0, inplace=inplace)
            assert_matching(self.index.labels, original_index.labels,
                            check_dtype=True)

    def test_set_labels(self):
        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        labels = self.index.labels
        major_labels, minor_labels = labels
        major_labels = [(x + 1) % 3 for x in major_labels]
        minor_labels = [(x + 1) % 1 for x in minor_labels]
        new_labels = [major_labels, minor_labels]

        def assert_matching(actual, expected):
            # avoid specifying internal representation
            # as much as possible
            self.assertEqual(len(actual), len(expected))
            for act, exp in zip(actual, expected):
                act = np.asarray(act)
                exp = np.asarray(exp, dtype=np.int8)
                tm.assert_numpy_array_equal(act, exp)

        # label changing [w/o mutation]
        ind2 = self.index.set_labels(new_labels)
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, new_labels)

        # label changing specific level [w/o mutation]
        ind2 = self.index.set_labels(new_labels[0], level=0)
        assert_matching(ind2.labels, [new_labels[0], labels[1]])
        assert_matching(self.index.labels, labels)

        ind2 = self.index.set_labels(new_labels[1], level=1)
        assert_matching(ind2.labels, [labels[0], new_labels[1]])
        assert_matching(self.index.labels, labels)

        # label changing multiple levels [w/o mutation]
        ind2 = self.index.set_labels(new_labels, level=[0, 1])
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing specific level [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels[0], level=0, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, [new_labels[0], labels[1]])
        assert_matching(self.index.labels, labels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels[1], level=1, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, [labels[0], new_labels[1]])
        assert_matching(self.index.labels, labels)

        # label changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, level=[0, 1],
                                         inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

    def test_set_levels_labels_names_bad_input(self):
        levels, labels = self.index.levels, self.index.labels
        names = self.index.names

        with tm.assertRaisesRegexp(ValueError, 'Length of levels'):
            self.index.set_levels([levels[0]])

        with tm.assertRaisesRegexp(ValueError, 'Length of labels'):
            self.index.set_labels([labels[0]])

        with tm.assertRaisesRegexp(ValueError, 'Length of names'):
            self.index.set_names([names[0]])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_names(names[0])

        # should have equal lengths
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_levels(levels, level=0)

        # should have equal lengths
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_labels(labels, level=0)

        # should have equal lengths
        with tm.assertRaisesRegexp(ValueError, 'Length of names'):
            self.index.set_names(names[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'string'):
            self.index.set_names(names, level=0)

    def test_set_levels_categorical(self):
        # GH13854
        index = MultiIndex.from_arrays([list("xyzx"), [0, 1, 2, 3]])
        for ordered in [False, True]:
            cidx = CategoricalIndex(list("bac"), ordered=ordered)
            result = index.set_levels(cidx, 0)
            expected = MultiIndex(levels=[cidx, [0, 1, 2, 3]],
                                  labels=index.labels)
            tm.assert_index_equal(result, expected)

            result_lvl = result.get_level_values(0)
            expected_lvl = CategoricalIndex(list("bacb"),
                                            categories=cidx.categories,
                                            ordered=cidx.ordered)
            tm.assert_index_equal(result_lvl, expected_lvl)

    def test_metadata_immutable(self):
        levels, labels = self.index.levels, self.index.labels
        # shouldn't be able to set at either the top level or base level
        mutable_regex = re.compile('does not support mutable operations')
        with assertRaisesRegexp(TypeError, mutable_regex):
            levels[0] = levels[0]
        with assertRaisesRegexp(TypeError, mutable_regex):
            levels[0][0] = levels[0][0]
        # ditto for labels
        with assertRaisesRegexp(TypeError, mutable_regex):
            labels[0] = labels[0]
        with assertRaisesRegexp(TypeError, mutable_regex):
            labels[0][0] = labels[0][0]
        # and for names
        names = self.index.names
        with assertRaisesRegexp(TypeError, mutable_regex):
            names[0] = names[0]

    def test_inplace_mutation_resets_values(self):
        levels = [['a', 'b', 'c'], [4]]
        levels2 = [[1, 2, 3], ['a']]
        labels = [[0, 1, 0, 2, 2, 0], [0, 0, 0, 0, 0, 0]]
        mi1 = MultiIndex(levels=levels, labels=labels)
        mi2 = MultiIndex(levels=levels2, labels=labels)
        vals = mi1.values.copy()
        vals2 = mi2.values.copy()
        self.assertIsNotNone(mi1._tuples)

        # make sure level setting works
        new_vals = mi1.set_levels(levels2).values
        assert_almost_equal(vals2, new_vals)
        # non-inplace doesn't kill _tuples [implementation detail]
        assert_almost_equal(mi1._tuples, vals)
        # and values is still same too
        assert_almost_equal(mi1.values, vals)

        # inplace should kill _tuples
        mi1.set_levels(levels2, inplace=True)
        assert_almost_equal(mi1.values, vals2)

        # make sure label setting works too
        labels2 = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        exp_values = np.empty((6, ), dtype=object)
        exp_values[:] = [(long(1), 'a')] * 6
        # must be 1d array of tuples
        self.assertEqual(exp_values.shape, (6, ))
        new_values = mi2.set_labels(labels2).values
        # not inplace shouldn't change
        assert_almost_equal(mi2._tuples, vals2)
        # should have correct values
        assert_almost_equal(exp_values, new_values)

        # and again setting inplace should kill _tuples, etc
        mi2.set_labels(labels2, inplace=True)
        assert_almost_equal(mi2.values, new_values)

    def test_copy_in_constructor(self):
        levels = np.array(["a", "b", "c"])
        labels = np.array([1, 1, 2, 0, 0, 1, 1])
        val = labels[0]
        mi = MultiIndex(levels=[levels, levels], labels=[labels, labels],
                        copy=True)
        self.assertEqual(mi.labels[0][0], val)
        labels[0] = 15
        self.assertEqual(mi.labels[0][0], val)
        val = levels[0]
        levels[0] = "PANDA"
        self.assertEqual(mi.levels[0][0], val)

    def test_set_value_keeps_names(self):
        # motivating example from #3742
        lev1 = ['hans', 'hans', 'hans', 'grethe', 'grethe', 'grethe']
        lev2 = ['1', '2', '3'] * 2
        idx = pd.MultiIndex.from_arrays([lev1, lev2], names=['Name', 'Number'])
        df = pd.DataFrame(
            np.random.randn(6, 4),
            columns=['one', 'two', 'three', 'four'],
            index=idx)
        df = df.sortlevel()
        self.assertIsNone(df.is_copy)
        self.assertEqual(df.index.names, ('Name', 'Number'))
        df = df.set_value(('grethe', '4'), 'one', 99.34)
        self.assertIsNone(df.is_copy)
        self.assertEqual(df.index.names, ('Name', 'Number'))

    def test_copy_names(self):
        # Check that adding a "names" parameter to the copy is honored
        # GH14302
        multi_idx = pd.Index([(1, 2), (3, 4)], names=['MyName1', 'MyName2'])
        multi_idx1 = multi_idx.copy()

        self.assertTrue(multi_idx.equals(multi_idx1))
        self.assertEqual(multi_idx.names, ['MyName1', 'MyName2'])
        self.assertEqual(multi_idx1.names, ['MyName1', 'MyName2'])

        multi_idx2 = multi_idx.copy(names=['NewName1', 'NewName2'])

        self.assertTrue(multi_idx.equals(multi_idx2))
        self.assertEqual(multi_idx.names, ['MyName1', 'MyName2'])
        self.assertEqual(multi_idx2.names, ['NewName1', 'NewName2'])

        multi_idx3 = multi_idx.copy(name=['NewName1', 'NewName2'])

        self.assertTrue(multi_idx.equals(multi_idx3))
        self.assertEqual(multi_idx.names, ['MyName1', 'MyName2'])
        self.assertEqual(multi_idx3.names, ['NewName1', 'NewName2'])

    def test_names(self):

        # names are assigned in __init__
        names = self.index_names
        level_names = [level.name for level in self.index.levels]
        self.assertEqual(names, level_names)

        # setting bad names on existing
        index = self.index
        assertRaisesRegexp(ValueError, "^Length of names", setattr, index,
                           "names", list(index.names) + ["third"])
        assertRaisesRegexp(ValueError, "^Length of names", setattr, index,
                           "names", [])

        # initializing with bad names (should always be equivalent)
        major_axis, minor_axis = self.index.levels
        major_labels, minor_labels = self.index.labels
        assertRaisesRegexp(ValueError, "^Length of names", MultiIndex,
                           levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels],
                           names=['first'])
        assertRaisesRegexp(ValueError, "^Length of names", MultiIndex,
                           levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels],
                           names=['first', 'second', 'third'])

        # names are assigned
        index.names = ["a", "b"]
        ind_names = list(index.names)
        level_names = [level.name for level in index.levels]
        self.assertEqual(ind_names, level_names)

    def test_reference_duplicate_name(self):
        idx = MultiIndex.from_tuples(
            [('a', 'b'), ('c', 'd')], names=['x', 'x'])
        self.assertTrue(idx._reference_duplicate_name('x'))

        idx = MultiIndex.from_tuples(
            [('a', 'b'), ('c', 'd')], names=['x', 'y'])
        self.assertFalse(idx._reference_duplicate_name('x'))

    def test_astype(self):
        expected = self.index.copy()
        actual = self.index.astype('O')
        assert_copy(actual.levels, expected.levels)
        assert_copy(actual.labels, expected.labels)
        self.check_level_names(actual, expected.names)

        with assertRaisesRegexp(TypeError, "^Setting.*dtype.*object"):
            self.index.astype(np.dtype(int))

    def test_constructor_single_level(self):
        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]], names=['first'])
        tm.assertIsInstance(single_level, Index)
        self.assertNotIsInstance(single_level, MultiIndex)
        self.assertEqual(single_level.name, 'first')

        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]])
        self.assertIsNone(single_level.name)

    def test_constructor_no_levels(self):
        assertRaisesRegexp(ValueError, "non-zero number of levels/labels",
                           MultiIndex, levels=[], labels=[])
        both_re = re.compile('Must pass both levels and labels')
        with tm.assertRaisesRegexp(TypeError, both_re):
            MultiIndex(levels=[])
        with tm.assertRaisesRegexp(TypeError, both_re):
            MultiIndex(labels=[])

    def test_constructor_mismatched_label_levels(self):
        labels = [np.array([1]), np.array([2]), np.array([3])]
        levels = ["a"]
        assertRaisesRegexp(ValueError, "Length of levels and labels must be"
                           " the same", MultiIndex, levels=levels,
                           labels=labels)
        length_error = re.compile('>= length of level')
        label_error = re.compile(r'Unequal label lengths: \[4, 2\]')

        # important to check that it's looking at the right thing.
        with tm.assertRaisesRegexp(ValueError, length_error):
            MultiIndex(levels=[['a'], ['b']],
                       labels=[[0, 1, 2, 3], [0, 3, 4, 1]])

        with tm.assertRaisesRegexp(ValueError, label_error):
            MultiIndex(levels=[['a'], ['b']], labels=[[0, 0, 0, 0], [0, 0]])

        # external API
        with tm.assertRaisesRegexp(ValueError, length_error):
            self.index.copy().set_levels([['a'], ['b']])

        with tm.assertRaisesRegexp(ValueError, label_error):
            self.index.copy().set_labels([[0, 0, 0, 0], [0, 0]])

        # deprecated properties
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            with tm.assertRaisesRegexp(ValueError, length_error):
                self.index.copy().levels = [['a'], ['b']]

            with tm.assertRaisesRegexp(ValueError, label_error):
                self.index.copy().labels = [[0, 0, 0, 0], [0, 0]]

    def assert_multiindex_copied(self, copy, original):
        # levels should be (at least, shallow copied)
        assert_copy(copy.levels, original.levels)

        assert_almost_equal(copy.labels, original.labels)

        # labels doesn't matter which way copied
        assert_almost_equal(copy.labels, original.labels)
        self.assertIsNot(copy.labels, original.labels)

        # names doesn't matter which way copied
        self.assertEqual(copy.names, original.names)
        self.assertIsNot(copy.names, original.names)

        # sort order should be copied
        self.assertEqual(copy.sortorder, original.sortorder)

    def test_copy(self):
        i_copy = self.index.copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_shallow_copy(self):
        i_copy = self.index._shallow_copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_view(self):
        i_view = self.index.view()

        self.assert_multiindex_copied(i_view, self.index)

    def check_level_names(self, index, names):
        self.assertEqual([level.name for level in index.levels], list(names))

    def test_changing_names(self):

        # names should be applied to levels
        level_names = [level.name for level in self.index.levels]
        self.check_level_names(self.index, self.index.names)

        view = self.index.view()
        copy = self.index.copy()
        shallow_copy = self.index._shallow_copy()

        # changing names should change level names on object
        new_names = [name + "a" for name in self.index.names]
        self.index.names = new_names
        self.check_level_names(self.index, new_names)

        # but not on copies
        self.check_level_names(view, level_names)
        self.check_level_names(copy, level_names)
        self.check_level_names(shallow_copy, level_names)

        # and copies shouldn't change original
        shallow_copy.names = [name + "c" for name in shallow_copy.names]
        self.check_level_names(self.index, new_names)

    def test_duplicate_names(self):
        self.index.names = ['foo', 'foo']
        assertRaisesRegexp(KeyError, 'Level foo not found',
                           self.index._get_level_number, 'foo')

    def test_get_level_number_integer(self):
        self.index.names = [1, 0]
        self.assertEqual(self.index._get_level_number(1), 0)
        self.assertEqual(self.index._get_level_number(0), 1)
        self.assertRaises(IndexError, self.index._get_level_number, 2)
        assertRaisesRegexp(KeyError, 'Level fourth not found',
                           self.index._get_level_number, 'fourth')

    def test_from_arrays(self):
        arrays = []
        for lev, lab in zip(self.index.levels, self.index.labels):
            arrays.append(np.asarray(lev).take(lab))

        result = MultiIndex.from_arrays(arrays)
        self.assertEqual(list(result), list(self.index))

        # infer correctly
        result = MultiIndex.from_arrays([[pd.NaT, Timestamp('20130101')],
                                         ['a', 'b']])
        self.assertTrue(result.levels[0].equals(Index([Timestamp('20130101')
                                                       ])))
        self.assertTrue(result.levels[1].equals(Index(['a', 'b'])))

    def test_from_arrays_index_series_datetimetz(self):
        idx1 = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                             tz='US/Eastern')
        idx2 = pd.date_range('2015-01-01 10:00', freq='H', periods=3,
                             tz='Asia/Tokyo')
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_timedelta(self):
        idx1 = pd.timedelta_range('1 days', freq='D', periods=3)
        idx2 = pd.timedelta_range('2 hours', freq='H', periods=3)
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_period(self):
        idx1 = pd.period_range('2011-01-01', freq='D', periods=3)
        idx2 = pd.period_range('2015-01-01', freq='H', periods=3)
        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_datetimelike_mixed(self):
        idx1 = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                             tz='US/Eastern')
        idx2 = pd.date_range('2015-01-01 10:00', freq='H', periods=3)
        idx3 = pd.timedelta_range('1 days', freq='D', periods=3)
        idx4 = pd.period_range('2011-01-01', freq='D', periods=3)

        result = pd.MultiIndex.from_arrays([idx1, idx2, idx3, idx4])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)
        tm.assert_index_equal(result.get_level_values(2), idx3)
        tm.assert_index_equal(result.get_level_values(3), idx4)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1),
                                             pd.Series(idx2),
                                             pd.Series(idx3),
                                             pd.Series(idx4)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)
        tm.assert_index_equal(result2.get_level_values(2), idx3)
        tm.assert_index_equal(result2.get_level_values(3), idx4)

        tm.assert_index_equal(result, result2)

    def test_from_arrays_index_series_categorical(self):
        # GH13743
        idx1 = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                   ordered=False)
        idx2 = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                   ordered=True)

        result = pd.MultiIndex.from_arrays([idx1, idx2])
        tm.assert_index_equal(result.get_level_values(0), idx1)
        tm.assert_index_equal(result.get_level_values(1), idx2)

        result2 = pd.MultiIndex.from_arrays([pd.Series(idx1), pd.Series(idx2)])
        tm.assert_index_equal(result2.get_level_values(0), idx1)
        tm.assert_index_equal(result2.get_level_values(1), idx2)

        result3 = pd.MultiIndex.from_arrays([idx1.values, idx2.values])
        tm.assert_index_equal(result3.get_level_values(0), idx1)
        tm.assert_index_equal(result3.get_level_values(1), idx2)

    def test_from_arrays_empty(self):
        # 0 levels
        with tm.assertRaisesRegexp(
                ValueError, "Must pass non-zero number of levels/labels"):
            MultiIndex.from_arrays(arrays=[])

        # 1 level
        result = MultiIndex.from_arrays(arrays=[[]], names=['A'])
        expected = Index([], name='A')
        tm.assert_index_equal(result, expected)

        # N levels
        for N in [2, 3]:
            arrays = [[]] * N
            names = list('ABC')[:N]
            result = MultiIndex.from_arrays(arrays=arrays, names=names)
            expected = MultiIndex(levels=[np.array([])] * N, labels=[[]] * N,
                                  names=names)
            tm.assert_index_equal(result, expected)

    def test_from_arrays_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        for i in invalid_inputs:
            tm.assertRaises(TypeError, MultiIndex.from_arrays, arrays=i)

    def test_from_arrays_different_lengths(self):
        # GH13599
        idx1 = [1, 2, 3]
        idx2 = ['a', 'b']
        assertRaisesRegexp(ValueError, '^all arrays must be same length$',
                           MultiIndex.from_arrays, [idx1, idx2])

        idx1 = []
        idx2 = ['a', 'b']
        assertRaisesRegexp(ValueError, '^all arrays must be same length$',
                           MultiIndex.from_arrays, [idx1, idx2])

        idx1 = [1, 2, 3]
        idx2 = []
        assertRaisesRegexp(ValueError, '^all arrays must be same length$',
                           MultiIndex.from_arrays, [idx1, idx2])

    def test_from_product(self):

        first = ['foo', 'bar', 'buz']
        second = ['a', 'b', 'c']
        names = ['first', 'second']
        result = MultiIndex.from_product([first, second], names=names)

        tuples = [('foo', 'a'), ('foo', 'b'), ('foo', 'c'), ('bar', 'a'),
                  ('bar', 'b'), ('bar', 'c'), ('buz', 'a'), ('buz', 'b'),
                  ('buz', 'c')]
        expected = MultiIndex.from_tuples(tuples, names=names)

        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, names)

    def test_from_product_empty(self):
        # 0 levels
        with tm.assertRaisesRegexp(
                ValueError, "Must pass non-zero number of levels/labels"):
            MultiIndex.from_product([])

        # 1 level
        result = MultiIndex.from_product([[]], names=['A'])
        expected = pd.Float64Index([], name='A')
        tm.assert_index_equal(result, expected)

        # 2 levels
        l1 = [[], ['foo', 'bar', 'baz'], []]
        l2 = [[], [], ['a', 'b', 'c']]
        names = ['A', 'B']
        for first, second in zip(l1, l2):
            result = MultiIndex.from_product([first, second], names=names)
            expected = MultiIndex(levels=[np.array(first), np.array(second)],
                                  labels=[[], []], names=names)
            tm.assert_index_equal(result, expected)

        # GH12258
        names = ['A', 'B', 'C']
        for N in range(4):
            lvl2 = lrange(N)
            result = MultiIndex.from_product([[], lvl2, []], names=names)
            expected = MultiIndex(levels=[np.array(A)
                                          for A in [[], lvl2, []]],
                                  labels=[[], [], []], names=names)
            tm.assert_index_equal(result, expected)

    def test_from_product_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        for i in invalid_inputs:
            tm.assertRaises(TypeError, MultiIndex.from_product, iterables=i)

    def test_from_product_datetimeindex(self):
        dt_index = date_range('2000-01-01', periods=2)
        mi = pd.MultiIndex.from_product([[1, 2], dt_index])
        etalon = pd.lib.list_to_object_array([(1, pd.Timestamp(
            '2000-01-01')), (1, pd.Timestamp('2000-01-02')), (2, pd.Timestamp(
                '2000-01-01')), (2, pd.Timestamp('2000-01-02'))])
        tm.assert_numpy_array_equal(mi.values, etalon)

    def test_from_product_index_series_categorical(self):
        # GH13743
        first = ['foo', 'bar']
        for ordered in [False, True]:
            idx = pd.CategoricalIndex(list("abcaab"), categories=list("bac"),
                                      ordered=ordered)
            expected = pd.CategoricalIndex(list("abcaab") + list("abcaab"),
                                           categories=list("bac"),
                                           ordered=ordered)

            for arr in [idx, pd.Series(idx), idx.values]:
                result = pd.MultiIndex.from_product([first, arr])
                tm.assert_index_equal(result.get_level_values(1), expected)

    def test_values_boxed(self):
        tuples = [(1, pd.Timestamp('2000-01-01')), (2, pd.NaT),
                  (3, pd.Timestamp('2000-01-03')),
                  (1, pd.Timestamp('2000-01-04')),
                  (2, pd.Timestamp('2000-01-02')),
                  (3, pd.Timestamp('2000-01-03'))]
        mi = pd.MultiIndex.from_tuples(tuples)
        tm.assert_numpy_array_equal(mi.values,
                                    pd.lib.list_to_object_array(tuples))
        # Check that code branches for boxed values produce identical results
        tm.assert_numpy_array_equal(mi.values[:4], mi[:4].values)

    def test_append(self):
        result = self.index[:3].append(self.index[3:])
        self.assertTrue(result.equals(self.index))

        foos = [self.index[:1], self.index[1:3], self.index[3:]]
        result = foos[0].append(foos[1:])
        self.assertTrue(result.equals(self.index))

        # empty
        result = self.index.append([])
        self.assertTrue(result.equals(self.index))

    def test_append_mixed_dtypes(self):
        # GH 13660
        dti = date_range('2011-01-01', freq='M', periods=3,)
        dti_tz = date_range('2011-01-01', freq='M', periods=3, tz='US/Eastern')
        pi = period_range('2011-01', freq='M', periods=3)

        mi = MultiIndex.from_arrays([[1, 2, 3],
                                     [1.1, np.nan, 3.3],
                                     ['a', 'b', 'c'],
                                     dti, dti_tz, pi])
        self.assertEqual(mi.nlevels, 6)

        res = mi.append(mi)
        exp = MultiIndex.from_arrays([[1, 2, 3, 1, 2, 3],
                                     [1.1, np.nan, 3.3, 1.1, np.nan, 3.3],
                                     ['a', 'b', 'c', 'a', 'b', 'c'],
                                     dti.append(dti),
                                     dti_tz.append(dti_tz),
                                     pi.append(pi)])
        tm.assert_index_equal(res, exp)

        other = MultiIndex.from_arrays([['x', 'y', 'z'], ['x', 'y', 'z'],
                                        ['x', 'y', 'z'], ['x', 'y', 'z'],
                                        ['x', 'y', 'z'], ['x', 'y', 'z']])

        res = mi.append(other)
        exp = MultiIndex.from_arrays([[1, 2, 3, 'x', 'y', 'z'],
                                     [1.1, np.nan, 3.3, 'x', 'y', 'z'],
                                     ['a', 'b', 'c', 'x', 'y', 'z'],
                                     dti.append(pd.Index(['x', 'y', 'z'])),
                                     dti_tz.append(pd.Index(['x', 'y', 'z'])),
                                     pi.append(pd.Index(['x', 'y', 'z']))])
        tm.assert_index_equal(res, exp)

    def test_get_level_values(self):
        result = self.index.get_level_values(0)
        expected = Index(['foo', 'foo', 'bar', 'baz', 'qux', 'qux'],
                         name='first')
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, 'first')

        result = self.index.get_level_values('first')
        expected = self.index.get_level_values(0)
        tm.assert_index_equal(result, expected)

        # GH 10460
        index = MultiIndex(levels=[CategoricalIndex(
            ['A', 'B']), CategoricalIndex([1, 2, 3])], labels=[np.array(
                [0, 0, 0, 1, 1, 1]), np.array([0, 1, 2, 0, 1, 2])])
        exp = CategoricalIndex(['A', 'A', 'A', 'B', 'B', 'B'])
        self.assert_index_equal(index.get_level_values(0), exp)
        exp = CategoricalIndex([1, 2, 3, 1, 2, 3])
        self.assert_index_equal(index.get_level_values(1), exp)

    def test_get_level_values_na(self):
        arrays = [['a', 'b', 'b'], [1, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = np.array([1, np.nan, 2])
        tm.assert_numpy_array_equal(values.values.astype(float), expected)

        arrays = [['a', 'b', 'b'], [np.nan, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = np.array([np.nan, np.nan, 2])
        tm.assert_numpy_array_equal(values.values.astype(float), expected)

        arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        expected = np.array([np.nan, np.nan, np.nan])
        tm.assert_numpy_array_equal(values.values.astype(float), expected)
        values = index.get_level_values(1)
        expected = np.array(['a', np.nan, 1], dtype=object)
        tm.assert_numpy_array_equal(values.values, expected)

        arrays = [['a', 'b', 'b'], pd.DatetimeIndex([0, 1, pd.NaT])]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = pd.DatetimeIndex([0, 1, pd.NaT])
        tm.assert_numpy_array_equal(values.values, expected.values)

        arrays = [[], []]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        self.assertEqual(values.shape, (0, ))

    def test_reorder_levels(self):
        # this blows up
        assertRaisesRegexp(IndexError, '^Too many levels',
                           self.index.reorder_levels, [2, 1, 0])

    def test_nlevels(self):
        self.assertEqual(self.index.nlevels, 2)

    def test_iter(self):
        result = list(self.index)
        expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                    ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
        self.assertEqual(result, expected)

    def test_legacy_pickle(self):
        if PY3:
            raise nose.SkipTest("testing for legacy pickles not "
                                "support on py3")

        path = tm.get_data_path('multiindex_v1.pickle')
        obj = pd.read_pickle(path)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assertTrue(obj.equals(obj2))

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj), dtype=np.intp)
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_legacy_v2_unpickle(self):

        # 0.7.3 -> 0.8.0 format manage
        path = tm.get_data_path('mindex_073.pickle')
        obj = pd.read_pickle(path)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assertTrue(obj.equals(obj2))

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj), dtype=np.intp)
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = MultiIndex.from_product(
            [[1, 2], ['a', 'b'], date_range('20130101', periods=3,
                                            tz='US/Eastern')
             ], names=['one', 'two', 'three'])
        unpickled = self.round_trip_pickle(index)
        self.assertTrue(index.equal_levels(unpickled))

    def test_from_tuples_index_values(self):
        result = MultiIndex.from_tuples(self.index)
        self.assertTrue((result.values == self.index.values).all())

    def test_contains(self):
        self.assertIn(('foo', 'two'), self.index)
        self.assertNotIn(('bar', 'two'), self.index)
        self.assertNotIn(None, self.index)

    def test_is_all_dates(self):
        self.assertFalse(self.index.is_all_dates)

    def test_is_numeric(self):
        # MultiIndex is never numeric
        self.assertFalse(self.index.is_numeric())

    def test_getitem(self):
        # scalar
        self.assertEqual(self.index[2], ('bar', 'one'))

        # slice
        result = self.index[2:5]
        expected = self.index[[2, 3, 4]]
        self.assertTrue(result.equals(expected))

        # boolean
        result = self.index[[True, False, True, False, True, True]]
        result2 = self.index[np.array([True, False, True, False, True, True])]
        expected = self.index[[0, 2, 4, 5]]
        self.assertTrue(result.equals(expected))
        self.assertTrue(result2.equals(expected))

    def test_getitem_group_select(self):
        sorted_idx, _ = self.index.sortlevel(0)
        self.assertEqual(sorted_idx.get_loc('baz'), slice(3, 4))
        self.assertEqual(sorted_idx.get_loc('foo'), slice(0, 2))

    def test_get_loc(self):
        self.assertEqual(self.index.get_loc(('foo', 'two')), 1)
        self.assertEqual(self.index.get_loc(('baz', 'two')), 3)
        self.assertRaises(KeyError, self.index.get_loc, ('bar', 'two'))
        self.assertRaises(KeyError, self.index.get_loc, 'quux')

        self.assertRaises(NotImplementedError, self.index.get_loc, 'foo',
                          method='nearest')

        # 3 levels
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])
        self.assertRaises(KeyError, index.get_loc, (1, 1))
        self.assertEqual(index.get_loc((2, 0)), slice(3, 5))

    def test_get_loc_duplicates(self):
        index = Index([2, 2, 2, 2])
        result = index.get_loc(2)
        expected = slice(0, 4)
        self.assertEqual(result, expected)
        # self.assertRaises(Exception, index.get_loc, 2)

        index = Index(['c', 'a', 'a', 'b', 'b'])
        rs = index.get_loc('c')
        xp = 0
        assert (rs == xp)

    def test_get_loc_level(self):
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        loc, new_index = index.get_loc_level((0, 1))
        expected = slice(1, 2)
        exp_index = index[expected].droplevel(0).droplevel(0)
        self.assertEqual(loc, expected)
        self.assertTrue(new_index.equals(exp_index))

        loc, new_index = index.get_loc_level((0, 1, 0))
        expected = 1
        self.assertEqual(loc, expected)
        self.assertIsNone(new_index)

        self.assertRaises(KeyError, index.get_loc_level, (2, 2))

        index = MultiIndex(levels=[[2000], lrange(4)], labels=[np.array(
            [0, 0, 0, 0]), np.array([0, 1, 2, 3])])
        result, new_index = index.get_loc_level((2000, slice(None, None)))
        expected = slice(None, None)
        self.assertEqual(result, expected)
        self.assertTrue(new_index.equals(index.droplevel(0)))

    def test_slice_locs(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index

        slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
        sliced = stacked[slob]
        expected = df[5:16].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

        slob = slice(*idx.slice_locs(df.index[5] + timedelta(seconds=30),
                                     df.index[15] - timedelta(seconds=30)))
        sliced = stacked[slob]
        expected = df[6:15].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

    def test_slice_locs_with_type_mismatch(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index
        assertRaisesRegexp(TypeError, '^Level type mismatch', idx.slice_locs,
                           (1, 3))
        assertRaisesRegexp(TypeError, '^Level type mismatch', idx.slice_locs,
                           df.index[5] + timedelta(seconds=30), (5, 2))
        df = tm.makeCustomDataframe(5, 5)
        stacked = df.stack()
        idx = stacked.index
        with assertRaisesRegexp(TypeError, '^Level type mismatch'):
            idx.slice_locs(timedelta(seconds=30))
        # TODO: Try creating a UnicodeDecodeError in exception message
        with assertRaisesRegexp(TypeError, '^Level type mismatch'):
            idx.slice_locs(df.index[1], (16, "a"))

    def test_slice_locs_not_sorted(self):
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        assertRaisesRegexp(KeyError, "[Kk]ey length.*greater than MultiIndex"
                           " lexsort depth", index.slice_locs, (1, 0, 1),
                           (2, 1, 0))

        # works
        sorted_index, _ = index.sortlevel(0)
        # should there be a test case here???
        sorted_index.slice_locs((1, 0, 1), (2, 1, 0))

    def test_slice_locs_partial(self):
        sorted_idx, _ = self.index.sortlevel(0)

        result = sorted_idx.slice_locs(('foo', 'two'), ('qux', 'one'))
        self.assertEqual(result, (1, 5))

        result = sorted_idx.slice_locs(None, ('qux', 'one'))
        self.assertEqual(result, (0, 5))

        result = sorted_idx.slice_locs(('foo', 'two'), None)
        self.assertEqual(result, (1, len(sorted_idx)))

        result = sorted_idx.slice_locs('bar', 'baz')
        self.assertEqual(result, (2, 4))

    def test_slice_locs_not_contained(self):
        # some searchsorted action

        index = MultiIndex(levels=[[0, 2, 4, 6], [0, 2, 4]],
                           labels=[[0, 0, 0, 1, 1, 2, 3, 3, 3],
                                   [0, 1, 2, 1, 2, 2, 0, 1, 2]], sortorder=0)

        result = index.slice_locs((1, 0), (5, 2))
        self.assertEqual(result, (3, 6))

        result = index.slice_locs(1, 5)
        self.assertEqual(result, (3, 6))

        result = index.slice_locs((2, 2), (5, 2))
        self.assertEqual(result, (3, 6))

        result = index.slice_locs(2, 5)
        self.assertEqual(result, (3, 6))

        result = index.slice_locs((1, 0), (6, 3))
        self.assertEqual(result, (3, 8))

        result = index.slice_locs(-1, 10)
        self.assertEqual(result, (0, len(index)))

    def test_consistency(self):
        # need to construct an overflow
        major_axis = lrange(70000)
        minor_axis = lrange(10)

        major_labels = np.arange(70000)
        minor_labels = np.repeat(lrange(10), 7000)

        # the fact that is works means it's consistent
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        # inconsistent
        major_labels = np.array([0, 0, 1, 1, 1, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1])
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        self.assertFalse(index.is_unique)

    def test_truncate(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        result = index.truncate(before=1)
        self.assertNotIn('foo', result.levels[0])
        self.assertIn(1, result.levels[0])

        result = index.truncate(after=1)
        self.assertNotIn(2, result.levels[0])
        self.assertIn(1, result.levels[0])

        result = index.truncate(before=1, after=2)
        self.assertEqual(len(result.levels[0]), 2)

        # after < before
        self.assertRaises(ValueError, index.truncate, 3, 1)

    def test_get_indexer(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3, 3], dtype=np.intp)
        minor_labels = np.array([0, 1, 0, 0, 1, 0, 1], dtype=np.intp)

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        idx1 = index[:5]
        idx2 = index[[1, 3, 5]]

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, np.array([1, 3, -1], dtype=np.intp))

        r1 = idx2.get_indexer(idx1, method='pad')
        e1 = np.array([-1, 0, 0, 1, 1], dtype=np.intp)
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='pad')
        assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        e1 = np.array([0, 0, 1, 1, 2], dtype=np.intp)
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='backfill')
        assert_almost_equal(r2, e1[::-1])

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        # pass non-MultiIndex
        r1 = idx1.get_indexer(idx2._tuple_index)
        rexp1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1, 2, 3])
        self.assertTrue((r1 == [-1, -1, -1]).all())

        # create index with duplicates
        idx1 = Index(lrange(10) + lrange(10))
        idx2 = Index(lrange(20))

        msg = "Reindexing only valid with uniquely valued Index objects"
        with assertRaisesRegexp(InvalidIndexError, msg):
            idx1.get_indexer(idx2)

    def test_get_indexer_nearest(self):
        midx = MultiIndex.from_tuples([('a', 1), ('b', 2)])
        with tm.assertRaises(NotImplementedError):
            midx.get_indexer(['a'], method='nearest')
        with tm.assertRaises(NotImplementedError):
            midx.get_indexer(['a'], method='pad', tolerance=2)

    def test_format(self):
        self.index.format()
        self.index[:0].format()

    def test_format_integer_names(self):
        index = MultiIndex(levels=[[0, 1], [0, 1]],
                           labels=[[0, 0, 1, 1], [0, 1, 0, 1]], names=[0, 1])
        index.format(names=True)

    def test_format_sparse_display(self):
        index = MultiIndex(levels=[[0, 1], [0, 1], [0, 1], [0]],
                           labels=[[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1],
                                   [0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0]])

        result = index.format()
        self.assertEqual(result[3], '1  0  0  0')

    def test_format_sparse_config(self):
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore', category=FutureWarning,
                                module=".*format")
        # GH1538
        pd.set_option('display.multi_sparse', False)

        result = self.index.format()
        self.assertEqual(result[1], 'foo  two')

        self.reset_display_options()

        warnings.filters = warn_filters

    def test_to_hierarchical(self):
        index = MultiIndex.from_tuples([(1, 'one'), (1, 'two'), (2, 'one'), (
            2, 'two')])
        result = index.to_hierarchical(3)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                      [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

        # K > 1
        result = index.to_hierarchical(3, 2)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                                      [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

        # non-sorted
        index = MultiIndex.from_tuples([(2, 'c'), (1, 'b'),
                                        (2, 'a'), (2, 'b')],
                                       names=['N1', 'N2'])

        result = index.to_hierarchical(2)
        expected = MultiIndex.from_tuples([(2, 'c'), (2, 'c'), (1, 'b'),
                                           (1, 'b'),
                                           (2, 'a'), (2, 'a'),
                                           (2, 'b'), (2, 'b')],
                                          names=['N1', 'N2'])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

    def test_bounds(self):
        self.index._bounds

    def test_equals_multi(self):
        self.assertTrue(self.index.equals(self.index))
        self.assertTrue(self.index.equal_levels(self.index))

        self.assertFalse(self.index.equals(self.index[:-1]))

        self.assertTrue(self.index.equals(self.index._tuple_index))

        # different number of levels
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        index2 = MultiIndex(levels=index.levels[:-1], labels=index.labels[:-1])
        self.assertFalse(index.equals(index2))
        self.assertFalse(index.equal_levels(index2))

        # levels are different
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assertFalse(self.index.equals(index))
        self.assertFalse(self.index.equal_levels(index))

        # some of the labels are different
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assertFalse(self.index.equals(index))

    def test_identical(self):
        mi = self.index.copy()
        mi2 = self.index.copy()
        self.assertTrue(mi.identical(mi2))

        mi = mi.set_names(['new1', 'new2'])
        self.assertTrue(mi.equals(mi2))
        self.assertFalse(mi.identical(mi2))

        mi2 = mi2.set_names(['new1', 'new2'])
        self.assertTrue(mi.identical(mi2))

        mi3 = Index(mi.tolist(), names=mi.names)
        mi4 = Index(mi.tolist(), names=mi.names, tupleize_cols=False)
        self.assertTrue(mi.identical(mi3))
        self.assertFalse(mi.identical(mi4))
        self.assertTrue(mi.equals(mi4))

    def test_is_(self):
        mi = MultiIndex.from_tuples(lzip(range(10), range(10)))
        self.assertTrue(mi.is_(mi))
        self.assertTrue(mi.is_(mi.view()))
        self.assertTrue(mi.is_(mi.view().view().view().view()))
        mi2 = mi.view()
        # names are metadata, they don't change id
        mi2.names = ["A", "B"]
        self.assertTrue(mi2.is_(mi))
        self.assertTrue(mi.is_(mi2))

        self.assertTrue(mi.is_(mi.set_names(["C", "D"])))
        mi2 = mi.view()
        mi2.set_names(["E", "F"], inplace=True)
        self.assertTrue(mi.is_(mi2))
        # levels are inherent properties, they change identity
        mi3 = mi2.set_levels([lrange(10), lrange(10)])
        self.assertFalse(mi3.is_(mi2))
        # shouldn't change
        self.assertTrue(mi2.is_(mi))
        mi4 = mi3.view()
        mi4.set_levels([[1 for _ in range(10)], lrange(10)], inplace=True)
        self.assertFalse(mi4.is_(mi3))
        mi5 = mi.view()
        mi5.set_levels(mi5.levels, inplace=True)
        self.assertFalse(mi5.is_(mi))

    def test_union(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_union = piece1 | piece2

        tups = sorted(self.index._tuple_index)
        expected = MultiIndex.from_tuples(tups)

        self.assertTrue(the_union.equals(expected))

        # corner case, pass self or empty thing:
        the_union = self.index.union(self.index)
        self.assertIs(the_union, self.index)

        the_union = self.index.union(self.index[:0])
        self.assertIs(the_union, self.index)

        # won't work in python 3
        # tuples = self.index._tuple_index
        # result = self.index[:4] | tuples[4:]
        # self.assertTrue(result.equals(tuples))

        # not valid for python 3
        # def test_union_with_regular_index(self):
        #     other = Index(['A', 'B', 'C'])

        #     result = other.union(self.index)
        #     self.assertIn(('foo', 'one'), result)
        #     self.assertIn('B', result)

        #     result2 = self.index.union(other)
        #     self.assertTrue(result.equals(result2))

    def test_intersection(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_int = piece1 & piece2
        tups = sorted(self.index[3:5]._tuple_index)
        expected = MultiIndex.from_tuples(tups)
        self.assertTrue(the_int.equals(expected))

        # corner case, pass self
        the_int = self.index.intersection(self.index)
        self.assertIs(the_int, self.index)

        # empty intersection: disjoint
        empty = self.index[:2] & self.index[2:]
        expected = self.index[:0]
        self.assertTrue(empty.equals(expected))

        # can't do in python 3
        # tuples = self.index._tuple_index
        # result = self.index & tuples
        # self.assertTrue(result.equals(tuples))

    def test_sub(self):

        first = self.index

        # - now raises (previously was set op difference)
        with tm.assertRaises(TypeError):
            first - self.index[-3:]
        with tm.assertRaises(TypeError):
            self.index[-3:] - first
        with tm.assertRaises(TypeError):
            self.index[-3:] - first.tolist()
        with tm.assertRaises(TypeError):
            first.tolist() - self.index[-3:]

    def test_difference(self):

        first = self.index
        result = first.difference(self.index[-3:])
        expected = MultiIndex.from_tuples(sorted(self.index[:-3].values),
                                          sortorder=0,
                                          names=self.index.names)

        tm.assertIsInstance(result, MultiIndex)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: reflexive
        result = self.index.difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: superset
        result = self.index[-3:].difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: degenerate
        result = self.index[:0].difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # names not the same
        chunklet = self.index[-3:]
        chunklet.names = ['foo', 'baz']
        result = first.difference(chunklet)
        self.assertEqual(result.names, (None, None))

        # empty, but non-equal
        result = self.index.difference(self.index.sortlevel(1)[0])
        self.assertEqual(len(result), 0)

        # raise Exception called with non-MultiIndex
        result = first.difference(first._tuple_index)
        self.assertTrue(result.equals(first[:0]))

        # name from empty array
        result = first.difference([])
        self.assertTrue(first.equals(result))
        self.assertEqual(first.names, result.names)

        # name from non-empty array
        result = first.difference([('foo', 'one')])
        expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), (
            'foo', 'two'), ('qux', 'one'), ('qux', 'two')])
        expected.names = first.names
        self.assertEqual(first.names, result.names)
        assertRaisesRegexp(TypeError, "other must be a MultiIndex or a list"
                           " of tuples", first.difference, [1, 2, 3, 4, 5])

    def test_from_tuples(self):
        assertRaisesRegexp(TypeError, 'Cannot infer number of levels from'
                           ' empty list', MultiIndex.from_tuples, [])

        idx = MultiIndex.from_tuples(((1, 2), (3, 4)), names=['a', 'b'])
        self.assertEqual(len(idx), 2)

    def test_argsort(self):
        result = self.index.argsort()
        expected = self.index._tuple_index.argsort()
        tm.assert_numpy_array_equal(result, expected)

    def test_sortlevel(self):
        import random

        tuples = list(self.index)
        random.shuffle(tuples)

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

    def test_sortlevel_not_sort_remaining(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        sorted_idx, _ = mi.sortlevel('A', sort_remaining=False)
        self.assertTrue(sorted_idx.equals(mi))

    def test_sortlevel_deterministic(self):
        tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
                  ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

    def test_dims(self):
        pass

    def test_drop(self):
        dropped = self.index.drop([('foo', 'two'), ('qux', 'one')])

        index = MultiIndex.from_tuples([('foo', 'two'), ('qux', 'one')])
        dropped2 = self.index.drop(index)

        expected = self.index[[0, 2, 3, 5]]
        self.assert_index_equal(dropped, expected)
        self.assert_index_equal(dropped2, expected)

        dropped = self.index.drop(['bar'])
        expected = self.index[[0, 1, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop('foo')
        expected = self.index[[2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        index = MultiIndex.from_tuples([('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, [('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, index)
        self.assertRaises(KeyError, self.index.drop, ['foo', 'two'])

        # partially correct argument
        mixed_index = MultiIndex.from_tuples([('qux', 'one'), ('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, mixed_index)

        # error='ignore'
        dropped = self.index.drop(index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop(['foo', 'two'], errors='ignore')
        expected = self.index[[2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        # mixed partial / full drop
        dropped = self.index.drop(['foo', ('qux', 'one')])
        expected = self.index[[2, 3, 5]]
        self.assert_index_equal(dropped, expected)

        # mixed partial / full drop / error='ignore'
        mixed_index = ['foo', ('qux', 'one'), 'two']
        self.assertRaises(KeyError, self.index.drop, mixed_index)
        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[2, 3, 5]]
        self.assert_index_equal(dropped, expected)

    def test_droplevel_with_names(self):
        index = self.index[self.index.get_loc('foo')]
        dropped = index.droplevel(0)
        self.assertEqual(dropped.name, 'second')

        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])],
            names=['one', 'two', 'three'])
        dropped = index.droplevel(0)
        self.assertEqual(dropped.names, ('two', 'three'))

        dropped = index.droplevel('two')
        expected = index.droplevel(1)
        self.assertTrue(dropped.equals(expected))

    def test_droplevel_multiple(self):
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])],
            names=['one', 'two', 'three'])

        dropped = index[:2].droplevel(['three', 'one'])
        expected = index[:2].droplevel(2).droplevel(0)
        self.assertTrue(dropped.equals(expected))

    def test_drop_not_lexsorted(self):
        # GH 12078

        # define the lexsorted version of the multi-index
        tuples = [('a', ''), ('b1', 'c1'), ('b2', 'c2')]
        lexsorted_mi = MultiIndex.from_tuples(tuples, names=['b', 'c'])
        self.assertTrue(lexsorted_mi.is_lexsorted())

        # and the not-lexsorted version
        df = pd.DataFrame(columns=['a', 'b', 'c', 'd'],
                          data=[[1, 'b1', 'c1', 3], [1, 'b2', 'c2', 4]])
        df = df.pivot_table(index='a', columns=['b', 'c'], values='d')
        df = df.reset_index()
        not_lexsorted_mi = df.columns
        self.assertFalse(not_lexsorted_mi.is_lexsorted())

        # compare the results
        self.assert_index_equal(lexsorted_mi, not_lexsorted_mi)
        with self.assert_produces_warning(PerformanceWarning):
            self.assert_index_equal(lexsorted_mi.drop('a'),
                                    not_lexsorted_mi.drop('a'))

    def test_insert(self):
        # key contained in all levels
        new_index = self.index.insert(0, ('bar', 'two'))
        self.assertTrue(new_index.equal_levels(self.index))
        self.assertEqual(new_index[0], ('bar', 'two'))

        # key not contained in all levels
        new_index = self.index.insert(0, ('abc', 'three'))

        exp0 = Index(list(self.index.levels[0]) + ['abc'], name='first')
        tm.assert_index_equal(new_index.levels[0], exp0)

        exp1 = Index(list(self.index.levels[1]) + ['three'], name='second')
        tm.assert_index_equal(new_index.levels[1], exp1)
        self.assertEqual(new_index[0], ('abc', 'three'))

        # key wrong length
        msg = "Item must have length equal to number of levels"
        with assertRaisesRegexp(ValueError, msg):
            self.index.insert(0, ('foo2', ))

        left = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1]],
                            columns=['1st', '2nd', '3rd'])
        left.set_index(['1st', '2nd'], inplace=True)
        ts = left['3rd'].copy(deep=True)

        left.loc[('b', 'x'), '3rd'] = 2
        left.loc[('b', 'a'), '3rd'] = -1
        left.loc[('b', 'b'), '3rd'] = 3
        left.loc[('a', 'x'), '3rd'] = 4
        left.loc[('a', 'w'), '3rd'] = 5
        left.loc[('a', 'a'), '3rd'] = 6

        ts.loc[('b', 'x')] = 2
        ts.loc['b', 'a'] = -1
        ts.loc[('b', 'b')] = 3
        ts.loc['a', 'x'] = 4
        ts.loc[('a', 'w')] = 5
        ts.loc['a', 'a'] = 6

        right = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1], ['b', 'x', 2],
                              ['b', 'a', -1], ['b', 'b', 3], ['a', 'x', 4],
                              ['a', 'w', 5], ['a', 'a', 6]],
                             columns=['1st', '2nd', '3rd'])
        right.set_index(['1st', '2nd'], inplace=True)
        # FIXME data types changes to float because
        # of intermediate nan insertion;
        tm.assert_frame_equal(left, right, check_dtype=False)
        tm.assert_series_equal(ts, right['3rd'])

        # GH9250
        idx = [('test1', i) for i in range(5)] + \
            [('test2', i) for i in range(6)] + \
            [('test', 17), ('test', 18)]

        left = pd.Series(np.linspace(0, 10, 11),
                         pd.MultiIndex.from_tuples(idx[:-2]))

        left.loc[('test', 17)] = 11
        left.ix[('test', 18)] = 12

        right = pd.Series(np.linspace(0, 12, 13),
                          pd.MultiIndex.from_tuples(idx))

        tm.assert_series_equal(left, right)

    def test_take_preserve_name(self):
        taken = self.index.take([3, 0, 1])
        self.assertEqual(taken.names, self.index.names)

    def test_take_fill_value(self):
        # GH 12631
        vals = [['A', 'B'],
                [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
        idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])

        result = idx.take(np.array([1, 0, -1]))
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    ('B', pd.Timestamp('2011-01-02'))]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    (np.nan, pd.NaT)]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    ('B', pd.Timestamp('2011-01-02'))]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def take_invalid_kwargs(self):
        vals = [['A', 'B'],
                [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
        idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])
        indices = [1, 2]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        tm.assertRaisesRegexp(TypeError, msg, idx.take,
                              indices, foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, mode='clip')

    def test_join_level(self):
        def _check_how(other, how):
            join_index, lidx, ridx = other.join(self.index, how=how,
                                                level='second',
                                                return_indexers=True)

            exp_level = other.join(self.index.levels[1], how=how)
            self.assertTrue(join_index.levels[0].equals(self.index.levels[0]))
            self.assertTrue(join_index.levels[1].equals(exp_level))

            # pare down levels
            mask = np.array(
                [x[1] in exp_level for x in self.index], dtype=bool)
            exp_values = self.index.values[mask]
            tm.assert_numpy_array_equal(join_index.values, exp_values)

            if how in ('outer', 'inner'):
                join_index2, ridx2, lidx2 = \
                    self.index.join(other, how=how, level='second',
                                    return_indexers=True)

                self.assertTrue(join_index.equals(join_index2))
                tm.assert_numpy_array_equal(lidx, lidx2)
                tm.assert_numpy_array_equal(ridx, ridx2)
                tm.assert_numpy_array_equal(join_index2.values, exp_values)

        def _check_all(other):
            _check_how(other, 'outer')
            _check_how(other, 'inner')
            _check_how(other, 'left')
            _check_how(other, 'right')

        _check_all(Index(['three', 'one', 'two']))
        _check_all(Index(['one']))
        _check_all(Index(['one', 'three']))

        # some corner cases
        idx = Index(['three', 'one', 'two'])
        result = idx.join(self.index, level='second')
        tm.assertIsInstance(result, MultiIndex)

        assertRaisesRegexp(TypeError, "Join.*MultiIndex.*ambiguous",
                           self.index.join, self.index, level=1)

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            res = self.index
            joined = res.join(res, how=kind)
            self.assertIs(res, joined)

    def test_join_multi(self):
        # GH 10665
        midx = pd.MultiIndex.from_product(
            [np.arange(4), np.arange(4)], names=['a', 'b'])
        idx = pd.Index([1, 2, 5], name='b')

        # inner
        jidx, lidx, ridx = midx.join(idx, how='inner', return_indexers=True)
        exp_idx = pd.MultiIndex.from_product(
            [np.arange(4), [1, 2]], names=['a', 'b'])
        exp_lidx = np.array([1, 2, 5, 6, 9, 10, 13, 14], dtype=np.intp)
        exp_ridx = np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=np.intp)
        self.assert_index_equal(jidx, exp_idx)
        self.assert_numpy_array_equal(lidx, exp_lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='inner', return_indexers=True)
        self.assert_index_equal(jidx, exp_idx)
        self.assert_numpy_array_equal(lidx, exp_lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)

        # keep MultiIndex
        jidx, lidx, ridx = midx.join(idx, how='left', return_indexers=True)
        exp_ridx = np.array([-1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1, -1, 0,
                             1, -1], dtype=np.intp)
        self.assert_index_equal(jidx, midx)
        self.assertIsNone(lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='right', return_indexers=True)
        self.assert_index_equal(jidx, midx)
        self.assertIsNone(lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)

    def test_reindex(self):
        result, indexer = self.index.reindex(list(self.index[:4]))
        tm.assertIsInstance(result, MultiIndex)
        self.check_level_names(result, self.index[:4].names)

        result, indexer = self.index.reindex(list(self.index))
        tm.assertIsInstance(result, MultiIndex)
        self.assertIsNone(indexer)
        self.check_level_names(result, self.index.names)

    def test_reindex_level(self):
        idx = Index(['one'])

        target, indexer = self.index.reindex(idx, level='second')
        target2, indexer2 = idx.reindex(self.index, level='second')

        exp_index = self.index.join(idx, level='second', how='right')
        exp_index2 = self.index.join(idx, level='second', how='left')

        self.assertTrue(target.equals(exp_index))
        exp_indexer = np.array([0, 2, 4])
        tm.assert_numpy_array_equal(indexer, exp_indexer, check_dtype=False)

        self.assertTrue(target2.equals(exp_index2))
        exp_indexer2 = np.array([0, -1, 0, -1, 0, -1])
        tm.assert_numpy_array_equal(indexer2, exp_indexer2, check_dtype=False)

        assertRaisesRegexp(TypeError, "Fill method not supported",
                           self.index.reindex, self.index, method='pad',
                           level='second')

        assertRaisesRegexp(TypeError, "Fill method not supported", idx.reindex,
                           idx, method='bfill', level='first')

    def test_duplicates(self):
        self.assertFalse(self.index.has_duplicates)
        self.assertTrue(self.index.append(self.index).has_duplicates)

        index = MultiIndex(levels=[[0, 1], [0, 1, 2]], labels=[
                           [0, 0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 0, 1, 2]])
        self.assertTrue(index.has_duplicates)

        # GH 9075
        t = [(u('x'), u('out'), u('z'), 5, u('y'), u('in'), u('z'), 169),
             (u('x'), u('out'), u('z'), 7, u('y'), u('in'), u('z'), 119),
             (u('x'), u('out'), u('z'), 9, u('y'), u('in'), u('z'), 135),
             (u('x'), u('out'), u('z'), 13, u('y'), u('in'), u('z'), 145),
             (u('x'), u('out'), u('z'), 14, u('y'), u('in'), u('z'), 158),
             (u('x'), u('out'), u('z'), 16, u('y'), u('in'), u('z'), 122),
             (u('x'), u('out'), u('z'), 17, u('y'), u('in'), u('z'), 160),
             (u('x'), u('out'), u('z'), 18, u('y'), u('in'), u('z'), 180),
             (u('x'), u('out'), u('z'), 20, u('y'), u('in'), u('z'), 143),
             (u('x'), u('out'), u('z'), 21, u('y'), u('in'), u('z'), 128),
             (u('x'), u('out'), u('z'), 22, u('y'), u('in'), u('z'), 129),
             (u('x'), u('out'), u('z'), 25, u('y'), u('in'), u('z'), 111),
             (u('x'), u('out'), u('z'), 28, u('y'), u('in'), u('z'), 114),
             (u('x'), u('out'), u('z'), 29, u('y'), u('in'), u('z'), 121),
             (u('x'), u('out'), u('z'), 31, u('y'), u('in'), u('z'), 126),
             (u('x'), u('out'), u('z'), 32, u('y'), u('in'), u('z'), 155),
             (u('x'), u('out'), u('z'), 33, u('y'), u('in'), u('z'), 123),
             (u('x'), u('out'), u('z'), 12, u('y'), u('in'), u('z'), 144)]

        index = pd.MultiIndex.from_tuples(t)
        self.assertFalse(index.has_duplicates)

        # handle int64 overflow if possible
        def check(nlevels, with_nulls):
            labels = np.tile(np.arange(500), 2)
            level = np.arange(500)

            if with_nulls:  # inject some null values
                labels[500] = -1  # common nan value
                labels = list(labels.copy() for i in range(nlevels))
                for i in range(nlevels):
                    labels[i][500 + i - nlevels // 2] = -1

                labels += [np.array([-1, 1]).repeat(500)]
            else:
                labels = [labels] * nlevels + [np.arange(2).repeat(500)]

            levels = [level] * nlevels + [[0, 1]]

            # no dups
            index = MultiIndex(levels=levels, labels=labels)
            self.assertFalse(index.has_duplicates)

            # with a dup
            if with_nulls:
                f = lambda a: np.insert(a, 1000, a[0])
                labels = list(map(f, labels))
                index = MultiIndex(levels=levels, labels=labels)
            else:
                values = index.values.tolist()
                index = MultiIndex.from_tuples(values + [values[0]])

            self.assertTrue(index.has_duplicates)

        # no overflow
        check(4, False)
        check(4, True)

        # overflow possible
        check(8, False)
        check(8, True)

        # GH 9125
        n, k = 200, 5000
        levels = [np.arange(n), tm.makeStringIndex(n), 1000 + np.arange(n)]
        labels = [np.random.choice(n, k * n) for lev in levels]
        mi = MultiIndex(levels=levels, labels=labels)

        for keep in ['first', 'last', False]:
            left = mi.duplicated(keep=keep)
            right = pd.hashtable.duplicated_object(mi.values, keep=keep)
            tm.assert_numpy_array_equal(left, right)

        # GH5873
        for a in [101, 102]:
            mi = MultiIndex.from_arrays([[101, a], [3.5, np.nan]])
            self.assertFalse(mi.has_duplicates)
            self.assertEqual(mi.get_duplicates(), [])
            tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
                2, dtype='bool'))

        for n in range(1, 6):  # 1st level shape
            for m in range(1, 5):  # 2nd level shape
                # all possible unique combinations, including nan
                lab = product(range(-1, n), range(-1, m))
                mi = MultiIndex(levels=[list('abcde')[:n], list('WXYZ')[:m]],
                                labels=np.random.permutation(list(lab)).T)
                self.assertEqual(len(mi), (n + 1) * (m + 1))
                self.assertFalse(mi.has_duplicates)
                self.assertEqual(mi.get_duplicates(), [])
                tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
                    len(mi), dtype='bool'))

    def test_duplicate_meta_data(self):
        # GH 10115
        index = MultiIndex(levels=[[0, 1], [0, 1, 2]], labels=[
                           [0, 0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 0, 1, 2]])
        for idx in [index,
                    index.set_names([None, None]),
                    index.set_names([None, 'Num']),
                    index.set_names(['Upper', 'Num']), ]:
            self.assertTrue(idx.has_duplicates)
            self.assertEqual(idx.drop_duplicates().names, idx.names)

    def test_get_unique_index(self):
        idx = self.index[[0, 1, 0, 1, 1, 0, 0]]
        expected = self.index._shallow_copy(idx[[0, 1]])

        for dropna in [False, True]:
            result = idx._get_unique_index(dropna=dropna)
            self.assertTrue(result.unique)
            self.assert_index_equal(result, expected)

    def test_unique(self):
        mi = pd.MultiIndex.from_arrays([[1, 2, 1, 2], [1, 1, 1, 2]])

        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([[1, 2, 2], [1, 1, 2]])
        tm.assert_index_equal(res, exp)

        mi = pd.MultiIndex.from_arrays([list('aaaa'), list('abab')])
        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([list('aa'), list('ab')])
        tm.assert_index_equal(res, exp)

        mi = pd.MultiIndex.from_arrays([list('aaaa'), list('aaaa')])
        res = mi.unique()
        exp = pd.MultiIndex.from_arrays([['a'], ['a']])
        tm.assert_index_equal(res, exp)

    def test_unique_datetimelike(self):
        idx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-01',
                                 '2015-01-01', 'NaT', 'NaT'])
        idx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-02',
                                 '2015-01-02', 'NaT', '2015-01-01'],
                                tz='Asia/Tokyo')
        result = pd.MultiIndex.from_arrays([idx1, idx2]).unique()

        eidx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', 'NaT', 'NaT'])
        eidx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-02',
                                  'NaT', '2015-01-01'],
                                 tz='Asia/Tokyo')
        exp = pd.MultiIndex.from_arrays([eidx1, eidx2])
        tm.assert_index_equal(result, exp)

    def test_tolist(self):
        result = self.index.tolist()
        exp = list(self.index.values)
        self.assertEqual(result, exp)

    def test_repr_with_unicode_data(self):
        with pd.core.config.option_context("display.encoding", 'UTF-8'):
            d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
            index = pd.DataFrame(d).set_index(["a", "b"]).index
            self.assertFalse("\\u" in repr(index)
                             )  # we don't want unicode-escaped

    def test_repr_roundtrip(self):

        mi = MultiIndex.from_product([list('ab'), range(3)],
                                     names=['first', 'second'])
        str(mi)

        if PY3:
            tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
        else:
            result = eval(repr(mi))
            # string coerces to unicode
            tm.assert_index_equal(result, mi, exact=False)
            self.assertEqual(
                mi.get_level_values('first').inferred_type, 'string')
            self.assertEqual(
                result.get_level_values('first').inferred_type, 'unicode')

        mi_u = MultiIndex.from_product(
            [list(u'ab'), range(3)], names=['first', 'second'])
        result = eval(repr(mi_u))
        tm.assert_index_equal(result, mi_u, exact=True)

        # formatting
        if PY3:
            str(mi)
        else:
            compat.text_type(mi)

        # long format
        mi = MultiIndex.from_product([list('abcdefg'), range(10)],
                                     names=['first', 'second'])
        result = str(mi)

        if PY3:
            tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
        else:
            result = eval(repr(mi))
            # string coerces to unicode
            tm.assert_index_equal(result, mi, exact=False)
            self.assertEqual(
                mi.get_level_values('first').inferred_type, 'string')
            self.assertEqual(
                result.get_level_values('first').inferred_type, 'unicode')

        mi = MultiIndex.from_product(
            [list(u'abcdefg'), range(10)], names=['first', 'second'])
        result = eval(repr(mi_u))
        tm.assert_index_equal(result, mi_u, exact=True)

    def test_str(self):
        # tested elsewhere
        pass

    def test_unicode_string_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if PY3:
            bytes(idx)
        else:
            str(idx)

    def test_slice_keep_name(self):
        x = MultiIndex.from_tuples([('a', 'b'), (1, 2), ('c', 'd')],
                                   names=['x', 'y'])
        self.assertEqual(x[1:].names, x.names)

    def test_isnull_behavior(self):
        # should not segfault GH5123
        # NOTE: if MI representation changes, may make sense to allow
        # isnull(MI)
        with tm.assertRaises(NotImplementedError):
            pd.isnull(self.index)

    def test_level_setting_resets_attributes(self):
        ind = MultiIndex.from_arrays([
            ['A', 'A', 'B', 'B', 'B'], [1, 2, 1, 2, 3]
        ])
        assert ind.is_monotonic
        ind.set_levels([['A', 'B', 'A', 'A', 'B'], [2, 1, 3, -2, 5]],
                       inplace=True)
        # if this fails, probably didn't reset the cache correctly.
        assert not ind.is_monotonic

    def test_isin(self):
        values = [('foo', 2), ('bar', 3), ('quux', 4)]

        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
            4)])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = MultiIndex.from_arrays([[], []])
        result = idx.isin(values)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.dtype, np.bool_)

    def test_isin_nan(self):
        idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
        tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                    np.array([False, False]))
        tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                    np.array([False, False]))

    def test_isin_level_kwarg(self):
        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
            4)])

        vals_0 = ['foo', 'bar', 'quux']
        vals_1 = [2, 3, 10]

        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=0))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=-2))

        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=1))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=-1))

        self.assertRaises(IndexError, idx.isin, vals_0, level=5)
        self.assertRaises(IndexError, idx.isin, vals_0, level=-5)

        self.assertRaises(KeyError, idx.isin, vals_0, level=1.0)
        self.assertRaises(KeyError, idx.isin, vals_1, level=-1.0)
        self.assertRaises(KeyError, idx.isin, vals_1, level='A')

        idx.names = ['A', 'B']
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level='A'))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level='B'))

        self.assertRaises(KeyError, idx.isin, vals_1, level='C')

    def test_reindex_preserves_names_when_target_is_list_or_ndarray(self):
        # GH6552
        idx = self.index.copy()
        target = idx.copy()
        idx.names = target.names = [None, None]

        other_dtype = pd.MultiIndex.from_product([[1, 2], [3, 4]])

        # list & ndarray cases
        self.assertEqual(idx.reindex([])[0].names, [None, None])
        self.assertEqual(idx.reindex(np.array([]))[0].names, [None, None])
        self.assertEqual(idx.reindex(target.tolist())[0].names, [None, None])
        self.assertEqual(idx.reindex(target.values)[0].names, [None, None])
        self.assertEqual(
            idx.reindex(other_dtype.tolist())[0].names, [None, None])
        self.assertEqual(
            idx.reindex(other_dtype.values)[0].names, [None, None])

        idx.names = ['foo', 'bar']
        self.assertEqual(idx.reindex([])[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(np.array([]))[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(target.tolist())[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(target.values)[0].names, ['foo', 'bar'])
        self.assertEqual(
            idx.reindex(other_dtype.tolist())[0].names, ['foo', 'bar'])
        self.assertEqual(
            idx.reindex(other_dtype.values)[0].names, ['foo', 'bar'])

    def test_reindex_lvl_preserves_names_when_target_is_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']],
                                         names=['foo', 'bar'])
        self.assertEqual(idx.reindex([], level=0)[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex([], level=1)[0].names, ['foo', 'bar'])

    def test_reindex_lvl_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
        self.assertEqual(idx.reindex([], level=0)[0].levels[0].dtype.type,
                         np.int64)
        self.assertEqual(idx.reindex([], level=1)[0].levels[1].dtype.type,
                         np.object_)

    def test_groupby(self):
        groups = self.index.groupby(np.array([1, 1, 1, 2, 2, 2]))
        labels = self.index.get_values().tolist()
        exp = {1: labels[:3], 2: labels[3:]}
        tm.assert_dict_equal(groups, exp)

        # GH5620
        groups = self.index.groupby(self.index)
        exp = dict((key, [key]) for key in self.index)
        tm.assert_dict_equal(groups, exp)

    def test_index_name_retained(self):
        # GH9857
        result = pd.DataFrame({'x': [1, 2, 6],
                               'y': [2, 2, 8],
                               'z': [-5, 0, 5]})
        result = result.set_index('z')
        result.loc[10] = [9, 10]
        df_expected = pd.DataFrame({'x': [1, 2, 6, 9],
                                    'y': [2, 2, 8, 10],
                                    'z': [-5, 0, 5, 10]})
        df_expected = df_expected.set_index('z')
        tm.assert_frame_equal(result, df_expected)

    def test_equals_operator(self):
        # GH9785
        self.assertTrue((self.index == self.index).all())

    def test_large_multiindex_error(self):
        # GH12527
        df_below_1000000 = pd.DataFrame(
            1, index=pd.MultiIndex.from_product([[1, 2], range(499999)]),
            columns=['dest'])
        with assertRaises(KeyError):
            df_below_1000000.loc[(-1, 0), 'dest']
        with assertRaises(KeyError):
            df_below_1000000.loc[(3, 0), 'dest']
        df_above_1000000 = pd.DataFrame(
            1, index=pd.MultiIndex.from_product([[1, 2], range(500001)]),
            columns=['dest'])
        with assertRaises(KeyError):
            df_above_1000000.loc[(-1, 0), 'dest']
        with assertRaises(KeyError):
            df_above_1000000.loc[(3, 0), 'dest']

    def test_partial_string_timestamp_multiindex(self):
        # GH10331
        dr = pd.date_range('2016-01-01', '2016-01-03', freq='12H')
        abc = ['a', 'b', 'c']
        ix = pd.MultiIndex.from_product([dr, abc])
        df = pd.DataFrame({'c1': range(0, 15)}, index=ix)
        idx = pd.IndexSlice

        #                        c1
        # 2016-01-01 00:00:00 a   0
        #                     b   1
        #                     c   2
        # 2016-01-01 12:00:00 a   3
        #                     b   4
        #                     c   5
        # 2016-01-02 00:00:00 a   6
        #                     b   7
        #                     c   8
        # 2016-01-02 12:00:00 a   9
        #                     b  10
        #                     c  11
        # 2016-01-03 00:00:00 a  12
        #                     b  13
        #                     c  14

        # partial string matching on a single index
        for df_swap in (df.swaplevel(),
                        df.swaplevel(0),
                        df.swaplevel(0, 1)):
            df_swap = df_swap.sort_index()
            just_a = df_swap.loc['a']
            result = just_a.loc['2016-01-01']
            expected = df.loc[idx[:, 'a'], :].iloc[0:2]
            expected.index = expected.index.droplevel(1)
            tm.assert_frame_equal(result, expected)

        # indexing with IndexSlice
        result = df.loc[idx['2016-01-01':'2016-02-01', :], :]
        expected = df
        tm.assert_frame_equal(result, expected)

        # match on secondary index
        result = df_swap.loc[idx[:, '2016-01-01':'2016-01-01'], :]
        expected = df_swap.iloc[[0, 1, 5, 6, 10, 11]]
        tm.assert_frame_equal(result, expected)

        # Even though this syntax works on a single index, this is somewhat
        # ambiguous and we don't want to extend this behavior forward to work
        # in multi-indexes. This would amount to selecting a scalar from a
        # column.
        with assertRaises(KeyError):
            df['2016-01-01']

        # partial string match on year only
        result = df.loc['2016']
        expected = df
        tm.assert_frame_equal(result, expected)

        # partial string match on date
        result = df.loc['2016-01-01']
        expected = df.iloc[0:6]
        tm.assert_frame_equal(result, expected)

        # partial string match on date and hour, from middle
        result = df.loc['2016-01-02 12']
        expected = df.iloc[9:12]
        tm.assert_frame_equal(result, expected)

        # partial string match on secondary index
        result = df_swap.loc[idx[:, '2016-01-02'], :]
        expected = df_swap.iloc[[2, 3, 7, 8, 12, 13]]
        tm.assert_frame_equal(result, expected)

        # tuple selector with partial string match on date
        result = df.loc[('2016-01-01', 'a'), :]
        expected = df.iloc[[0, 3]]
        tm.assert_frame_equal(result, expected)

        # Slicing date on first level should break (of course)
        with assertRaises(KeyError):
            df_swap.loc['2016-01-01']

        # GH12685 (partial string with daily resolution or below)
        dr = date_range('2013-01-01', periods=100, freq='D')
        ix = MultiIndex.from_product([dr, ['a', 'b']])
        df = DataFrame(np.random.randn(200, 1), columns=['A'], index=ix)

        result = df.loc[idx['2013-03':'2013-03', :], :]
        expected = df.iloc[118:180]
        tm.assert_frame_equal(result, expected)

    def test_rangeindex_fallback_coercion_bug(self):
        # GH 12893
        foo = pd.DataFrame(np.arange(100).reshape((10, 10)))
        bar = pd.DataFrame(np.arange(100).reshape((10, 10)))
        df = pd.concat({'foo': foo.stack(), 'bar': bar.stack()}, axis=1)
        df.index.names = ['fizz', 'buzz']

        str(df)
        expected = pd.DataFrame({'bar': np.arange(100),
                                 'foo': np.arange(100)},
                                index=pd.MultiIndex.from_product(
                                    [range(10), range(10)],
                                    names=['fizz', 'buzz']))
        tm.assert_frame_equal(df, expected, check_like=True)

        result = df.index.get_level_values('fizz')
        expected = pd.Int64Index(np.arange(10), name='fizz').repeat(10)
        tm.assert_index_equal(result, expected)

        result = df.index.get_level_values('buzz')
        expected = pd.Int64Index(np.tile(np.arange(10), 10), name='buzz')
        tm.assert_index_equal(result, expected)

    def test_dropna(self):
        # GH 6194
        idx = pd.MultiIndex.from_arrays([[1, np.nan, 3, np.nan, 5],
                                         [1, 2, np.nan, np.nan, 5],
                                         ['a', 'b', 'c', np.nan, 'e']])

        exp = pd.MultiIndex.from_arrays([[1, 5],
                                         [1, 5],
                                         ['a', 'e']])
        tm.assert_index_equal(idx.dropna(), exp)
        tm.assert_index_equal(idx.dropna(how='any'), exp)

        exp = pd.MultiIndex.from_arrays([[1, np.nan, 3, 5],
                                         [1, 2, np.nan, 5],
                                         ['a', 'b', 'c', 'e']])
        tm.assert_index_equal(idx.dropna(how='all'), exp)

        msg = "invalid how option: xxx"
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.dropna(how='xxx')
