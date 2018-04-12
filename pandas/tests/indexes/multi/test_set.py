# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd

from pandas import (CategoricalIndex, MultiIndex)
from pandas.compat import range

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestSet(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_set_name_methods(self):
        # so long as these are synonyms, we don't need to test set_names
        assert self.index.rename == self.index.set_names
        new_names = [name + "SUFFIX" for name in self.index_names]
        ind = self.index.set_names(new_names)
        assert self.index.names == self.index_names
        assert ind.names == new_names
        with tm.assert_raises_regex(ValueError, "^Length"):
            ind.set_names(new_names + new_names)
        new_names2 = [name + "SUFFIX2" for name in new_names]
        res = ind.set_names(new_names2, inplace=True)
        assert res is None
        assert ind.names == new_names2

        # set names for specific level (# GH7792)
        ind = self.index.set_names(new_names[0], level=0)
        assert self.index.names == self.index_names
        assert ind.names == [new_names[0], self.index_names[1]]

        res = ind.set_names(new_names2[0], level=0, inplace=True)
        assert res is None
        assert ind.names == [new_names2[0], self.index_names[1]]

        # set names for multiple levels
        ind = self.index.set_names(new_names, level=[0, 1])
        assert self.index.names == self.index_names
        assert ind.names == new_names

        res = ind.set_names(new_names2, level=[0, 1], inplace=True)
        assert res is None
        assert ind.names == new_names2

    def test_set_levels_labels_directly(self):
        # setting levels/labels directly raises AttributeError

        levels = self.index.levels
        new_levels = [[lev + 'a' for lev in level] for level in levels]

        labels = self.index.labels
        major_labels, minor_labels = labels
        major_labels = [(x + 1) % 3 for x in major_labels]
        minor_labels = [(x + 1) % 1 for x in minor_labels]
        new_labels = [major_labels, minor_labels]

        with pytest.raises(AttributeError):
            self.index.levels = new_levels

        with pytest.raises(AttributeError):
            self.index.labels = new_labels

    def test_set_levels(self):
        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        levels = self.index.levels
        new_levels = [[lev + 'a' for lev in level] for level in levels]

        def assert_matching(actual, expected, check_dtype=False):
            # avoid specifying internal representation
            # as much as possible
            assert len(actual) == len(expected)
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
        assert inplace_return is None
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
        assert inplace_return is None
        assert_matching(ind2.levels, [new_levels[0], levels[1]])
        assert_matching(self.index.levels, levels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels[1], level=1, inplace=True)
        assert inplace_return is None
        assert_matching(ind2.levels, [levels[0], new_levels[1]])
        assert_matching(self.index.levels, levels)

        # level changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, level=[0, 1],
                                         inplace=True)
        assert inplace_return is None
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # illegal level changing should not change levels
        # GH 13754
        original_index = self.index.copy()
        for inplace in [True, False]:
            with tm.assert_raises_regex(ValueError, "^On"):
                self.index.set_levels(['c'], level=0, inplace=inplace)
            assert_matching(self.index.levels, original_index.levels,
                            check_dtype=True)

            with tm.assert_raises_regex(ValueError, "^On"):
                self.index.set_labels([0, 1, 2, 3, 4, 5], level=0,
                                      inplace=inplace)
            assert_matching(self.index.labels, original_index.labels,
                            check_dtype=True)

            with tm.assert_raises_regex(TypeError, "^Levels"):
                self.index.set_levels('c', level=0, inplace=inplace)
            assert_matching(self.index.levels, original_index.levels,
                            check_dtype=True)

            with tm.assert_raises_regex(TypeError, "^Labels"):
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
            assert len(actual) == len(expected)
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
        assert inplace_return is None
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
        assert inplace_return is None
        assert_matching(ind2.labels, [new_labels[0], labels[1]])
        assert_matching(self.index.labels, labels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels[1], level=1, inplace=True)
        assert inplace_return is None
        assert_matching(ind2.labels, [labels[0], new_labels[1]])
        assert_matching(self.index.labels, labels)

        # label changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, level=[0, 1],
                                         inplace=True)
        assert inplace_return is None
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing for levels of different magnitude of categories
        ind = pd.MultiIndex.from_tuples([(0, i) for i in range(130)])
        new_labels = range(129, -1, -1)
        expected = pd.MultiIndex.from_tuples(
            [(0, i) for i in new_labels])

        # [w/o mutation]
        result = ind.set_labels(labels=new_labels, level=1)
        assert result.equals(expected)

        # [w/ mutation]
        result = ind.copy()
        result.set_labels(labels=new_labels, level=1, inplace=True)
        assert result.equals(expected)

    def test_set_levels_labels_names_bad_input(self):
        levels, labels = self.index.levels, self.index.labels
        names = self.index.names

        with tm.assert_raises_regex(ValueError, 'Length of levels'):
            self.index.set_levels([levels[0]])

        with tm.assert_raises_regex(ValueError, 'Length of labels'):
            self.index.set_labels([labels[0]])

        with tm.assert_raises_regex(ValueError, 'Length of names'):
            self.index.set_names([names[0]])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assert_raises_regex(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assert_raises_regex(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assert_raises_regex(TypeError, 'list-like'):
            self.index.set_names(names[0])

        # should have equal lengths
        with tm.assert_raises_regex(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0], level=[0, 1])

        with tm.assert_raises_regex(TypeError, 'list-like'):
            self.index.set_levels(levels, level=0)

        # should have equal lengths
        with tm.assert_raises_regex(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0], level=[0, 1])

        with tm.assert_raises_regex(TypeError, 'list-like'):
            self.index.set_labels(labels, level=0)

        # should have equal lengths
        with tm.assert_raises_regex(ValueError, 'Length of names'):
            self.index.set_names(names[0], level=[0, 1])

        with tm.assert_raises_regex(TypeError, 'string'):
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

    def test_set_value_keeps_names(self):
        # motivating example from #3742
        lev1 = ['hans', 'hans', 'hans', 'grethe', 'grethe', 'grethe']
        lev2 = ['1', '2', '3'] * 2
        idx = pd.MultiIndex.from_arrays([lev1, lev2], names=['Name', 'Number'])
        df = pd.DataFrame(
            np.random.randn(6, 4),
            columns=['one', 'two', 'three', 'four'],
            index=idx)
        df = df.sort_index()
        assert df._is_copy is None
        assert df.index.names == ('Name', 'Number')
        df.at[('grethe', '4'), 'one'] = 99.34
        assert df._is_copy is None
        assert df.index.names == ('Name', 'Number')
