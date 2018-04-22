# -*- coding: utf-8 -*-

from itertools import product

import pytest

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)
from pandas.compat import range, u

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestDuplicates(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setup_method(self, method):
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

    @pytest.mark.parametrize('names', [['a', 'b', 'a'], [1, 1, 2],
                                       [1, 'a', 1]])
    def test_duplicate_level_names(self, names):
        # GH18872
        pytest.raises(ValueError, pd.MultiIndex.from_product,
                      [[0, 1]] * 3, names=names)

        # With .rename()
        mi = pd.MultiIndex.from_product([[0, 1]] * 3)
        tm.assert_raises_regex(ValueError, "Duplicated level name:",
                               mi.rename, names)

        # With .rename(., level=)
        mi.rename(names[0], level=1, inplace=True)
        tm.assert_raises_regex(ValueError, "Duplicated level name:",
                               mi.rename, names[:2], level=[0, 2])

    def test_duplicates(self):
        assert not self.index.has_duplicates
        assert self.index.append(self.index).has_duplicates

        index = MultiIndex(levels=[[0, 1], [0, 1, 2]], labels=[
            [0, 0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 0, 1, 2]])
        assert index.has_duplicates

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
        assert not index.has_duplicates

        # handle int64 overflow if possible
        def check(nlevels, with_nulls):
            labels = np.tile(np.arange(500), 2)
            level = np.arange(500)

            if with_nulls:  # inject some null values
                labels[500] = -1  # common nan value
                labels = [labels.copy() for i in range(nlevels)]
                for i in range(nlevels):
                    labels[i][500 + i - nlevels // 2] = -1

                labels += [np.array([-1, 1]).repeat(500)]
            else:
                labels = [labels] * nlevels + [np.arange(2).repeat(500)]

            levels = [level] * nlevels + [[0, 1]]

            # no dups
            index = MultiIndex(levels=levels, labels=labels)
            assert not index.has_duplicates

            # with a dup
            if with_nulls:
                f = lambda a: np.insert(a, 1000, a[0])
                labels = list(map(f, labels))
                index = MultiIndex(levels=levels, labels=labels)
            else:
                values = index.values.tolist()
                index = MultiIndex.from_tuples(values + [values[0]])

            assert index.has_duplicates

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
            right = pd._libs.hashtable.duplicated_object(mi.values, keep=keep)
            tm.assert_numpy_array_equal(left, right)

        # GH5873
        for a in [101, 102]:
            mi = MultiIndex.from_arrays([[101, a], [3.5, np.nan]])
            assert not mi.has_duplicates
            assert mi.get_duplicates() == []
            tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
                2, dtype='bool'))

        for n in range(1, 6):  # 1st level shape
            for m in range(1, 5):  # 2nd level shape
                # all possible unique combinations, including nan
                lab = product(range(-1, n), range(-1, m))
                mi = MultiIndex(levels=[list('abcde')[:n], list('WXYZ')[:m]],
                                labels=np.random.permutation(list(lab)).T)
                assert len(mi) == (n + 1) * (m + 1)
                assert not mi.has_duplicates
                assert mi.get_duplicates() == []
                tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
                    len(mi), dtype='bool'))

    def test_duplicate_meta_data(self):
        # GH 10115
        index = MultiIndex(
            levels=[[0, 1], [0, 1, 2]],
            labels=[[0, 0, 0, 0, 1, 1, 1],
                    [0, 1, 2, 0, 0, 1, 2]])

        for idx in [index,
                    index.set_names([None, None]),
                    index.set_names([None, 'Num']),
                    index.set_names(['Upper', 'Num']), ]:
            assert idx.has_duplicates
            assert idx.drop_duplicates().names == idx.names

    def test_duplicate_multiindex_labels(self):
        # GH 17464
        # Make sure that a MultiIndex with duplicate levels throws a ValueError
        with pytest.raises(ValueError):
            ind = pd.MultiIndex([['A'] * 10, range(10)], [[0] * 10, range(10)])

        # And that using set_levels with duplicate levels fails
        ind = MultiIndex.from_arrays([['A', 'A', 'B', 'B', 'B'],
                                      [1, 2, 1, 2, 3]])
        with pytest.raises(ValueError):
            ind.set_levels([['A', 'B', 'A', 'A', 'B'], [2, 1, 3, -2, 5]],
                           inplace=True)
