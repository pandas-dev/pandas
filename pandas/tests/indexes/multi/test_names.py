# -*- coding: utf-8 -*-

import pytest

import pandas as pd

from pandas import MultiIndex

import pandas.util.testing as tm

from .common import Base


class TestNames(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_copy_names(self):
        # Check that adding a "names" parameter to the copy is honored
        # GH14302
        multi_idx = pd.Index([(1, 2), (3, 4)], names=['MyName1', 'MyName2'])
        multi_idx1 = multi_idx.copy()

        assert multi_idx.equals(multi_idx1)
        assert multi_idx.names == ['MyName1', 'MyName2']
        assert multi_idx1.names == ['MyName1', 'MyName2']

        multi_idx2 = multi_idx.copy(names=['NewName1', 'NewName2'])

        assert multi_idx.equals(multi_idx2)
        assert multi_idx.names == ['MyName1', 'MyName2']
        assert multi_idx2.names == ['NewName1', 'NewName2']

        multi_idx3 = multi_idx.copy(name=['NewName1', 'NewName2'])

        assert multi_idx.equals(multi_idx3)
        assert multi_idx.names == ['MyName1', 'MyName2']
        assert multi_idx3.names == ['NewName1', 'NewName2']

    def test_names(self):

        # names are assigned in setup
        names = self.index_names
        level_names = [level.name for level in self.index.levels]
        assert names == level_names

        # setting bad names on existing
        index = self.index
        tm.assert_raises_regex(ValueError, "^Length of names",
                               setattr, index, "names",
                               list(index.names) + ["third"])
        tm.assert_raises_regex(ValueError, "^Length of names",
                               setattr, index, "names", [])

        # initializing with bad names (should always be equivalent)
        major_axis, minor_axis = self.index.levels
        major_labels, minor_labels = self.index.labels
        tm.assert_raises_regex(ValueError, "^Length of names", MultiIndex,
                               levels=[major_axis, minor_axis],
                               labels=[major_labels, minor_labels],
                               names=['first'])
        tm.assert_raises_regex(ValueError, "^Length of names", MultiIndex,
                               levels=[major_axis, minor_axis],
                               labels=[major_labels, minor_labels],
                               names=['first', 'second', 'third'])

        # names are assigned
        index.names = ["a", "b"]
        ind_names = list(index.names)
        level_names = [level.name for level in index.levels]
        assert ind_names == level_names

    def check_level_names(self, index, names):
        assert [level.name for level in index.levels] == list(names)

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

    def test_take_preserve_name(self):
        taken = self.index.take([3, 0, 1])
        assert taken.names == self.index.names

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

    def test_tuples_with_name_string(self):
        # GH 15110 and GH 14848

        li = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
        with pytest.raises(ValueError):
            pd.Index(li, name='abc')
        with pytest.raises(ValueError):
            pd.Index(li, name='a')
