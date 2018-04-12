# -*- coding: utf-8 -*-

import re

import numpy as np

from pandas import (Index, MultiIndex)

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestConstructor(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_copy_in_constructor(self):
        levels = np.array(["a", "b", "c"])
        labels = np.array([1, 1, 2, 0, 0, 1, 1])
        val = labels[0]
        mi = MultiIndex(levels=[levels, levels], labels=[labels, labels],
                        copy=True)
        assert mi.labels[0][0] == val
        labels[0] = 15
        assert mi.labels[0][0] == val
        val = levels[0]
        levels[0] = "PANDA"
        assert mi.levels[0][0] == val

    def test_constructor_single_level(self):
        result = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                            labels=[[0, 1, 2, 3]], names=['first'])
        assert isinstance(result, MultiIndex)
        expected = Index(['foo', 'bar', 'baz', 'qux'], name='first')
        tm.assert_index_equal(result.levels[0], expected)
        assert result.names == ['first']

    def test_constructor_no_levels(self):
        tm.assert_raises_regex(ValueError, "non-zero number "
                                           "of levels/labels",
                               MultiIndex, levels=[], labels=[])
        both_re = re.compile('Must pass both levels and labels')
        with tm.assert_raises_regex(TypeError, both_re):
            MultiIndex(levels=[])
        with tm.assert_raises_regex(TypeError, both_re):
            MultiIndex(labels=[])

    def test_constructor_mismatched_label_levels(self):
        labels = [np.array([1]), np.array([2]), np.array([3])]
        levels = ["a"]
        tm.assert_raises_regex(ValueError, "Length of levels and labels "
                                           "must be the same", MultiIndex,
                               levels=levels, labels=labels)
        length_error = re.compile('>= length of level')
        label_error = re.compile(r'Unequal label lengths: \[4, 2\]')

        # important to check that it's looking at the right thing.
        with tm.assert_raises_regex(ValueError, length_error):
            MultiIndex(levels=[['a'], ['b']],
                       labels=[[0, 1, 2, 3], [0, 3, 4, 1]])

        with tm.assert_raises_regex(ValueError, label_error):
            MultiIndex(levels=[['a'], ['b']], labels=[[0, 0, 0, 0], [0, 0]])

        # external API
        with tm.assert_raises_regex(ValueError, length_error):
            self.index.copy().set_levels([['a'], ['b']])

        with tm.assert_raises_regex(ValueError, label_error):
            self.index.copy().set_labels([[0, 0, 0, 0], [0, 0]])
