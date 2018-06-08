# -*- coding: utf-8 -*-


import pandas.util.testing as tm


def assert_multiindex_copied(copy, original):
    # Levels should be (at least, shallow copied)
    tm.assert_copy(copy.levels, original.levels)
    tm.assert_almost_equal(copy.labels, original.labels)

    # Labels doesn't matter which way copied
    tm.assert_almost_equal(copy.labels, original.labels)
    assert copy.labels is not original.labels

    # Names doesn't matter which way copied
    assert copy.names == original.names
    assert copy.names is not original.names

    # Sort order should be copied
    assert copy.sortorder == original.sortorder


def test_copy(_index):
    i_copy = _index.copy()

    assert_multiindex_copied(i_copy, _index)


def test_shallow_copy(_index):
    i_copy = _index._shallow_copy()

    assert_multiindex_copied(i_copy, _index)


def test_view(_index):
    i_view = _index.view()
    assert_multiindex_copied(i_view, _index)
