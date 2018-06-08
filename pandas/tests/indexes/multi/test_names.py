# -*- coding: utf-8 -*-


import pandas as pd
import pandas.util.testing as tm
from pandas import MultiIndex


def check_level_names(index, names):
    assert [level.name for level in index.levels] == list(names)


def test_slice_keep_name():
    x = MultiIndex.from_tuples([('a', 'b'), (1, 2), ('c', 'd')],
                               names=['x', 'y'])
    assert x[1:].names == x.names


def test_index_name_retained():
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


def test_changing_names(_index):

    # names should be applied to levels
    level_names = [level.name for level in _index.levels]
    check_level_names(_index, _index.names)

    view = _index.view()
    copy = _index.copy()
    shallow_copy = _index._shallow_copy()

    # changing names should change level names on object
    new_names = [name + "a" for name in _index.names]
    _index.names = new_names
    check_level_names(_index, new_names)

    # but not on copies
    check_level_names(view, level_names)
    check_level_names(copy, level_names)
    check_level_names(shallow_copy, level_names)

    # and copies shouldn't change original
    shallow_copy.names = [name + "c" for name in shallow_copy.names]
    check_level_names(_index, new_names)


def test_take_preserve_name(_index):
    taken = _index.take([3, 0, 1])
    assert taken.names == _index.names


def test_copy_names():
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


def test_names(_index, index_names):

    # names are assigned in setup
    names = index_names
    level_names = [level.name for level in _index.levels]
    assert names == level_names

    # setting bad names on existing
    index = _index
    tm.assert_raises_regex(ValueError, "^Length of names",
                           setattr, index, "names",
                           list(index.names) + ["third"])
    tm.assert_raises_regex(ValueError, "^Length of names",
                           setattr, index, "names", [])

    # initializing with bad names (should always be equivalent)
    major_axis, minor_axis = _index.levels
    major_labels, minor_labels = _index.labels
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
