import itertools
from warnings import catch_warnings

import numpy as np
import pytest

from pandas import DataFrame, Index, MultiIndex, Series
from pandas.util import testing as tm


@pytest.fixture
def single_level_multiindex():
    """single level MultiIndex"""
    return MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                      codes=[[0, 1, 2, 3]], names=['first'])


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
class TestMultiIndexLoc(object):

    def test_loc_getitem_series(self):
        # GH14730
        # passing a series as a key with a MultiIndex
        index = MultiIndex.from_product([[1, 2, 3], ['A', 'B', 'C']])
        x = Series(index=index, data=range(9), dtype=np.float64)
        y = Series([1, 3])
        expected = Series(
            data=[0, 1, 2, 6, 7, 8],
            index=MultiIndex.from_product([[1, 3], ['A', 'B', 'C']]),
            dtype=np.float64)
        result = x.loc[y]
        tm.assert_series_equal(result, expected)

        result = x.loc[[1, 3]]
        tm.assert_series_equal(result, expected)

        # GH15424
        y1 = Series([1, 3], index=[1, 2])
        result = x.loc[y1]
        tm.assert_series_equal(result, expected)

        empty = Series(data=[], dtype=np.float64)
        expected = Series([], index=MultiIndex(
            levels=index.levels, codes=[[], []], dtype=np.float64))
        result = x.loc[empty]
        tm.assert_series_equal(result, expected)

    def test_loc_getitem_array(self):
        # GH15434
        # passing an array as a key with a MultiIndex
        index = MultiIndex.from_product([[1, 2, 3], ['A', 'B', 'C']])
        x = Series(index=index, data=range(9), dtype=np.float64)
        y = np.array([1, 3])
        expected = Series(
            data=[0, 1, 2, 6, 7, 8],
            index=MultiIndex.from_product([[1, 3], ['A', 'B', 'C']]),
            dtype=np.float64)
        result = x.loc[y]
        tm.assert_series_equal(result, expected)

        # empty array:
        empty = np.array([])
        expected = Series([], index=MultiIndex(
            levels=index.levels, codes=[[], []], dtype=np.float64))
        result = x.loc[empty]
        tm.assert_series_equal(result, expected)

        # 0-dim array (scalar):
        scalar = np.int64(1)
        expected = Series(
            data=[0, 1, 2],
            index=['A', 'B', 'C'],
            dtype=np.float64)
        result = x.loc[scalar]
        tm.assert_series_equal(result, expected)

    def test_loc_multiindex(self):

        mi_labels = DataFrame(np.random.randn(3, 3),
                              columns=[['i', 'i', 'j'], ['A', 'A', 'B']],
                              index=[['i', 'i', 'j'], ['X', 'X', 'Y']])

        mi_int = DataFrame(np.random.randn(3, 3),
                           columns=[[2, 2, 4], [6, 8, 10]],
                           index=[[4, 4, 8], [8, 10, 12]])

        # the first row
        rs = mi_labels.loc['i']
        with catch_warnings(record=True):
            xp = mi_labels.ix['i']
        tm.assert_frame_equal(rs, xp)

        # 2nd (last) columns
        rs = mi_labels.loc[:, 'j']
        with catch_warnings(record=True):
            xp = mi_labels.ix[:, 'j']
        tm.assert_frame_equal(rs, xp)

        # corner column
        rs = mi_labels.loc['j'].loc[:, 'j']
        with catch_warnings(record=True):
            xp = mi_labels.ix['j'].ix[:, 'j']
        tm.assert_frame_equal(rs, xp)

        # with a tuple
        rs = mi_labels.loc[('i', 'X')]
        with catch_warnings(record=True):
            xp = mi_labels.ix[('i', 'X')]
        tm.assert_frame_equal(rs, xp)

        rs = mi_int.loc[4]
        with catch_warnings(record=True):
            xp = mi_int.ix[4]
        tm.assert_frame_equal(rs, xp)

        # missing label
        pytest.raises(KeyError, lambda: mi_int.loc[2])
        with catch_warnings(record=True):
            # GH 21593
            pytest.raises(KeyError, lambda: mi_int.ix[2])

    def test_loc_multiindex_indexer_none(self):

        # GH6788
        # multi-index indexer is None (meaning take all)
        attributes = ['Attribute' + str(i) for i in range(1)]
        attribute_values = ['Value' + str(i) for i in range(5)]

        index = MultiIndex.from_product([attributes, attribute_values])
        df = 0.1 * np.random.randn(10, 1 * 5) + 0.5
        df = DataFrame(df, columns=index)
        result = df[attributes]
        tm.assert_frame_equal(result, df)

        # GH 7349
        # loc with a multi-index seems to be doing fallback
        df = DataFrame(np.arange(12).reshape(-1, 1),
                       index=MultiIndex.from_product([[1, 2, 3, 4],
                                                      [1, 2, 3]]))

        expected = df.loc[([1, 2], ), :]
        result = df.loc[[1, 2]]
        tm.assert_frame_equal(result, expected)

    def test_loc_multiindex_incomplete(self):

        # GH 7399
        # incomplete indexers
        s = Series(np.arange(15, dtype='int64'),
                   MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.loc[:, 'a':'c']

        result = s.loc[0:4, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        result = s.loc[:4, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        result = s.loc[0:, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        # GH 7400
        # multiindexer gettitem with list of indexers skips wrong element
        s = Series(np.arange(15, dtype='int64'),
                   MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.iloc[[6, 7, 8, 12, 13, 14]]
        result = s.loc[2:4:2, 'a':'c']
        tm.assert_series_equal(result, expected)

    def test_get_loc_single_level(self, single_level_multiindex):
        single_level = single_level_multiindex
        s = Series(np.random.randn(len(single_level)),
                   index=single_level)
        for k in single_level.values:
            s[k]

    def test_loc_getitem_int_slice(self):
        # GH 3053
        # loc should treat integer slices like label slices

        index = MultiIndex.from_tuples([t for t in itertools.product(
            [6, 7, 8], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[6:8, :]
        expected = df
        tm.assert_frame_equal(result, expected)

        index = MultiIndex.from_tuples([t
                                        for t in itertools.product(
                                            [10, 20, 30], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[20:30, :]
        expected = df.iloc[2:]
        tm.assert_frame_equal(result, expected)

        # doc examples
        result = df.loc[10, :]
        expected = df.iloc[0:2]
        expected.index = ['a', 'b']
        tm.assert_frame_equal(result, expected)

        result = df.loc[:, 10]
        # expected = df.ix[:,10] (this fails)
        expected = df[10]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        'indexer_type_1',
        (list, tuple, set, slice, np.ndarray, Series, Index))
    @pytest.mark.parametrize(
        'indexer_type_2',
        (list, tuple, set, slice, np.ndarray, Series, Index))
    def test_loc_getitem_nested_indexer(self, indexer_type_1, indexer_type_2):
        # GH #19686
        # .loc should work with nested indexers which can be
        # any list-like objects (see `pandas.api.types.is_list_like`) or slices

        def convert_nested_indexer(indexer_type, keys):
            if indexer_type == np.ndarray:
                return np.array(keys)
            if indexer_type == slice:
                return slice(*keys)
            return indexer_type(keys)

        a = [10, 20, 30]
        b = [1, 2, 3]
        index = MultiIndex.from_product([a, b])
        df = DataFrame(
            np.arange(len(index), dtype='int64'),
            index=index, columns=['Data'])

        keys = ([10, 20], [2, 3])
        types = (indexer_type_1, indexer_type_2)

        # check indexers with all the combinations of nested objects
        # of all the valid types
        indexer = tuple(
            convert_nested_indexer(indexer_type, k)
            for indexer_type, k in zip(types, keys))

        result = df.loc[indexer, 'Data']
        expected = Series(
            [1, 2, 4, 5], name='Data',
            index=MultiIndex.from_product(keys))

        tm.assert_series_equal(result, expected)
