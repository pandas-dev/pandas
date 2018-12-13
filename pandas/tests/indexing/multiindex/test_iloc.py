import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Series
from pandas.util import testing as tm


@pytest.mark.parametrize('indexer, expected', [
    (lambda df: df.iloc[0],
     lambda arr, df: Series(arr[0], index=df.columns, name=(4, 8))),
    (lambda df: df.iloc[2],
     lambda arr, df: Series(arr[2], index=df.columns, name=(8, 12))),
    (lambda df: df.iloc[:, 2],
     lambda arr, df: Series(arr[:, 2], index=df.index, name=(4, 10))),
    (lambda df: df.iloc[2, 2],
     lambda arr, df: arr[2, 2]),
    (lambda df: df.iloc[[0, 1]],
     lambda arr, df: df.xs(4, drop_level=False))
])
def test_iloc_getitem(indexer, expected):
    arr = np.random.randn(3, 3)
    df = DataFrame(arr, columns=[[2, 2, 4], [6, 8, 10]],
                   index=[[4, 4, 8], [8, 10, 12]])

    result = indexer(df)
    expected = expected(arr, df)

    try:
        tm.assert_equal(result, expected)
    except NotImplementedError:
        assert result == expected


def test_iloc_getitem_multiple_items():
    # GH 5528
    tup = zip(*[['a', 'a', 'b', 'b'], ['x', 'y', 'x', 'y']])
    index = MultiIndex.from_tuples(tup)
    df = DataFrame(np.random.randn(4, 4), index=index)
    result = df.iloc[[2, 3]]
    expected = df.xs('b', drop_level=False)
    tm.assert_frame_equal(result, expected)


def test_iloc_getitem_labels():
    # this is basically regular indexing
    arr = np.random.randn(4, 3)
    df = DataFrame(arr,
                   columns=[['i', 'i', 'j'], ['A', 'A', 'B']],
                   index=[['i', 'i', 'j', 'k'], ['X', 'X', 'Y', 'Y']])
    result = df.iloc[2, 2]
    expected = arr[2, 2]
    assert result == expected


def test_frame_getitem_slice(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    result = frame.iloc[:4]
    expected = frame[:4]
    tm.assert_frame_equal(result, expected)


def test_frame_setitem_slice(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data
    df.iloc[:4] = 0

    assert (df.values[:4] == 0).all()
    assert (df.values[4:] != 0).all()


def test_indexing_ambiguity_bug_1678():
    columns = MultiIndex.from_tuples(
        [('Ohio', 'Green'), ('Ohio', 'Red'), ('Colorado', 'Green')])
    index = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1), ('b', 2)])

    df = DataFrame(np.arange(12).reshape((4, 3)), index=index, columns=columns)

    result = df.iloc[:, 1]
    expected = df.loc[:, ('Ohio', 'Red')]
    tm.assert_series_equal(result, expected)


def test_iloc_mi():
    # GH 13797
    # Test if iloc can handle integer locations in MultiIndexed DataFrame

    data = [['str00', 'str01'], ['str10', 'str11'], ['str20', 'srt21'],
            ['str30', 'str31'], ['str40', 'str41']]

    index = MultiIndex.from_tuples(
        [('CC', 'A'), ('CC', 'B'), ('CC', 'B'), ('BB', 'a'), ('BB', 'b')])

    expected = DataFrame(data)
    df = DataFrame(data, index=index)

    result = DataFrame([[df.iloc[r, c] for c in range(2)] for r in range(5)])

    tm.assert_frame_equal(result, expected)
