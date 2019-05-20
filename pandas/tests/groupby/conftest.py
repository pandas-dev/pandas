import numpy as np
import pytest

from pandas import DataFrame, CategoricalIndex, Index, MultiIndex
from pandas.util import testing as tm


@pytest.fixture
def mframe():
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                              'three']],
                       codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                              [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
    return DataFrame(np.random.randn(10, 3), index=index,
                     columns=['A', 'B', 'C'])


@pytest.fixture
def df():
    return DataFrame(
        {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
         'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
         'C': np.random.randn(8),
         'D': np.random.randn(8)})


@pytest.fixture
def ts():
    return tm.makeTimeSeries()


@pytest.fixture
def seriesd():
    return tm.getSeriesData()


@pytest.fixture
def tsd():
    return tm.getTimeSeriesData()


@pytest.fixture
def frame(seriesd):
    return DataFrame(seriesd)


@pytest.fixture
def tsframe(tsd):
    return DataFrame(tsd)


@pytest.fixture
def df_mixed_floats():
    return DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'foo', 'foo'],
                      'B': ['one', 'one', 'two', 'three',
                            'two', 'two', 'one', 'three'],
                      'C': np.random.randn(8),
                      'D': np.array(
                          np.random.randn(8), dtype='float32')})


@pytest.fixture
def three_group():
    return DataFrame({'A': ['foo', 'foo', 'foo',
                            'foo', 'bar', 'bar',
                            'bar', 'bar',
                            'foo', 'foo', 'foo'],
                      'B': ['one', 'one', 'one',
                            'two', 'one', 'one', 'one', 'two',
                            'two', 'two', 'one'],
                      'C': ['dull', 'dull', 'shiny',
                            'dull', 'dull', 'shiny', 'shiny',
                            'dull', 'shiny', 'shiny', 'shiny'],
                      'D': np.random.randn(11),
                      'E': np.random.randn(11),
                      'F': np.random.randn(11)})


@pytest.fixture
def df_cat():
    df = DataFrame({'a': ['x', 'x', 'x', 'y'],
                    'b': ['a', 'a', 'b', 'a'],
                    'c': [1, 2, 3, 4]})
    df['a'] = df['a'].astype('category')
    df['b'] = df['b'].astype('category')
    return df


@pytest.fixture
def multi_index_cat_complete():
    lvls = [CategoricalIndex(['x', 'y'], categories=['x', 'y'], ordered=False),
            CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=False)]
    index = MultiIndex.from_product(lvls, names=['a', 'b'])
    return index


@pytest.fixture
def multi_index_cat_partial(df_cat):
    return MultiIndex.from_frame(df_cat[['a', 'b']].drop_duplicates())


@pytest.fixture
def multi_index_non_cat_partial():
    return MultiIndex.from_tuples([('x', 'a'), ('x', 'b'), ('y', 'a')],
                                  names=('a', 'b'))


@pytest.fixture
def multi_index_cat_compl_dict():
    lvls = [CategoricalIndex(['x', 'y'], categories=['x', 'y'], ordered=False),
            CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=False),
            Index(['min', 'max'])]
    index = MultiIndex.from_product(lvls, names=['a', 'b', None])
    return index


@pytest.fixture
def multi_index_non_cat_partial_dict():
    return MultiIndex.from_tuples([('x', 'a', 'min'), ('x', 'a', 'max'),
                                   ('x', 'b', 'min'), ('x', 'b', 'max'),
                                   ('y', 'a', 'min'), ('y', 'a', 'max')],
                                  names=('a', 'b', None))
