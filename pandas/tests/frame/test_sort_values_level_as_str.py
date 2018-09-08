import numpy as np
import pytest

from pandas import DataFrame, Categorical
from pandas.errors import PerformanceWarning
from pandas.util import testing as tm
from pandas.util.testing import assert_frame_equal


@pytest.fixture
def df_none():
    return DataFrame({
        'outer': ['a', 'a', 'a', 'b', 'b', 'b'],
        'inner': [1, 2, 2, 2, 1, 1],
        'A': np.arange(6, 0, -1),
        ('B', 5): ['one', 'one', 'two', 'two', 'one', 'one']})


@pytest.fixture
def df_category_with_nan():
    return DataFrame({
        'c': Categorical(['A', np.nan, 'B'],
                         categories=['A', 'B'],
                         ordered=True)})


@pytest.fixture(params=[
    ['outer'],
    ['outer', 'inner']
])
def df_idx(request, df_none):
    levels = request.param
    return df_none.set_index(levels)


@pytest.fixture(params=[
    'inner',     # index level
    ['outer'],   # list of index level
    'A',         # column
    [('B', 5)],  # list of column
    ['inner', 'outer'],   # two index levels
    [('B', 5), 'outer'],  # index level and column
    ['A', ('B', 5)],      # Two columns
    ['inner', 'outer']    # two index levels and column
])
def sort_names(request):
    return request.param


@pytest.fixture(params=[True, False])
def ascending(request):
    return request.param


def test_sort_index_level_and_column_label(
        df_none, df_idx, sort_names, ascending):

    # GH 14353

    # Get index levels from df_idx
    levels = df_idx.index.names

    # Compute expected by sorting on columns and the setting index
    expected = df_none.sort_values(by=sort_names,
                                   ascending=ascending,
                                   axis=0).set_index(levels)

    # Compute result sorting on mix on columns and index levels
    result = df_idx.sort_values(by=sort_names,
                                ascending=ascending,
                                axis=0)

    assert_frame_equal(result, expected)


def test_sort_column_level_and_index_label(
        df_none, df_idx, sort_names, ascending):

    # GH 14353

    # Get levels from df_idx
    levels = df_idx.index.names

    # Compute expected by sorting on axis=0, setting index levels, and then
    # transposing. For some cases this will result in a frame with
    # multiple column levels
    expected = df_none.sort_values(by=sort_names,
                                   ascending=ascending,
                                   axis=0).set_index(levels).T

    # Compute result by transposing and sorting on axis=1.
    result = df_idx.T.sort_values(by=sort_names,
                                  ascending=ascending,
                                  axis=1)

    if len(levels) > 1:
        # Accessing multi-level columns that are not lexsorted raises a
        # performance warning
        with tm.assert_produces_warning(PerformanceWarning,
                                        check_stacklevel=False):
            assert_frame_equal(result, expected)
    else:
        assert_frame_equal(result, expected)


def test_sort_index_na_position_with_categories(
        df_category_with_nan, sort_names):
    result = df_category_with_nan.sort_values(by='c', na_position='first')
    expected = DataFrame({
        'c': Categorical([np.nan, 'A', 'B'],
                         categories=['A', 'B'],
                         ordered=True)})
    
    assert_frame_equal(result, expected)

    result = df_category_with_nan.sort_values(by='c', na_position='last')
    expected = DataFrame({
        'c': Categorical(['A', 'B', np.nan],
                         categories=['A', 'B'],
                         ordered=True)})
    
    assert_frame_equal(result, expected)