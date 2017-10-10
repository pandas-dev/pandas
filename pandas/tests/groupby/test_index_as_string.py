import pytest
import pandas as pd
import numpy as np

from pandas.util.testing import assert_frame_equal, assert_series_equal
import pandas.util.testing as tm


def build_df_multi():
    idx = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('a', 3),
                                     ('b', 1), ('b', 2), ('b', 3)])
    idx.names = ['outer', 'inner']
    df_multi = pd.DataFrame({"A": np.arange(6),
                             'B': ['one', 'one', 'two',
                                   'two', 'one', 'one']},
                            index=idx)
    return df_multi


def build_df_single():
    df_single = build_df_multi().reset_index('outer')
    return df_single


def build_test_series():
    series_multi = build_df_multi().set_index('B', append=True)['A']
    return series_multi


class TestGroupByIndexAsString(object):

    @pytest.mark.parametrize('frame', [build_df_multi(), build_df_single()])
    def test_grouper_index_level_as_string(self, frame):
        # Column and Index
        result = frame.groupby(['B', 'inner']).mean()
        expected = frame.groupby(['B', pd.Grouper(level='inner')]).mean()
        assert_frame_equal(result, expected)

        # Index and Column
        result = frame.groupby(['inner', 'B']).mean()
        expected = frame.groupby([pd.Grouper(level='inner'), 'B']).mean()
        assert_frame_equal(result, expected)

        # Single element list of Index
        result = frame.groupby(['inner']).mean()
        expected = frame.groupby(pd.Grouper(level='inner')).mean()
        assert_frame_equal(result, expected)

        # Index name
        result = frame.groupby('inner').mean()
        expected = frame.groupby(pd.Grouper(level='inner')).mean()
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('levels', [
        'inner', 'outer', 'B',
        ['inner'], ['outer'], ['B'],
        ['inner', 'outer'], ['inner', 'outer', 'B']
    ])
    def test_grouper_index_level_as_string_series(self, levels):
        s = build_test_series()

        # Compute expected result
        if isinstance(levels, list):
            groupers = [pd.Grouper(level=lv) for lv in levels]
        else:
            groupers = pd.Grouper(level=levels)

        expected = s.groupby(groupers).mean()

        # Compute and check result
        result = s.groupby(levels).mean()
        assert_series_equal(result, expected)

    @pytest.mark.parametrize('frame', [build_df_multi(), build_df_single()])
    def test_grouper_column_index_level_precedence(self, frame):
        # GH 5677, when a string passed as the `by` parameter
        # matches a column and an index level the column takes
        # precedence

        # Add 'inner' column to frame
        # (frame already has an 'inner' index)
        frame['inner'] = [1, 1, 1, 1, 1, 1]

        # Group by single key
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frame.groupby('inner').mean()

        with tm.assert_produces_warning(False):
            expected = frame.groupby(pd.Grouper(key='inner')).mean()

        assert_frame_equal(result, expected)

        with tm.assert_produces_warning(False):
            not_expected = frame.groupby(pd.Grouper(level='inner')).mean()

        assert not result.index.equals(not_expected.index)

        # Group by single key list
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frame.groupby(['inner']).mean()

        with tm.assert_produces_warning(False):
            expected = frame.groupby([pd.Grouper(key='inner')]).mean()

        assert_frame_equal(result, expected)

        with tm.assert_produces_warning(False):
            not_expected = frame.groupby(pd.Grouper(level='inner')).mean()

        assert not result.index.equals(not_expected.index)

        # Group by two keys ('B', 'inner')
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frame.groupby(['B', 'inner']).mean()

        with tm.assert_produces_warning(False):
            expected = frame.groupby(['B',
                                      pd.Grouper(key='inner')]).mean()

        assert_frame_equal(result, expected)

        with tm.assert_produces_warning(False):
            not_expected = frame.groupby(['B',
                                          pd.Grouper(level='inner')
                                          ]).mean()

        assert not result.index.equals(not_expected.index)

        # Group by two keys ('inner', 'B')
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frame.groupby(['inner', 'B']).mean()

        with tm.assert_produces_warning(False):
            expected = frame.groupby([pd.Grouper(key='inner'),
                                      'B']).mean()

        assert_frame_equal(result, expected)

        with tm.assert_produces_warning(False):
            not_expected = frame.groupby([pd.Grouper(level='inner'),
                                          'B']).mean()
        assert not result.index.equals(not_expected.index)
