import pytest

import pandas as pd
from pandas.util import testing as tm


class TestTimedeltaIndexing(object):
    def test_boolean_indexing(self):
        # GH 14946
        df = pd.DataFrame({'x': range(10)})
        df.index = pd.to_timedelta(range(10), unit='s')
        conditions = [df['x'] > 3, df['x'] == 3, df['x'] < 3]
        expected_data = [[0, 1, 2, 3, 10, 10, 10, 10, 10, 10],
                         [0, 1, 2, 10, 4, 5, 6, 7, 8, 9],
                         [10, 10, 10, 3, 4, 5, 6, 7, 8, 9]]
        for cond, data in zip(conditions, expected_data):
            result = df.assign(x=df.mask(cond, 10).astype('int64'))
            expected = pd.DataFrame(data,
                                    index=pd.to_timedelta(range(10), unit='s'),
                                    columns=['x'],
                                    dtype='int64')
            tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize(
        "indexer, expected",
        [(0, [20, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
         (slice(4, 8), [0, 1, 2, 3, 20, 20, 20, 20, 8, 9]),
         ([3, 5], [0, 1, 2, 20, 4, 20, 6, 7, 8, 9])])
    def test_list_like_indexing(self, indexer, expected):
        # GH 16637
        df = pd.DataFrame({'x': range(10)}, dtype="int64")
        df.index = pd.to_timedelta(range(10), unit='s')

        df.loc[df.index[indexer], 'x'] = 20

        expected = pd.DataFrame(expected,
                                index=pd.to_timedelta(range(10), unit='s'),
                                columns=['x'],
                                dtype="int64")

        tm.assert_frame_equal(expected, df)

    def test_string_indexing(self):
        # GH 16896
        df = pd.DataFrame({'x': range(3)},
                          index=pd.to_timedelta(range(3), unit='days'))
        expected = df.iloc[0]
        sliced = df.loc['0 days']
        tm.assert_series_equal(sliced, expected)
