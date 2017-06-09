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

    def test_list_like_indexing(self):
        # GH 16637
        df = pd.DataFrame({'x': range(10)})
        df.index = pd.to_timedelta(range(10), unit='s')

        conditions = [df.index[0], df.index[4:8], df.index[[3, 5]]]
        expected_data = [[20, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                         [0, 1, 2, 3, 20, 20, 20, 20, 8, 9],
                         [0, 1, 2, 20, 4, 20, 6, 7, 8, 9]]
        for cond, data in zip(conditions, expected_data):
            result = df.copy()
            result.loc[cond, 'x'] = 20
            expected = pd.DataFrame(data,
                                    index=pd.to_timedelta(range(10), unit='s'),
                                    columns=['x'])
            tm.assert_frame_equal(expected, result)
