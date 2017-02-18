import pandas as pd
from pandas.util import testing as tm


class TestTimedeltaIndexing(tm.TestCase):

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
                                    columns=['x'])
            tm.assert_frame_equal(expected, result)
