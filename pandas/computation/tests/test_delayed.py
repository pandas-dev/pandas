import pandas as pd
from pandas import X

import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal


class TestDelayedApi(tm.TestCase):
    def test_basics(self):
        df = pd.DataFrame({'a': [1, 2, 3], 'b': [1.5, 2.5, 3.4],
                           'c': ['abc', 'def', 'efg'],
                           'd': pd.to_datetime(['2014-01-01',
                                                '2014-01-02', '2014-01-03'])})
        assert_frame_equal(df[df['a'] > 1], df[X.a > 1])
        assert_frame_equal(df[df['a'] == 1], df[X.a == 1])
        assert_frame_equal(df.assign(e=lambda x: x['b'] + 1),
                           df.assign(e=X.b + 1))
        assert_frame_equal(df.assign(e=lambda x: x['d'].dt.day),
                           df.assign(e=X.d.dt.day))
        assert_frame_equal(df.assign(e=lambda x: x['c'].str.upper()),
                           df.assign(e=X.c.str.upper()))
