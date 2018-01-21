# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd
import pandas.util.testing as tm


class TestPeriodFrameArithmetic(object):

    def test_ops_frame_period(self):
        # GH 13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        assert df['A'].dtype == object
        assert df['B'].dtype == object

        p = pd.Period('2015-03', freq='M')
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2, 1], dtype=object),
                            'B': np.array([14, 13], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        assert df2['A'].dtype == object
        assert df2['B'].dtype == object

        exp = pd.DataFrame({'A': np.array([4, 4], dtype=object),
                            'B': np.array([16, 16], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -exp)
