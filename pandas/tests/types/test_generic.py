# -*- coding: utf-8 -*-

import nose
import numpy as np
import pandas as pd
import pandas.core.common as com
import pandas.util.testing as tm

_multiprocess_can_split_ = True


class TestABCClasses(tm.TestCase):
    tuples = [[1, 2, 2], ['red', 'blue', 'red']]
    multi_index = pd.MultiIndex.from_arrays(tuples, names=('number', 'color'))
    datetime_index = pd.to_datetime(['2000/1/1', '2010/1/1'])
    timedelta_index = pd.to_timedelta(np.arange(5), unit='s')
    period_index = pd.period_range('2000/1/1', '2010/1/1/', freq='M')
    categorical = pd.Categorical([1, 2, 3], categories=[2, 3, 1])
    categorical_df = pd.DataFrame({"values": [1, 2, 3]}, index=categorical)
    df = pd.DataFrame({'names': ['a', 'b', 'c']}, index=multi_index)
    sparse_series = pd.Series([1, 2, 3]).to_sparse()
    sparse_array = pd.SparseArray(np.random.randn(10))

    def test_abc_types(self):
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), com.ABCIndex)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), com.ABCInt64Index)
        self.assertIsInstance(pd.Float64Index([1, 2, 3]), com.ABCFloat64Index)
        self.assertIsInstance(self.multi_index, com.ABCMultiIndex)
        self.assertIsInstance(self.datetime_index, com.ABCDatetimeIndex)
        self.assertIsInstance(self.timedelta_index, com.ABCTimedeltaIndex)
        self.assertIsInstance(self.period_index, com.ABCPeriodIndex)
        self.assertIsInstance(self.categorical_df.index,
                              com.ABCCategoricalIndex)
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), com.ABCIndexClass)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), com.ABCIndexClass)
        self.assertIsInstance(pd.Series([1, 2, 3]), com.ABCSeries)
        self.assertIsInstance(self.df, com.ABCDataFrame)
        self.assertIsInstance(self.df.to_panel(), com.ABCPanel)
        self.assertIsInstance(self.sparse_series, com.ABCSparseSeries)
        self.assertIsInstance(self.sparse_array, com.ABCSparseArray)
        self.assertIsInstance(self.categorical, com.ABCCategorical)
        self.assertIsInstance(pd.Period('2012', freq='A-DEC'), com.ABCPeriod)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
