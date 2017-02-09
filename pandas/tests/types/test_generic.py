# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas.types import generic as gt


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
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), gt.ABCIndex)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), gt.ABCInt64Index)
        self.assertIsInstance(pd.UInt64Index([1, 2, 3]), gt.ABCUInt64Index)
        self.assertIsInstance(pd.Float64Index([1, 2, 3]), gt.ABCFloat64Index)
        self.assertIsInstance(self.multi_index, gt.ABCMultiIndex)
        self.assertIsInstance(self.datetime_index, gt.ABCDatetimeIndex)
        self.assertIsInstance(self.timedelta_index, gt.ABCTimedeltaIndex)
        self.assertIsInstance(self.period_index, gt.ABCPeriodIndex)
        self.assertIsInstance(self.categorical_df.index,
                              gt.ABCCategoricalIndex)
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), gt.ABCIndexClass)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), gt.ABCIndexClass)
        self.assertIsInstance(pd.Series([1, 2, 3]), gt.ABCSeries)
        self.assertIsInstance(self.df, gt.ABCDataFrame)
        self.assertIsInstance(self.df.to_panel(), gt.ABCPanel)
        self.assertIsInstance(self.sparse_series, gt.ABCSparseSeries)
        self.assertIsInstance(self.sparse_array, gt.ABCSparseArray)
        self.assertIsInstance(self.categorical, gt.ABCCategorical)
        self.assertIsInstance(pd.Period('2012', freq='A-DEC'), gt.ABCPeriod)
