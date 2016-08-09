import numpy as np
import pandas as pd
import pandas.util.testing as tm


class TestSparseDataFrameMultitype(tm.TestCase):
    def setUp(self):
        super(TestSparseDataFrameMultitype, self).setUp()
        self.string_series = pd.SparseSeries(['a', 'b', 'c'])
        self.int_series = pd.SparseSeries([1, 2, 3])
        self.float_series = pd.SparseSeries([1.1, 1.2, 1.3])
        self.object_series = pd.SparseSeries([[], {}, set()])
        self.sdf = pd.SparseDataFrame({
            'string': self.string_series,
            'int': self.int_series,
            'float': self.float_series,
            'object': self.object_series,
        })
        self.cols = ['string', 'int', 'float', 'object']
        self.sdf = self.sdf[self.cols]

    def test_basic_dtypes(self):
        for _, row in self.sdf.iterrows():
            self.assertEqual(row.dtype, object)
        tm.assert_sp_series_equal(self.sdf['string'], self.string_series,
                                  check_names=False)
        tm.assert_sp_series_equal(self.sdf['int'], self.int_series,
                                  check_names=False)
        tm.assert_sp_series_equal(self.sdf['float'], self.float_series,
                                  check_names=False)
        tm.assert_sp_series_equal(self.sdf['object'], self.object_series,
                                  check_names=False)

    def test_indexing_single(self):
        tm.assert_sp_series_equal(self.sdf.iloc[0],
                                  pd.SparseSeries(['a', 1, 1.1, []],
                                                  index=self.cols),
                                  check_names=False)
        tm.assert_sp_series_equal(self.sdf.iloc[1],
                                  pd.SparseSeries(['b', 2, 1.2, {}],
                                                  index=self.cols),
                                  check_names=False)
        tm.assert_sp_series_equal(self.sdf.iloc[2],
                                  pd.SparseSeries(['c', 3, 1.3, set()],
                                                  index=self.cols),
                                  check_names=False)

    def test_indexing_multiple(self):
        tm.assert_sp_frame_equal(self.sdf, self.sdf[:])
        tm.assert_sp_frame_equal(self.sdf, self.sdf.loc[:])
        tm.assert_sp_frame_equal(self.sdf.iloc[[1, 2]],
                                 pd.SparseDataFrame({
                                     'string': self.string_series.iloc[[1, 2]],
                                     'int': self.int_series.iloc[[1, 2]],
                                     'float': self.float_series.iloc[[1, 2]],
                                     'object': self.object_series.iloc[[1, 2]]
                                 }, index=[1, 2])[self.cols])
        tm.assert_sp_frame_equal(self.sdf[['int', 'string']],
                                 pd.SparseDataFrame({
                                     'int': self.int_series,
                                     'string': self.string_series,
                                 }))


class TestSparseSeriesMultitype(tm.TestCase):
    def setUp(self):
        super(TestSparseSeriesMultitype, self).setUp()
        self.index = ['string', 'int', 'float', 'object']
        self.ss = pd.SparseSeries(['a', 1, 1.1, []],
                                  index=self.index)

    def test_indexing_single(self):
        for i, idx in enumerate(self.index):
            self.assertEqual(self.ss.iloc[i], self.ss[idx])
            self.assertEqual(type(self.ss.iloc[i]),
                             type(self.ss[idx]))
        self.assertEqual(self.ss['string'], 'a')
        self.assertEqual(self.ss['int'], 1)
        self.assertEqual(self.ss['float'], 1.1)
        self.assertEqual(self.ss['object'], [])

    def test_indexing_multiple(self):
        tm.assert_sp_series_equal(self.ss.loc[['string', 'int']],
                                  pd.SparseSeries(['a', 1],
                                                  index=['string', 'int']))
        tm.assert_sp_series_equal(self.ss.loc[['string', 'object']],
                                  pd.SparseSeries(['a', []],
                                                  index=['string', 'object']))
