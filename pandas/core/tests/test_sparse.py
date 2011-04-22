from unittest import TestCase
import operator

from numpy import nan
import numpy as np

from pandas.util.testing import assert_almost_equal, assert_series_equal
from numpy.testing import assert_equal

from pandas.core.sparse import (IntIndex, BlockIndex, SparseSeries)
import pandas.core.sparse as sparsem

"""
Testing TODO


"""
def _test_data1():
    # nan-based
    arr = np.arange(20, dtype=float)
    index = np.arange(20)
    arr[:2] = nan
    arr[5:10] = nan
    arr[-3:] = nan

    return arr, index

def _test_data2():
    # nan-based
    arr = np.arange(15, dtype=float)
    index = np.arange(15)
    arr[7:12] = nan
    arr[-1:] = nan
    return arr, index

def _test_data1_zero():
    # zero-based
    arr, index = _test_data1()
    arr[np.isnan(arr)] = 0
    return arr, index

def _test_data2_zero():
    # zero-based
    arr, index = _test_data2()
    arr[np.isnan(arr)] = 0
    return arr, index

def assert_sp_series_equal(a, b):
    assert_equal(a.sp_values, b.sp_values)
    assert(a.sp_index.equals(b.sp_index))
    if np.isnan(a.fill_value):
        assert(np.isnan(b.fill_value))
    else:
        assert(a.fill_value == b.fill_value)

def assert_sp_frame_equal(left, right):
    for col, series in left.iteritems():
        assert(col in right)
        assert_sp_series_equal(series, right[col])

    for col in right:
        assert(col in left)

class TestSparseSeries(TestCase):

    def setUp(self):
        arr, index = _test_data1()
        self.bseries = SparseSeries(arr, index=index, kind='block')
        self.iseries = SparseSeries(arr, index=index, kind='integer')

        arr, index = _test_data2()
        self.bseries2 = SparseSeries(arr, index=index, kind='block')
        self.iseries2 = SparseSeries(arr, index=index, kind='integer')

        arr, index = _test_data1_zero()
        self.zbseries = SparseSeries(arr, index=index, kind='block',
                                     fill_value=0)
        self.ziseries = SparseSeries(arr, index=index, kind='integer',
                                     fill_value=0)

        arr, index = _test_data2_zero()
        self.zbseries2 = SparseSeries(arr, index=index, kind='block',
                                      fill_value=0)
        self.ziseries2 = SparseSeries(arr, index=index, kind='integer',
                                      fill_value=0)

    def test_to_dense(self):
        arr, index = _test_data1()
        series = self.bseries.to_dense()
        assert_equal(series, arr)

        series = self.iseries.to_dense()
        assert_equal(series, arr)

    def test_to_sparse_series(self):
        series = self.bseries.to_dense()
        bseries = sparsem.to_sparse_series(series, kind='block')
        iseries = sparsem.to_sparse_series(series, kind='integer')
        assert_sp_series_equal(bseries, self.bseries)
        assert_sp_series_equal(iseries, self.iseries)

        # non-NaN fill value
        series = self.zbseries.to_dense()
        zbseries = sparsem.to_sparse_series(series, kind='block',
                                            fill_value=0)
        ziseries = sparsem.to_sparse_series(series, kind='integer',
                                            fill_value=0)
        assert_sp_series_equal(zbseries, self.zbseries)
        assert_sp_series_equal(ziseries, self.ziseries)

    def test_constructor(self):
        # test setup guys
        self.assert_(np.isnan(self.bseries.fill_value))
        self.assert_(isinstance(self.bseries.sp_index, BlockIndex))
        self.assert_(np.isnan(self.iseries.fill_value))
        self.assert_(isinstance(self.iseries.sp_index, IntIndex))

        self.assertEquals(self.zbseries.fill_value, 0)
        assert_equal(self.zbseries.values, self.bseries.to_dense().fillna(0))

        # pass SparseSeries

        # pass Series
        series = self.bseries.to_dense()
        bseries2 = SparseSeries(series)
        assert_equal(self.bseries.sp_values, bseries2.sp_values)

        # pass dict

    def test_constructor_nonnan(self):
        arr = [0, 0, 0, nan, nan]
        sp_series = SparseSeries(arr, fill_value=0)
        assert_equal(sp_series.values, arr)

    def test_getitem(self):
        pass

    def test_setitem(self):
        pass

    def test_slice(self):
        pass

    def test_setslice(self):
        pass

    def test_serialize(self):
        pass

    def test_operators(self):
        def _check_op(a, b, op):
            sp_result = op(a, b)
            adense = a.to_dense() if isinstance(a, SparseSeries) else a
            bdense = b.to_dense() if isinstance(b, SparseSeries) else b
            dense_result = op(adense, bdense)
            assert_almost_equal(sp_result.to_dense(), dense_result)

        def check(a, b):
            _check_op(a, b, operator.add)
            _check_op(a, b, operator.sub)
            _check_op(a, b, operator.div)
            _check_op(a, b, operator.mul)

            _check_op(a, b, lambda x, y: operator.add(y, x))
            _check_op(a, b, lambda x, y: operator.sub(y, x))
            _check_op(a, b, lambda x, y: operator.div(y, x))
            _check_op(a, b, lambda x, y: operator.mul(y, x))

            # NaN ** 0 = 1 in C?
            # _check_op(a, b, operator.pow)
            # _check_op(a, b, lambda x, y: operator.pow(y, x))

        check(self.bseries, self.bseries)
        check(self.iseries, self.iseries)
        check(self.bseries, self.iseries)

        check(self.bseries, self.bseries2)
        check(self.bseries, self.iseries2)
        check(self.iseries, self.iseries2)

        # scalar value
        check(self.bseries, 5)

        # zero-based
        check(self.zbseries, self.zbseries * 2)
        check(self.zbseries, self.zbseries2)
        check(self.ziseries, self.ziseries2)

    def test_reindex(self):
        pass

    def test_repr(self):
        pass

    def test_iter(self):
        pass

    def test_truncate(self):
        pass

    def test_fillna(self):
        pass

    def test_groupby(self):
        pass

class TestSparseTimeSeries(TestCase):
    pass


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

