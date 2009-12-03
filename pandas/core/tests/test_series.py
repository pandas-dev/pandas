from copy import deepcopy
import cPickle as pickle
import os
import unittest

import numpy as np

from pandas.core.series import Series
import pandas.core.tests.common as common

#-------------------------------------------------------------------------------
# Series test cases

class TestSeries(unittest.TestCase):
    def setUp(self):
        self.ts = common.makeTimeSeries()
        self.series = common.makeStringSeries()
        self.objSeries = common.makeObjectSeries()

    def test_constructor(self):
        pass

    def test_fromDict(self):
        pass

    def test_fromValue(self):
        pass

    def test_getitem(self):
        idx1 = self.series.index[5]
        idx2 = self.objSeries.index[5]
        self.assertEqual(self.series[idx1], self.series.get(idx1))
        self.assertEqual(self.objSeries[idx2], self.objSeries.get(idx2))

    def test_fancy(self):
        slice1 = self.series[[1,2,3]]
        slice2 = self.objSeries[[1,2,3]]
        self.assertEqual(self.series.index[2], slice1.index[1])
        self.assertEqual(self.objSeries.index[2], slice2.index[1])
        self.assertEqual(self.series[2], slice1[1])
        self.assertEqual(self.objSeries[2], slice2[1])

    def test_setitem(self):
        self.ts[self.ts.index[5]] = np.NaN
        self.ts[[1,2,17]] = np.NaN
        self.ts[6] = np.NaN
        self.assert_(np.isnan(self.ts[6]))
        self.assert_(np.isnan(self.ts[2]))
        self.ts[np.isnan(self.ts)] = 5
        self.assert_(not np.isnan(self.ts[2]))

    def test_setslice(self):
        slice = self.ts[5:20]
        self.assertEqual(len(slice), len(slice.index))
        self.assertEqual(len(slice.index.indexMap), len(slice.index))

    def test_repr(self):
        pass

    def test_iter(self):
        for i, val in enumerate(self.series):
            self.assertEqual(val, self.series[i])
        for i, val in enumerate(self.ts):
            self.assertEqual(val, self.ts[i])
        for idx, val in self.series.iteritems():
            self.assertEqual(val, self.series[idx])
        for idx, val in self.ts.iteritems():
            self.assertEqual(val, self.ts[idx])

    def test_stats(self):
        self.series[5:15] = np.NaN

        s1 = np.array(self.series)
        s1 = s1[-np.isnan(s1)]
        self.assertEquals(np.mean(s1), self.series.mean())
        self.assertEquals(np.std(s1, ddof=1), self.series.std())
        self.assertEquals(np.var(s1, ddof=1), self.series.var())
        self.assertEquals(np.sum(s1), self.series.sum())
        self.assert_(not np.isnan(np.sum(self.series)))
        self.assert_(not np.isnan(np.mean(self.series)))
        self.assert_(not np.isnan(np.std(self.series)))
        self.assert_(not np.isnan(np.var(self.series)))

    def test_keys(self):
        pass

    def test_values(self):
        pass

    def test_iteritems(self):
        pass

    def test_get(self):
        self.assertEqual(self.series[5], self.series.get(self.series.index[5]))

    def test_append(self):
        appendedSeries = self.series.append(self.ts)
        for idx, value in appendedSeries.iteritems():
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.ts.index:
                self.assertEqual(value, self.ts[idx])
            else:
                self.fail("orphaned index!")

    def test_combineFunc(self):
        shiftedSum = self.ts + self.ts.shift(5)
        idSum = self.ts + self.ts
        self.assert_(np.isnan(shiftedSum[0]))
        for idx, val in idSum.iteritems():
            self.assertAlmostEqual(self.ts[idx] + self.ts[idx], val)
        multiplied = self.ts * 5
        for idx, val in multiplied.iteritems():
            self.assertEqual(self.ts[idx] * 5, val)

    def test_combineFirst(self):
        pass

    def test_argsort(self):
        pass

    def test_cumsum(self):
        pass

    def test_cumprod(self):
        pass

    def test_copy(self):
        pass

    def test_corr(self):
        pass

    def test_count(self):
        pass

    def test_median(self):
        pass

    def test_sort(self):
        pass

    def test_order(self):
        pass

    def test_map(self):
        pass

    def test_toCSV(self):
        pass

    def test_toDict(self):
        pass

    def test_cap(self):
        pass

    def test_floor(self):
        pass

    def test_valid(self):
        pass

    def test_truncate(self):
        pass

    def test_meta(self):
        wrapped = Series(self.series)
        self.assert_(common.equalContents(wrapped.index, self.series.index))
        # Ensure new index is not created
        self.assertEquals(id(self.series.index), id(wrapped.index))

    def test_asOf(self):
        self.ts[5:10] = np.NaN
        self.ts[15:20] = np.NaN
        val1 = self.ts.asOf(self.ts.index[7])
        val2 = self.ts.asOf(self.ts.index[19])
        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])

    def test_pickle(self):
        f = open('tmp1', 'wb')
        h = open('tmp3', 'wb')
        pickle.dump(self.series, f)
        pickle.dump(self.ts, h)
        f.close()
        h.close()
        f = open('tmp1', 'rb')
        h = open('tmp3', 'rb')
        unPickledf = pickle.load(f)
        unPickledh = pickle.load(h)
        f.close()
        h.close()
        os.remove('tmp1')
        os.remove('tmp3')
        self.assert_(isinstance(unPickledf, Series))
        self.assert_(isinstance(unPickledh, Series))
        self.assert_(common.equalContents(unPickledf, self.series))
        self.assert_(common.equalContents(unPickledh, self.ts))
        for idx in self.series.index:
            self.assert_(idx in unPickledf.index)
            self.assertEqual(unPickledf[idx], self.series[idx])
        for idx in self.ts.index:
            self.assert_(idx in unPickledh.index)
            self.assertEqual(unPickledh[idx], self.ts[idx])

    def test_reindex(self):
        identity = self.series.reindex(self.series.index)
        self.assertEqual(id(self.series.index), id(identity.index))
        subIndex = self.series.index[10:20]
        subSeries = self.series.reindex(subIndex)
        for idx, val in subSeries.iteritems():
            self.assertEqual(val, self.series[idx])
        subIndex2 = self.ts.index[10:20]
        subTS = self.ts.reindex(subIndex2)
        for idx, val in subTS.iteritems():
            self.assertEqual(val, self.ts[idx])
        crapSeries = self.ts.reindex(subIndex)
        self.assert_(np.isnan(crapSeries).all())

        # This is extremely important for the Cython code to not screw up
        nonContigIndex = self.ts.index[::2]
        subNonContig = self.ts.reindex(nonContigIndex)
        for idx, val in subNonContig.iteritems():
            self.assertEqual(val, self.ts[idx])

    def test_operators(self):
        newSeries = deepcopy(self.series)
        newSeries[5:10] = np.NaN
        newSeries[10:20] = newSeries[10:20] + 1
        newSeries[20:30] = newSeries[20:30] - 1
        eqSeries = (newSeries == self.series)
        ltSeries = (newSeries < self.series)
        gtSeries = (newSeries > self.series)
        self.assertTrue(ltSeries[20])
        self.assertFalse(ltSeries[10])
        self.assertTrue(gtSeries[10])
        self.assertFalse(gtSeries[20])

    def test_preserveRefs(self):
        slice = self.ts[5:10]
        seq = self.ts[[5,10,15]]
        slice[4] = np.NaN
        seq[1] = np.NaN
        self.assertFalse(np.isnan(self.ts[9]))
        self.assertFalse(np.isnan(self.ts[10]))

class TestTimeSeries(unittest.TestCase):
    def setUp(self):
        self.ts = common.makeTimeSeries()
        self.series = common.makeStringSeries()
        self.objSeries = common.makeObjectSeries()

    def test_shift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)
        idxMap = self.ts.index.indexMap
        for k, v in unshifted.iteritems():
            self.assertEqual(self.ts[idxMap[k]], v)

    def test_slice(self):
        numSlice = self.series[10:20]
        numSliceEnd = self.series[-10:]
        objSlice = self.objSeries[10:20]
        self.assert_(self.series.index[9] not in numSlice.index)
        self.assert_(self.objSeries.index[9] not in objSlice.index)
        self.assertEqual(len(numSlice), len(numSlice.index))
        self.assertEqual(self.series[numSlice.index[0]], numSlice[numSlice.index[0]])
        self.assertEqual(numSlice.index[1], self.series.index[11])
        self.assert_(common.equalContents(numSliceEnd, np.array(self.series)[-10:]))

    def test_fill(self):
        pass

    def test_asfreq(self):
        pass

    def test_interpolate(self):
        pass

    def test_weekday(self):
        pass

    def test_diff(self):
        pass

    def test_autocorr(self):
        pass

    def test_firstValid(self):
        pass

    def test_lastValid(self):
        pass

    def test_reindex(self):
        pass

if __name__ == '__main__':
    unittest.main()
