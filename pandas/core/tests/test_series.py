from pandas.core.daterange import DateRange
from pandas.core.datetools import bday
from pandas.core.index import Index
from pandas.core.series import Series
from pandas.lib.tseries import map_indices
from copy import deepcopy
from datetime import datetime
from numpy import isnan, array, NaN, alltrue
from numpy import random
from random import choice
import numpy as np
import os
import pickle
import string
import sys
import unittest

def rands(n):
    return ''.join([choice(string.letters + string.digits) for i in range(n)])
def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)

#-------------------------------------------------------------------------------
# Series test cases

class TestSeries(unittest.TestCase):
    def setUp(self):
        index = Index([rands(10) for i in range(50)])
        dateIndex = DateRange(datetime(2008,4,22), periods=50, offset=bday)
        self.ts = Series(random.random(50), index=dateIndex)
        self.series = Series(random.random(50), index=index)
        self.objSeries = Series(dateIndex, index=index)
        #self.plainSeries = Series(random.random(50))
    
    def testSlice(self):
        numSlice = self.series[10:20]
        numSliceEnd = self.series[-10:]
        objSlice = self.objSeries[10:20]
        self.assert_(self.series.index[9] not in numSlice.index)
        self.assert_(self.objSeries.index[9] not in objSlice.index)
        self.assertEqual(len(numSlice), len(numSlice.index))
        self.assertEqual(self.series[numSlice.index[0]], numSlice[numSlice.index[0]])
        self.assertEqual(numSlice.index[1], self.series.index[11])
        self.assert_(equalContents(numSliceEnd, array(self.series)[-10:]))
    
    def testGet(self):
        self.assertEqual(self.series[5], self.series.get(self.series.index[5]))
    
    def testGetItem(self):
        idx1 = self.series.index[5]
        idx2 = self.objSeries.index[5]
        self.assertEqual(self.series[idx1], self.series.get(idx1))
        self.assertEqual(self.objSeries[idx2], self.objSeries.get(idx2))
    
    def testSetItem(self):
        self.ts[self.ts.index[5]] = NaN
        self.ts[[1,2,27]] = NaN
        self.ts[6] = NaN
        self.assert_(isnan(self.ts[6]))
        self.assert_(isnan(self.ts[2]))
        self.ts[isnan(self.ts)] = 5
        self.assert_(not isnan(self.ts[2]))
        
    def testSetSlice(self):
        slice = self.ts[5:20]
        self.assertEqual(len(slice), len(slice.index))
        self.assertEqual(len(slice.index.indexMap), len(slice.index))
        
    def testGetSequence(self):
        slice1 = self.series[[1,2,3]]
        slice2 = self.objSeries[[1,2,3]]
        self.assertEqual(self.series.index[2], slice1.index[1])
        self.assertEqual(self.objSeries.index[2], slice2.index[1])
        self.assertEqual(self.series[2], slice1[1])
        self.assertEqual(self.objSeries[2], slice2[1])

        #slice3 = self.plainSeries[[1,2,3]]
        #self.assertEqual(self.plainSeries[2], slice3[1])
        
    def testMeta(self):
        wrapped = Series(self.series)
        self.assert_(equalContents(wrapped.index, self.series.index))
        # Ensure new index is not created
        self.assertEquals(id(self.series.index), id(wrapped.index))
    
    def testStatistics(self):
        self.series[5:15] = NaN
        
        s1 = array(self.series)
        s1 = s1[-isnan(s1)]
        self.assertEquals(np.mean(s1), self.series.mean())
        self.assertEquals(np.std(s1, ddof=1), self.series.std())
        self.assertEquals(np.var(s1, ddof=1), self.series.var())
        self.assertEquals(np.sum(s1), self.series.sum())
        self.assert_(not isnan(np.sum(self.series)))
        self.assert_(not isnan(np.mean(self.series)))
        self.assert_(not isnan(np.std(self.series)))
        self.assert_(not isnan(np.var(self.series)))

        #self.plainSeries[20:25] = NaN
        #s2 = array(self.plainSeries)
        #s2 = s2[-isnan(s2)]
        #self.assertEquals(np.mean(s2), self.plainSeries.mean())
        #self.assertEquals(np.std(s2, ddof=1), self.plainSeries.std())
        #self.assertEquals(np.var(s2, ddof=1), self.plainSeries.var())
        #self.assertEquals(np.sum(s2), self.plainSeries.sum())

        #self.assert_(isnan(self.series.mean(removeNA=False)))
        #self.assert_(isnan(self.series.std(removeNA=False)))
        #self.assert_(isnan(self.series.var(removeNA=False)))
        #self.assert_(isnan(self.series.sum(removeNA=False)))
        
    def testPickle(self):
        f = open('tmp1', 'wb')
        #g = open('tmp2', 'wb')
        h = open('tmp3', 'wb')
        pickle.dump(self.series, f)
        #pickle.dump(self.plainSeries, g)
        pickle.dump(self.ts, h)
        f.close()
        #g.close()
        h.close()
        f = open('tmp1', 'rb')
        #g = open('tmp2', 'rb')
        h = open('tmp3', 'rb')
        unPickledf = pickle.load(f)
        #unPickledg = pickle.load(g)
        unPickledh = pickle.load(h)
        f.close()
        #g.close()
        h.close()
        os.remove('tmp1')
        #os.remove('tmp2')
        os.remove('tmp3')
        self.assert_(isinstance(unPickledf, Series))
        #self.assert_(isinstance(unPickledg, Series))
        self.assert_(isinstance(unPickledh, Series))
        self.assert_(equalContents(unPickledf, self.series))
        #self.assert_(equalContents(unPickledg, self.plainSeries))
        self.assert_(equalContents(unPickledh, self.ts))
        for idx in self.series.index:
            self.assert_(idx in unPickledf.index)
            self.assertEqual(unPickledf[idx], self.series[idx])    
        for idx in self.ts.index:
            self.assert_(idx in unPickledh.index)
            self.assertEqual(unPickledh[idx], self.ts[idx])    
            
    def testIter(self):
        for i, val in enumerate(self.series):
            self.assertEqual(val, self.series[i])
        for i, val in enumerate(self.ts):
            self.assertEqual(val, self.ts[i])
        for idx, val in self.series.iteritems():
            self.assertEqual(val, self.series[idx])
        for idx, val in self.ts.iteritems():
            self.assertEqual(val, self.ts[idx])            
        #self.assertRaises(Exception, self.plainSeries.iteritems)
        
    def testFromIndex(self):
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
        self.assert_(alltrue(isnan(crapSeries)))
        
        # This is extremely important for the Cython code to not screw up
        nonContigIndex = self.ts.index[::2]
        subNonContig = self.ts.reindex(nonContigIndex)
        for idx, val in subNonContig.iteritems():
            self.assertEqual(val, self.ts[idx])            
    
    def testCombineFunc(self):
        shiftedSum = self.ts + self.ts.shift(5)
        idSum = self.ts + self.ts
        self.assert_(isnan(shiftedSum[0]))
        for idx, val in idSum.iteritems():
            self.assertAlmostEqual(self.ts[idx] + self.ts[idx], val)
        multiplied = self.ts * 5
        for idx, val in multiplied.iteritems():
            self.assertEqual(self.ts[idx] * 5, val)

    def testOperators(self):
        newSeries = deepcopy(self.series)
        newSeries[5:10] = NaN
        newSeries[10:20] = newSeries[10:20] + 1
        newSeries[20:30] = newSeries[20:30] - 1
        eqSeries = (newSeries == self.series)
        ltSeries = (newSeries < self.series)
        gtSeries = (newSeries > self.series)
        self.assertTrue(ltSeries[20])
        self.assertFalse(ltSeries[10])
        self.assertTrue(gtSeries[10])
        self.assertFalse(gtSeries[20])
        
        
    def testShift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)
        #self.assert_(equalContents(self.ts.index, unshifted.index))
        idxMap = self.ts.index.indexMap
        for k, v in unshifted.iteritems():
            self.assertEqual(self.ts[idxMap[k]], v)

    def testAsOf(self):
        self.ts[5:10] = NaN
        self.ts[15:20] = NaN
        val1 = self.ts.asOf(self.ts.index[7])
        val2 = self.ts.asOf(self.ts.index[19])
        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])
        
    def testPreserveReferences(self):
        slice = self.ts[5:10]
        seq = self.ts[[5,10,15]]
        slice[4] = NaN
        seq[1] = NaN
        self.assertFalse(isnan(self.ts[9]))
        self.assertFalse(isnan(self.ts[10]))
        
    def testAppend(self):
        appendedSeries = self.series.append(self.ts)
        for idx, value in appendedSeries.iteritems():
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.ts.index:
                self.assertEqual(value, self.ts[idx])
            else:
                self.fail("orphaned index!")

if __name__ == '__main__':
    unittest.main()
