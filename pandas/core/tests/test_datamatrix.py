from pandas.core.daterange import DateRange
from pandas.core.datetools import bday
from pandas.core.index import Index, NULL_INDEX
from pandas.core.matrix import DataMatrix
from pandas.core.series import Series
from pandas.lib.tseries import map_indices
from copy import deepcopy
from datetime import datetime
from numpy import isnan, array, NaN, alltrue
from numpy import random
from random import choice
import numpy as np
import pickle
import string
import unittest

def rands(n):
    return ''.join([choice(string.letters + string.digits) for i in range(n)])
def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)

#-------------------------------------------------------------------------------
# DataMatrix test cases

class TestDataMatrix(unittest.TestCase):
    def setUp(self):
        index1 = DateRange(datetime(2008,4,22), periods=50)
        index2 = DateRange(datetime(2008,4,29), periods=50)
        index3 = DateRange(datetime(2008,4,28), periods=50)
        ts1 = Series(random.random(50), index=index1)
        ts2 = Series(random.random(50), index=index2)
        ts3 = Series(random.random(50), index=index3)
        ts4 = Series(random.random(50), index=index1)
        data = {'col1' : ts1,'col2' : ts2,'col3' : ts3, 'col4' : ts4}
        self.frame = DataMatrix(data=data, index=index3)
        self.ts1 = ts1
        self.ts2 = ts2
        self.ts3 = ts3
        self.ts4 = ts4

    def testReindex(self):
        newFrame = self.frame.reindex(self.ts1.index)
        for col in newFrame.cols():
            for idx, val in newFrame[col].iteritems():
                if idx in self.frame.index:
                    if isnan(val):
                        self.assert_(isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(isnan(val))
        for col, series in newFrame.iteritems():
            self.assert_(equalContents(series.index, newFrame.index))    
        emptyFrame = self.frame.reindex(Index([]))
        self.assert_(len(emptyFrame.index) == 0)

        nonContigFrame = self.frame.reindex(self.ts1.index[::2])
        for col in nonContigFrame.cols():
            for idx, val in nonContigFrame[col].iteritems():
                if idx in self.frame.index:
                    if isnan(val):
                        self.assert_(isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(isnan(val))
        for col, series in nonContigFrame.iteritems():
            self.assert_(equalContents(series.index, nonContigFrame.index))    

    def testShift(self):
        shiftedFrame = self.frame.shift(5)
        for i, idx in enumerate(shiftedFrame.index):
            self.assert_(idx-5*bday == self.frame.index[i])
        series = shiftedFrame['col1']
        for i, idx in enumerate(series.index):
            self.assert_(idx-5*bday == self.frame.index[i])        

    def testOperators(self):
        garbage = random.random(4)
        colSeries = Series(garbage, index=array(self.frame.cols()))
        idSum = self.frame + self.frame
        seriesSum = self.frame + colSeries
        for col, series in idSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] * 2
                if not isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(isnan(origVal))
        for col, series in seriesSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] + colSeries[col]
                if not isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(isnan(origVal))

    def testSlice(self):
        """Slicing NOT intended for production code"""
        slice = self.frame[:20]
        self.assertEqual(20, len(slice.index))
        for col, series in slice.iteritems():
            self.assertEqual(20, len(series.index))
            self.assert_(equalContents(series.index, slice.index))

    def testGetItem(self):
        for key, value in self.frame._series.iteritems():
            self.assert_(self.frame[key] is not None)
        self.assert_('random' not in self.frame)
        
#    def testGetRow(self):
#        rowFrame = self.frame.getRow(self.frame.index[5])
#        idx = rowFrame.index[0]
#        self.assertEquals(idx, self.frame.index[5])
#        for key, values in rowFrame.iteritems():
#            self.assertEquals(self.frame[key][idx], values[0])
#            self.assertEquals(self.frame[key][idx], values[idx])
        
    def testStack(self):
        frameSlice = self.frame.getTS(fromDate=self.frame.index[0], nPeriods=5)
        stacked = frameSlice.stack()
        for idx, value in stacked.iteritems():
            date, col = idx.split(';')
            date = datetime.fromordinal(int(date))
            if isnan(value):
                self.assert_(isnan(frameSlice[col][date]))
            else:
                self.assertEquals(value, frameSlice[col][date])
        
        unstacked = stacked.unstack().T
        
        for i, idx in enumerate(unstacked.index):
            self.assertEquals(idx, frameSlice.index[i])
        for col, series in unstacked.iteritems():
            for idx, value in series.iteritems():
                if isnan(value):
                    self.assert_(isnan(frameSlice[col][idx]))
                else:                
                    self.assertEquals(value, frameSlice[col][idx])

    def testSetItem(self):
        # not sure what else to do here
        series = self.frame['col1']
        self.frame['col5'] = series
        self.assert_('col5' in self.frame)

    def testStatistics(self):
        sumFrame = self.frame.apply(np.sum)
        for col, series in self.frame.iteritems():
            val = sumFrame[col]
            if isnan(val):
                print self.frame[col]
            self.assertEqual(val, series.sum())
        
    def testDelItem(self):
        del self.frame['col1']
        self.assert_('col1' not in self.frame)
        self.assert_('col1' not in self.frame._series)

    def testGetXS(self):
        idx = self.frame.index[5]
        xs = self.frame.getXS(idx)
        for item, value in xs.iteritems():
            if isnan(value):
                self.assert_(isnan(self.frame[item][idx]))
            else:
                self.assertEqual(value, self.frame[item][idx])

    def testGetTS(self):
        frame = self.frame
        tsFrame = frame.getTS(fromDate=frame.index[5], nPeriods=5)
        for i, idx in enumerate(tsFrame.index):
            self.assertEqual(idx, frame.index[5+i])
            for col, series in tsFrame.iteritems():
                self.assertEqual(idx, series.index[i])
        for col, series in frame.iteritems():
            for idx, value in series.iteritems():
                if isnan(value):
                    self.assert_(isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

    def testTranspose(self):
        frame = self.frame
        dft = frame.T
        for idx, series in dft.iteritems():
            for col, value in series.iteritems():
                if isnan(value):
                    self.assert_(isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

    def testAsMatrix(self):
        frame = self.frame
        mat = frame.asMatrix()
        smallerCols = ['col3', 'col1']
        smallerMat = frame.asMatrix(smallerCols)
        frameCols = frame.cols()
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if isnan(value):
                    self.assert_(isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])
        
    def testDeepcopy(self):
        cp = deepcopy(self.frame)
        series = cp['col1']
        series[:] = 10
        for idx, value in series.iteritems():
            self.assertNotEqual(self.frame['col1'][idx], value)

    def testFilterItems(self):
        pass
    
    def testGroupBy(self):
        
        pass

    def testApply(self):
        pass    

    def testSort(self):
        pass    

    def testToCSV(self):
        pass

    def testPickle(self):
        pass

    def testToDictList(self):
        pass
    
    def testDictToDataFrame(self):
        pass

    def testDataFrameToDict(self):
        pass

    def testFromDict(self):
        newFrame = DataMatrix.fromDict(col1=self.ts1, col2 = self.ts2)
        for idx in newFrame.index:
            if idx in self.ts1.index:
                self.assertEqual(newFrame['col1'][idx], self.ts1[idx])
            if idx in self.ts2.index:
                self.assertEqual(newFrame['col2'][idx], self.ts2[idx])
            
    
    def testPreserveReferences(self):        
        pass
    
    def testCleanNaN(self):
        pass

if __name__ == '__main__':
    unittest.main()
