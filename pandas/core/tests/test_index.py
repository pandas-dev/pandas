from pandas.core.index import Index
import pandas.core.tests.common as common
import pandas.lib.tseries as tseries
import numpy as np
import os
import pickle
import unittest

class TestIndex(unittest.TestCase):
    def setUp(self):
        self.strIndex = common.makeStringIndex(100)
        self.dateIndex = common.makeStringIndex(100)
        self.intIndex = common.makeIntIndex(100)

    def testSlice(self):
        strSlice = self.strIndex[10:20]
        dateSlice = self.dateIndex[10:20]
        intSlice = self.intIndex[10:20]
        strMap = tseries.map_indices(np.array(strSlice))
        dateMap = tseries.map_indices(np.array(dateSlice))
        intMap = tseries.map_indices(np.array(intSlice))

        common.assert_dict_equal(strSlice.indexMap, strMap)
        common.assert_dict_equal(dateSlice.indexMap, dateMap)
        common.assert_dict_equal(intSlice.indexMap, intMap)

    def testGetItem(self):
        sl = self.strIndex[[1,2,3]]
        for i in sl:
            self.assertEqual(i, sl[sl.indexMap[i]])
        boolIdx = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIdx[5:30:2] = False
        subIndex = self.strIndex[boolIdx]
        strMap = tseries.map_indices(subIndex)
        for key, value in strMap.iteritems():
            self.assert_(subIndex.indexMap[key] == value)

    def testAdd(self):
        firstCat = self.strIndex + self.dateIndex
        secondCat = self.strIndex + self.strIndex
        self.assert_(common.equalContents(np.append(self.strIndex, self.dateIndex), firstCat))
        self.assert_(common.equalContents(secondCat, self.strIndex))
        for key in self.strIndex:
            self.assert_(key in firstCat.indexMap)
            self.assert_(key in secondCat.indexMap)
        for key in self.dateIndex:
            self.assert_(key in firstCat.indexMap)

    def testContains(self):
        self.assert_(self.strIndex[10] in self.strIndex)
        self.assert_(self.dateIndex[10] in self.dateIndex)
        self.assert_(self.intIndex[10] in self.intIndex)
        strSlice = self.strIndex[10:20]
        dateSlice = self.dateIndex[10:20]
        intSlice = self.intIndex[10:20]
        self.assert_(self.strIndex[9] not in strSlice)
        self.assert_(self.dateIndex[9] not in dateSlice)
        self.assert_(self.intIndex[9] not in intSlice)

    def testMutability(self):
        self.assertRaises(Exception, self.strIndex.__setitem__, 5, 0)
        self.assertRaises(Exception, self.strIndex.__setitem__, slice(1,5), 0)

    def testPickle(self):
        f = open('tmp', 'wb')
        pickle.dump(self.strIndex, f)
        f.close()
        f = open('tmp', 'rb')
        unPickled = pickle.load(f)
        f.close()
        os.remove('tmp')
        self.assert_(isinstance(unPickled, Index))
        self.assert_(common.equalContents(unPickled, self.strIndex))
        for k, v in self.strIndex.indexMap.iteritems():
            self.assert_(k in unPickled.indexMap)
            self.assertEqual(unPickled.indexMap[k], v)
