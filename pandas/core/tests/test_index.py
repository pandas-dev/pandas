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

class TestIndex(unittest.TestCase):
    def setUp(self):
        self.strIndex = Index([rands(10) for i in range(50)])
        self.dateIndex = DateRange(datetime(2008,4,22), periods=50, offset=bday)
        self.intIndex = Index(np.arange(50))
        
    def testSlice(self):
        strSlice = self.strIndex[10:20]
        dateSlice = self.dateIndex[10:20]
        intSlice = self.intIndex[10:20]
        strMap = map_indices(array(strSlice))
        dateMap = map_indices(array(dateSlice))
        intMap = map_indices(array(intSlice))
        for key, value in strMap.iteritems():
            self.assert_(strSlice.indexMap[key] == value)
        for key, value in dateMap.iteritems():
            self.assert_(dateSlice.indexMap[key] == value)        
        for key, value in intMap.iteritems():
            self.assert_(intSlice.indexMap[key] == value)        
    
    def testGetItem(self):
        sl = self.strIndex[[1,2,3]]
        for i in sl:
            self.assertEqual(i, sl[sl.indexMap[i]])
        boolIdx = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIdx[5:30:2] = False
        subIndex = self.strIndex[boolIdx]
        strMap = map_indices(subIndex)
        for key, value in strMap.iteritems():
            self.assert_(subIndex.indexMap[key] == value)
        
    def testAdd(self):
        firstCat = self.strIndex + self.dateIndex
        secondCat = self.strIndex + self.strIndex
        self.assert_(equalContents(np.append(self.strIndex, self.dateIndex), firstCat))
        self.assert_(equalContents(secondCat, self.strIndex))
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
        self.assert_(equalContents(unPickled, self.strIndex))
        for k, v in self.strIndex.indexMap.iteritems():
            self.assert_(k in unPickled.indexMap)
            self.assertEqual(unPickled.indexMap[k], v)
