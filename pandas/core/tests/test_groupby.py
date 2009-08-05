import unittest

from pandas.core.daterange import DateRange
from pandas.core.index import Index
from pandas.core.groupby import GroupBy
from pandas.core.pytools import rands, groupby
from pandas.core.frame import DataFrame
from pandas.core.matrix import DataMatrix
from pandas.core.series import Series
import pandas.core.datetools as dt
import pandas.lib.tseries as tseries
import numpy as np

# unittest.TestCase

def commonSetUp(self):    
    self.dateRange = DateRange('1/1/2005', periods=250, offset=dt.bday)
    self.stringIndex = Index([rands(8).upper() for x in xrange(250)])

    self.groupId = Series([x[0] for x in self.stringIndex],
                              index=self.stringIndex)
    self.groupDict = dict((k, v) for k, v in self.groupId.iteritems())
    
    self.columnIndex = Index(['A', 'B', 'C', 'D', 'E'])
    
    randMat = np.random.randn(250, 5)
    self.stringMatrix = DataMatrix(randMat, columns=self.columnIndex,
                                  index=self.stringIndex)
    
    self.timeMatrix = DataMatrix(randMat, columns=self.columnIndex,
                                 index=self.dateRange)


class GroupByTestCase(unittest.TestCase):
    setUp = commonSetUp
    
    def testPythonGrouper(self):
        groupFunc = self.groupDict.get

        groups = groupby(self.stringIndex, groupFunc)
        
        setDict = dict((k, set(v)) for k, v in groups.iteritems())
        for idx in self.stringIndex:
            key = groupFunc(idx)
            groupSet = setDict[key]        
            self.assert_(idx in groupSet)
        
    def testCythonGrouper(self):
        pass
        
    def testNaNGrouping(self):
        pass
    
    def testMembership(self):
        pass
        
    def testByColumnName(self):
        pass
        
class TestAggregate(unittest.TestCase):
    setUp = commonSetUp
    
class TestTransform(unittest.TestCase):
    setUp = commonSetUp
