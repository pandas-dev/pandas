from copy import deepcopy
from datetime import datetime
import unittest

from numpy.random import randn
import numpy as np

from pandas.core.api import DataMatrix
import pandas.core.tests.test_frame as test_frame
import pandas.core.tests.common as common

#-------------------------------------------------------------------------------
# DataMatrix test cases

class TestDataMatrix(test_frame.TestDataFrame):
    klass = DataMatrix

    def test_more_constructor(self):
        arr = randn(10)
        dm = self.klass(arr, columns=['A'], index=np.arange(10))
        self.assertEqual(dm.values.ndim, 2)

        arr = randn(0)
        dm = self.klass(arr)
        self.assertEqual(dm.values.ndim, 2)
        self.assertEqual(dm.values.ndim, 2)

        # no data specified
        dm = self.klass(columns=['A', 'B'], index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 2))

        dm = self.klass(columns=['A', 'B'])
        self.assertEqual(dm.values.shape, (0, 2))

        dm = self.klass(index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 0))

        # corner, silly
        self.assertRaises(Exception, self.klass, (1, 2, 3))

    def test_setitem_corner(self):
        # corner case
        df = self.klass({'B' : [1., 2., 3.],
                         'C' : ['a', 'b', 'c']},
                        index=np.arange(3))
        del df['B']
        df['B'] = [1., 2., 3.]
        self.assert_('B' in df)
        self.assertEqual(len(df.columns), 1)

        df['A'] = 'beginning'
        df['E'] = 'foo'
        df['D'] = 'bar'
        df[datetime.now()] = 'date'
        df[datetime.now()] = 5.

    def test_more_fromDict(self):
        pass


if __name__ == '__main__':
    unittest.main()
