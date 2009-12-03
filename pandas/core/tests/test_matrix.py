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
        pass

    def test_more_fromDict(self):
        pass


if __name__ == '__main__':
    unittest.main()
