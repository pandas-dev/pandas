from cStringIO import StringIO
from datetime import datetime
import os
import unittest

from numpy import nan
import numpy as np

from pandas import DataFrame
from pandas.io.parsers import read_csv, read_table, ExcelFile
from pandas.util.testing import assert_almost_equal, assert_frame_equal

class TestReadTable(unittest.TestCase):

    def setUp(self):
        self.dirpath = curpath()

    def test_read_csv(self):
        pass

    def test_custom_na_values(self):
        data = """A,B,C
ignore,this,row
1,NA,3
-1.#IND,5,baz
7,8,NaN
"""
        expected = [[1., nan, 3],
                    [nan, 5, nan],
                    [7, 8, nan]]

        df = read_csv(StringIO(data), index_col=None, na_values=['baz'],
                      skiprows=[1])
        assert_almost_equal(df.values, expected)

        df2 = read_table(StringIO(data), sep=',', index_col=None,
                         na_values=['baz'], skiprows=[1])
        assert_almost_equal(df2.values, expected)

    def test_unnamed_columns(self):
        data = """A,B,C,,
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        expected = [[1,2,3,4,5.],
                    [6,7,8,9,10],
                    [11,12,13,14,15]]
        df = read_table(StringIO(data), sep=',', index_col=None)
        assert_almost_equal(df.values, expected)
        self.assert_(np.array_equal(df.columns,
                                    ['A', 'B', 'C', 'Unnamed: 3',
                                     'Unnamed: 4']))

    def test_duplicate_columns(self):
        data = """A,A,B,B,B
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df = read_table(StringIO(data), sep=',', index_col=None)
        self.assert_(np.array_equal(df.columns,
                                    ['A', 'A.1', 'B', 'B.1', 'B.2']))

    def test_no_header(self):
        data = """1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df = read_table(StringIO(data), sep=',', index_col=None,
                        header=None)
        names = ['foo', 'bar', 'baz', 'quux', 'panda']
        df2 = read_table(StringIO(data), sep=',', index_col=None,
                        header=None, names=names)
        expected = [[1,2,3,4,5.],
                    [6,7,8,9,10],
                    [11,12,13,14,15]]
        assert_almost_equal(df.values, expected)
        self.assert_(np.array_equal(df.columns,
                                    ['X.1', 'X.2', 'X.3', 'X.4', 'X.5']))
        self.assert_(np.array_equal(df2.columns, names))

    def test_read_csv_dataframe(self):
        pth = os.path.join(self.dirpath, 'test1.csv')
        df = read_csv(pth)
        df2 = read_table(pth, sep=',')
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D']))
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_read_csv_no_index_name(self):
        pth = os.path.join(self.dirpath, 'test2.csv')
        df = read_csv(pth)
        df2 = read_table(pth, sep=',')
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D']))
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.values.dtype == np.float64)
        assert_frame_equal(df, df2)

class TestExcelFile(unittest.TestCase):
    pass

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

