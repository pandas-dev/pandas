from cStringIO import StringIO
from datetime import datetime
import os
import unittest

from numpy import nan
import numpy as np

from pandas import DataFrame
from pandas.io.parsers import read_csv, read_table, ExcelFile
from pandas.util.testing import assert_almost_equal, assert_frame_equal

class TestParsers(unittest.TestCase):

    def setUp(self):
        self.dirpath = curpath()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')

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

    def test_csv_mixed_type(self):
        data = """A,B,C
a,1,2
b,3,4
c,4,5
"""
        df = read_csv(StringIO(data), index_col=None)
        # TODO

    def test_csv_custom_parser(self):
        data = """A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        df = read_csv(StringIO(data),
                      date_parser=lambda x: datetime.strptime(x, '%Y%m%d'))
        expected = read_csv(StringIO(data))
        assert_frame_equal(df, expected)

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
        df = read_csv(self.csv1)
        df2 = read_table(self.csv1, sep=',')
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D']))
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_read_csv_no_index_name(self):
        df = read_csv(self.csv2)
        df2 = read_table(self.csv2, sep=',')
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D', 'E']))
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.ix[:, ['A', 'B', 'C', 'D']].values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_excel_table(self):
        pth = os.path.join(self.dirpath, 'test.xls')
        xls = ExcelFile(pth)
        df = xls.parse('Sheet1')
        df2 = read_csv(self.csv1)
        df3 = xls.parse('Sheet2', skiprows=[1])
        assert_frame_equal(df, df2)
        assert_frame_equal(df3, df2)

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

