import unittest

import numpy as np

from pandas import DataFrame
from pandas.tools.pivot import pivot_table
from pandas.util.testing import assert_frame_equal

class TestPivotTable(unittest.TestCase):

    def setUp(self):
        self.data = DataFrame({'A' : ['foo', 'foo', 'foo', 'foo',
                                      'bar', 'bar', 'bar', 'bar',
                                      'foo', 'foo', 'foo'],
                               'B' : ['one', 'one', 'one', 'two',
                                      'one', 'one', 'one', 'two',
                                      'two', 'two', 'one'],
                               'C' : ['dull', 'dull', 'shiny', 'dull',
                                      'dull', 'shiny', 'shiny', 'dull',
                                      'shiny', 'shiny', 'shiny'],
                               'D' : np.random.randn(11),
                               'E' : np.random.randn(11)})

    def test_pivot_table(self):
        xby = ['A', 'B']
        yby=  ['C']
        table = pivot_table(self.data, values='D', xby=xby, yby=yby)

        if len(xby) > 1:
            self.assertEqual(table.index.names, xby)
        else:
            self.assertEqual(table.index.name, xby[0])

        if len(yby) > 1:
            self.assertEqual(table.columns.names, yby)
        else:
            self.assertEqual(table.columns.name, yby[0])

        expected = self.data.groupby(xby + yby)['D'].agg(np.mean).unstack()
        assert_frame_equal(table, expected)

    def test_pivot_table_multiple(self):
        xby = ['A', 'B']
        yby=  ['C']
        table = pivot_table(self.data, xby=xby, yby=yby)
        expected = self.data.groupby(xby + yby).agg(np.mean).unstack()
        assert_frame_equal(table, expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


