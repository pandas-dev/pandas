import unittest

import numpy as np

from pandas import DataFrame, concat
from pandas.tools.pivot import pivot_table
import pandas.util.testing as tm

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
                               'E' : np.random.randn(11),
                               'F' : np.random.randn(11)})

    def test_pivot_table(self):
        rows = ['A', 'B']
        cols=  'C'
        table = pivot_table(self.data, values='D', rows=rows, cols=cols)

        table2 = self.data.pivot_table(values='D', rows=rows, cols=cols)
        tm.assert_frame_equal(table, table2)

        # this works
        pivot_table(self.data, values='D', rows=rows)

        if len(rows) > 1:
            self.assertEqual(table.index.names, rows)
        else:
            self.assertEqual(table.index.name, rows[0])

        if len(cols) > 1:
            self.assertEqual(table.columns.names, cols)
        else:
            self.assertEqual(table.columns.name, cols[0])

        expected = self.data.groupby(rows + [cols])['D'].agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_table_multiple(self):
        rows = ['A', 'B']
        cols=  'C'
        table = pivot_table(self.data, rows=rows, cols=cols)
        expected = self.data.groupby(rows + [cols]).agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_multi_values(self):
        result = pivot_table(self.data, values=['D', 'E'],
                             rows='A', cols=['B', 'C'], fill_value=0)
        expected = pivot_table(self.data.drop(['F'], axis=1),
                               rows='A', cols=['B', 'C'], fill_value=0)
        tm.assert_frame_equal(result, expected)

    def test_pivot_multi_functions(self):
        f = lambda func: pivot_table(self.data, values=['D', 'E'],
                                     rows=['A', 'B'], cols='C',
                                     aggfunc=func)
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=['mean', 'std'], axis=1)
        tm.assert_frame_equal(result, expected)

        # margins not supported??
        f = lambda func: pivot_table(self.data, values=['D', 'E'],
                                     rows=['A', 'B'], cols='C',
                                     aggfunc=func, margins=True)
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=['mean', 'std'], axis=1)
        tm.assert_frame_equal(result, expected)


    def test_margins(self):
        def _check_output(res, col, rows=['A', 'B'], cols=['C']):
            cmarg = res['All'][:-1]
            exp = self.data.groupby(rows)[col].mean()
            tm.assert_series_equal(cmarg, exp)

            rmarg = res.xs(('All', ''))[:-1]
            exp = self.data.groupby(cols)[col].mean()
            tm.assert_series_equal(rmarg, exp)

            gmarg = res['All']['All', '']
            exp = self.data[col].mean()
            self.assertEqual(gmarg, exp)

        # column specified
        table = self.data.pivot_table('D', rows=['A', 'B'], cols='C',
                                      margins=True, aggfunc=np.mean)
        _check_output(table, 'D')

        # no column specified
        table = self.data.pivot_table(rows=['A', 'B'], cols='C',
                                      margins=True, aggfunc=np.mean)
        for valcol in table.columns.levels[0]:
            _check_output(table[valcol], valcol)

        # no col

        # to help with a buglet
        self.data.columns = [k * 2 for k in self.data.columns]
        table = self.data.pivot_table(rows=['AA', 'BB'], margins=True,
                                      aggfunc=np.mean)
        for valcol in table.columns:
            gmarg = table[valcol]['All', '']
            self.assertEqual(gmarg, self.data[valcol].mean())

        # doesn't quite work yet

        # # no rows
        # table = self.data.pivot_table(cols=['A', 'B'], margins=True,
        #                               aggfunc=np.mean)
        # for valcol in table.columns:
        #     gmarg = table[valcol]['All', '']
        #     self.assertEqual(gmarg, self.data[valcol].mean())

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


