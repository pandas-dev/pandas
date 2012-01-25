import unittest

import numpy as np

from pandas import DataFrame, Series
from pandas.tools.merge import concat
from pandas.tools.pivot import pivot_table, crosstab
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

        # this is OK
        table = self.data.pivot_table(rows=['AA', 'BB'], margins=True,
                                      aggfunc='mean')

        # no rows
        rtable = self.data.pivot_table(cols=['AA', 'BB'], margins=True,
                                      aggfunc=np.mean)
        self.assert_(isinstance(rtable, Series))
        for item in ['DD', 'EE', 'FF']:
            gmarg = table[item]['All', '']
            self.assertEqual(gmarg, self.data[item].mean())

    def test_pivot_integer_columns(self):
        # caused by upstream bug in unstack
        from pandas.util.compat import product
        import datetime
        import pandas

        d = datetime.date.min
        data = list(product(['foo', 'bar'], ['A', 'B', 'C'], ['x1', 'x2'],
                [d + datetime.timedelta(i) for i in xrange(20)], [1.0]))
        df = pandas.DataFrame(data)
        table = df.pivot_table(values=4, rows=[0,1,3],cols=[2])

        df2 = df.rename(columns=str)
        table2 = df2.pivot_table(values='4', rows=['0','1','3'], cols=['2'])

        tm.assert_frame_equal(table, table2)

class TestCrosstab(unittest.TestCase):

    def setUp(self):
        df = DataFrame({'A' : ['foo', 'foo', 'foo', 'foo',
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

        self.df = df.append(df, ignore_index=True)

    def test_crosstab_single(self):
        df = self.df
        result = crosstab(df['A'], df['C'])
        expected = df.groupby(['A', 'C']).size().unstack()
        tm.assert_frame_equal(result, expected.fillna(0).astype(np.int64))

    def test_crosstab_multiple(self):
        df = self.df

        result = crosstab(df['A'], [df['B'], df['C']])
        expected = df.groupby(['A', 'B', 'C']).size()
        expected = expected.unstack('B').unstack('C').fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

        result = crosstab([df['B'], df['C']], df['A'])
        expected = df.groupby(['B', 'C', 'A']).size()
        expected = expected.unstack('A').fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_ndarray(self):
        a = np.random.randint(0, 5, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 10, size=100)

        df = DataFrame({'a': a, 'b': b, 'c': c})

        result = crosstab(a, [b, c], rownames=['a'], colnames=('b', 'c'))
        expected = crosstab(df['a'], [df['b'], df['c']])
        tm.assert_frame_equal(result, expected)

        result = crosstab([b, c], a, colnames=['a'], rownames=('b', 'c'))
        expected = crosstab([df['b'], df['c']], df['a'])
        tm.assert_frame_equal(result, expected)

        # assign arbitrary names
        result = crosstab(self.df['A'].values, self.df['C'].values)
        self.assertEqual(result.index.name, 'row_0')
        self.assertEqual(result.columns.name, 'col_0')

    def test_crosstab_margins(self):
        a = np.random.randint(0, 7, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 5, size=100)

        df = DataFrame({'a': a, 'b': b, 'c': c})

        result = crosstab(a, [b, c], rownames=['a'], colnames=('b', 'c'),
                          margins=True)

        self.assertEqual(result.index.names, ['a'])
        self.assertEqual(result.columns.names, ['b', 'c'])

        all_cols = result['All', '']
        exp_cols = df.groupby(['a']).size().astype('i8')
        exp_cols = exp_cols.append(Series([len(df)], index=['All']))

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.ix['All']
        exp_rows = df.groupby(['b', 'c']).size().astype('i8')
        exp_rows = exp_rows.append(Series([len(df)], index=[('All', '')]))

        exp_rows = exp_rows.reindex(all_rows.index)
        exp_rows = exp_rows.fillna(0).astype(np.int64)
        tm.assert_series_equal(all_rows, exp_rows)

    def test_crosstab_pass_values(self):
        a = np.random.randint(0, 7, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 5, size=100)
        values = np.random.randn(100)

        table = crosstab([a, b], c, values, aggfunc=np.sum,
                         rownames=['foo', 'bar'], colnames=['baz'])

        df = DataFrame({'foo': a, 'bar': b, 'baz': c, 'values' : values})

        expected = df.pivot_table('values', rows=['foo', 'bar'], cols='baz',
                                  aggfunc=np.sum)
        tm.assert_frame_equal(table, expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--ipdb', '--ipdb-failure'],
                   exit=False)


