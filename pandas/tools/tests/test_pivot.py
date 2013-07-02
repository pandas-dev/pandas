import unittest

import numpy as np
from numpy.testing import assert_equal

from pandas import DataFrame, Series, Index, MultiIndex
from pandas.tools.merge import concat
from pandas.tools.pivot import pivot_table, crosstab
import pandas.util.testing as tm


class TestPivotTable(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.data = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                     'bar', 'bar', 'bar', 'bar',
                                     'foo', 'foo', 'foo'],
                               'B': ['one', 'one', 'one', 'two',
                                     'one', 'one', 'one', 'two',
                                     'two', 'two', 'one'],
                               'C': ['dull', 'dull', 'shiny', 'dull',
                                     'dull', 'shiny', 'shiny', 'dull',
                                     'shiny', 'shiny', 'shiny'],
                               'D': np.random.randn(11),
                               'E': np.random.randn(11),
                               'F': np.random.randn(11)})

    def test_pivot_table(self):
        rows = ['A', 'B']
        cols = 'C'
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

    def test_pivot_table_nocols(self):
        df = DataFrame({'rows': ['a', 'b', 'c'],
                        'cols': ['x', 'y', 'z'],
                        'values': [1,2,3]})
        rs = df.pivot_table(cols='cols', aggfunc=np.sum)
        xp = df.pivot_table(rows='cols', aggfunc=np.sum).T
        tm.assert_frame_equal(rs, xp)

        rs = df.pivot_table(cols='cols', aggfunc={'values': 'mean'})
        xp = df.pivot_table(rows='cols', aggfunc={'values': 'mean'}).T
        tm.assert_frame_equal(rs, xp)

    def test_pivot_table_dropna(self):
        df = DataFrame({'amount': {0: 60000, 1: 100000, 2: 50000, 3: 30000},
                        'customer': {0: 'A', 1: 'A', 2: 'B', 3: 'C'},
                        'month': {0: 201307, 1: 201309, 2: 201308, 3: 201310},
                        'product': {0: 'a', 1: 'b', 2: 'c', 3: 'd'},
                        'quantity': {0: 2000000, 1: 500000, 2: 1000000, 3: 1000000}})
        pv_col = df.pivot_table('quantity', 'month', ['customer', 'product'], dropna=False)
        pv_ind = df.pivot_table('quantity', ['customer', 'product'], 'month', dropna=False)

        m = MultiIndex.from_tuples([(u'A', u'a'), (u'A', u'b'), (u'A', u'c'), (u'A', u'd'), 
                                   (u'B', u'a'), (u'B', u'b'), (u'B', u'c'), (u'B', u'd'),
                                   (u'C', u'a'), (u'C', u'b'), (u'C', u'c'), (u'C', u'd')])

        assert_equal(pv_col.columns.values, m.values)
        assert_equal(pv_ind.index.values, m.values)


    def test_pass_array(self):
        result = self.data.pivot_table('D', rows=self.data.A, cols=self.data.C)
        expected = self.data.pivot_table('D', rows='A', cols='C')
        tm.assert_frame_equal(result, expected)

    def test_pass_function(self):
        result = self.data.pivot_table('D', rows=lambda x: x // 5,
                                       cols=self.data.C)
        expected = self.data.pivot_table('D', rows=self.data.index // 5,
                                         cols='C')
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_multiple(self):
        rows = ['A', 'B']
        cols = 'C'
        table = pivot_table(self.data, rows=rows, cols=cols)
        expected = self.data.groupby(rows + [cols]).agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_dtypes(self):

        # can convert dtypes
        f = DataFrame({'a' : ['cat', 'bat', 'cat', 'bat'], 'v' : [1,2,3,4], 'i' : ['a','b','a','b']})
        self.assert_(f.dtypes['v'] == 'int64')

        z = pivot_table(f, values='v', rows=['a'], cols=['i'], fill_value=0, aggfunc=np.sum)
        result = z.get_dtype_counts()
        expected = Series(dict(int64 = 2))
        tm.assert_series_equal(result, expected)

        # cannot convert dtypes
        f = DataFrame({'a' : ['cat', 'bat', 'cat', 'bat'], 'v' : [1.5,2.5,3.5,4.5], 'i' : ['a','b','a','b']})
        self.assert_(f.dtypes['v'] == 'float64')

        z = pivot_table(f, values='v', rows=['a'], cols=['i'], fill_value=0, aggfunc=np.mean)
        result = z.get_dtype_counts()
        expected = Series(dict(float64 = 2))
        tm.assert_series_equal(result, expected)

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

    def test_pivot_index_with_nan(self):
        # GH 3588
        nan = np.nan
        df = DataFrame({"a":['R1', 'R2', nan, 'R4'], 'b':["C1", "C2", "C3" , "C4"], "c":[10, 15, nan , 20]})
        result = df.pivot('a','b','c')
        expected = DataFrame([[nan,nan,nan,nan],[nan,10,nan,nan], 
                              [nan,nan,nan,nan],[nan,nan,15,20]],
                             index = Index(['R1','R2',nan,'R4'],name='a'),
                             columns = Index(['C1','C2','C3','C4'],name='b'))
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
        table = df.pivot_table(values=4, rows=[0, 1, 3], cols=[2])

        df2 = df.rename(columns=str)
        table2 = df2.pivot_table(values='4', rows=['0', '1', '3'], cols=['2'])

        tm.assert_frame_equal(table, table2, check_names=False)

    def test_pivot_no_level_overlap(self):
        # GH #1181

        data = DataFrame({'a': ['a', 'a', 'a', 'a', 'b', 'b', 'b', 'b'] * 2,
                          'b': [0, 0, 0, 0, 1, 1, 1, 1] * 2,
                          'c': (['foo'] * 4 + ['bar'] * 4) * 2,
                          'value': np.random.randn(16)})

        table = data.pivot_table('value', rows='a', cols=['b', 'c'])

        grouped = data.groupby(['a', 'b', 'c'])['value'].mean()
        expected = grouped.unstack('b').unstack('c').dropna(axis=1, how='all')
        tm.assert_frame_equal(table, expected)

    def test_pivot_columns_lexsorted(self):
        import datetime
        import numpy as np
        import pandas

        n = 10000

        dtype = np.dtype([
            ("Index", object),
            ("Symbol", object),
            ("Year", int),
            ("Month", int),
            ("Day", int),
            ("Quantity", int),
            ("Price", float),
        ])

        products = np.array([
            ('SP500', 'ADBE'),
            ('SP500', 'NVDA'),
            ('SP500', 'ORCL'),
            ('NDQ100', 'AAPL'),
            ('NDQ100', 'MSFT'),
            ('NDQ100', 'GOOG'),
            ('FTSE', 'DGE.L'),
            ('FTSE', 'TSCO.L'),
            ('FTSE', 'GSK.L'),
        ], dtype=[('Index', object), ('Symbol', object)])
        items = np.empty(n, dtype=dtype)
        iproduct = np.random.randint(0, len(products), n)
        items['Index'] = products['Index'][iproduct]
        items['Symbol'] = products['Symbol'][iproduct]
        dr = pandas.date_range(datetime.date(2000, 1, 1),
                               datetime.date(2010, 12, 31))
        dates = dr[np.random.randint(0, len(dr), n)]
        items['Year'] = dates.year
        items['Month'] = dates.month
        items['Day'] = dates.day
        items['Price'] = np.random.lognormal(4.0, 2.0, n)

        df = DataFrame(items)

        pivoted = df.pivot_table('Price', rows=['Month', 'Day'],
                                 cols=['Index', 'Symbol', 'Year'],
                                 aggfunc='mean')

        self.assert_(pivoted.columns.is_monotonic)

    def test_pivot_complex_aggfunc(self):
        f = {'D': ['std'], 'E': ['sum']}
        expected = self.data.groupby(['A', 'B']).agg(f).unstack('B')
        result = self.data.pivot_table(rows='A', cols='B', aggfunc=f)

        tm.assert_frame_equal(result, expected)


class TestCrosstab(unittest.TestCase):

    def setUp(self):
        df = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                              'bar', 'bar', 'bar', 'bar',
                              'foo', 'foo', 'foo'],
                        'B': ['one', 'one', 'one', 'two',
                              'one', 'one', 'one', 'two',
                              'two', 'two', 'one'],
                        'C': ['dull', 'dull', 'shiny', 'dull',
                              'dull', 'shiny', 'shiny', 'dull',
                              'shiny', 'shiny', 'shiny'],
                        'D': np.random.randn(11),
                        'E': np.random.randn(11),
                        'F': np.random.randn(11)})

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
        expected = expected.unstack(
            'B').unstack('C').fillna(0).astype(np.int64)
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

        df = DataFrame({'foo': a, 'bar': b, 'baz': c, 'values': values})

        expected = df.pivot_table('values', rows=['foo', 'bar'], cols='baz',
                                  aggfunc=np.sum)
        tm.assert_frame_equal(table, expected)

    def test_crosstab_dropna(self):
        # GH 3820
        a = np.array(['foo', 'foo', 'foo', 'bar', 'bar', 'foo', 'foo'], dtype=object)
        b = np.array(['one', 'one', 'two', 'one', 'two', 'two', 'two'], dtype=object)
        c = np.array(['dull', 'dull', 'dull', 'dull', 'dull', 'shiny', 'shiny'], dtype=object)
        res = crosstab(a, [b, c], rownames=['a'], colnames=['b', 'c'], dropna=False)
        m = MultiIndex.from_tuples([('one', 'dull'), ('one', 'shiny'),
                                    ('two', 'dull'), ('two', 'shiny')])
        assert_equal(res.columns.values, m.values)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
