from datetime import datetime, date, timedelta

import numpy as np
from numpy.testing import assert_equal

import pandas as pd
from pandas import DataFrame, Series, Index, MultiIndex, Grouper
from pandas.tools.merge import concat
from pandas.tools.pivot import pivot_table, crosstab
from pandas.compat import range, u, product
import pandas.util.testing as tm


class TestPivotTable(tm.TestCase):

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
        index = ['A', 'B']
        columns = 'C'
        table = pivot_table(self.data, values='D',
                            index=index, columns=columns)

        table2 = self.data.pivot_table(
            values='D', index=index, columns=columns)
        tm.assert_frame_equal(table, table2)

        # this works
        pivot_table(self.data, values='D', index=index)

        if len(index) > 1:
            self.assertEqual(table.index.names, tuple(index))
        else:
            self.assertEqual(table.index.name, index[0])

        if len(columns) > 1:
            self.assertEqual(table.columns.names, columns)
        else:
            self.assertEqual(table.columns.name, columns[0])

        expected = self.data.groupby(
            index + [columns])['D'].agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_table_nocols(self):
        df = DataFrame({'rows': ['a', 'b', 'c'],
                        'cols': ['x', 'y', 'z'],
                        'values': [1, 2, 3]})
        rs = df.pivot_table(columns='cols', aggfunc=np.sum)
        xp = df.pivot_table(index='cols', aggfunc=np.sum).T
        tm.assert_frame_equal(rs, xp)

        rs = df.pivot_table(columns='cols', aggfunc={'values': 'mean'})
        xp = df.pivot_table(index='cols', aggfunc={'values': 'mean'}).T
        tm.assert_frame_equal(rs, xp)

    def test_pivot_table_dropna(self):
        df = DataFrame({'amount': {0: 60000, 1: 100000, 2: 50000, 3: 30000},
                        'customer': {0: 'A', 1: 'A', 2: 'B', 3: 'C'},
                        'month': {0: 201307, 1: 201309, 2: 201308, 3: 201310},
                        'product': {0: 'a', 1: 'b', 2: 'c', 3: 'd'},
                        'quantity': {0: 2000000, 1: 500000,
                                     2: 1000000, 3: 1000000}})
        pv_col = df.pivot_table('quantity', 'month', [
                                'customer', 'product'], dropna=False)
        pv_ind = df.pivot_table(
            'quantity', ['customer', 'product'], 'month', dropna=False)

        m = MultiIndex.from_tuples([(u('A'), u('a')),
                                    (u('A'), u('b')),
                                    (u('A'), u('c')),
                                    (u('A'), u('d')),
                                    (u('B'), u('a')),
                                    (u('B'), u('b')),
                                    (u('B'), u('c')),
                                    (u('B'), u('d')),
                                    (u('C'), u('a')),
                                    (u('C'), u('b')),
                                    (u('C'), u('c')),
                                    (u('C'), u('d'))])

        assert_equal(pv_col.columns.values, m.values)
        assert_equal(pv_ind.index.values, m.values)

    def test_pass_array(self):
        result = self.data.pivot_table(
            'D', index=self.data.A, columns=self.data.C)
        expected = self.data.pivot_table('D', index='A', columns='C')
        tm.assert_frame_equal(result, expected)

    def test_pass_function(self):
        result = self.data.pivot_table('D', index=lambda x: x // 5,
                                       columns=self.data.C)
        expected = self.data.pivot_table('D', index=self.data.index // 5,
                                         columns='C')
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_multiple(self):
        index = ['A', 'B']
        columns = 'C'
        table = pivot_table(self.data, index=index, columns=columns)
        expected = self.data.groupby(index + [columns]).agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_dtypes(self):

        # can convert dtypes
        f = DataFrame({'a': ['cat', 'bat', 'cat', 'bat'], 'v': [
                      1, 2, 3, 4], 'i': ['a', 'b', 'a', 'b']})
        self.assertEqual(f.dtypes['v'], 'int64')

        z = pivot_table(f, values='v', index=['a'], columns=[
                        'i'], fill_value=0, aggfunc=np.sum)
        result = z.get_dtype_counts()
        expected = Series(dict(int64=2))
        tm.assert_series_equal(result, expected)

        # cannot convert dtypes
        f = DataFrame({'a': ['cat', 'bat', 'cat', 'bat'], 'v': [
                      1.5, 2.5, 3.5, 4.5], 'i': ['a', 'b', 'a', 'b']})
        self.assertEqual(f.dtypes['v'], 'float64')

        z = pivot_table(f, values='v', index=['a'], columns=[
                        'i'], fill_value=0, aggfunc=np.mean)
        result = z.get_dtype_counts()
        expected = Series(dict(float64=2))
        tm.assert_series_equal(result, expected)

    def test_pivot_multi_values(self):
        result = pivot_table(self.data, values=['D', 'E'],
                             index='A', columns=['B', 'C'], fill_value=0)
        expected = pivot_table(self.data.drop(['F'], axis=1),
                               index='A', columns=['B', 'C'], fill_value=0)
        tm.assert_frame_equal(result, expected)

    def test_pivot_multi_functions(self):
        f = lambda func: pivot_table(self.data, values=['D', 'E'],
                                     index=['A', 'B'], columns='C',
                                     aggfunc=func)
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=['mean', 'std'], axis=1)
        tm.assert_frame_equal(result, expected)

        # margins not supported??
        f = lambda func: pivot_table(self.data, values=['D', 'E'],
                                     index=['A', 'B'], columns='C',
                                     aggfunc=func, margins=True)
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=['mean', 'std'], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_pivot_index_with_nan(self):
        # GH 3588
        nan = np.nan
        df = DataFrame({'a': ['R1', 'R2', nan, 'R4'],
                        'b': ['C1', 'C2', 'C3', 'C4'],
                        'c': [10, 15, 17, 20]})
        result = df.pivot('a', 'b', 'c')
        expected = DataFrame([[nan, nan, 17, nan], [10, nan, nan, nan],
                              [nan, 15, nan, nan], [nan, nan, nan, 20]],
                             index=Index([nan, 'R1', 'R2', 'R4'], name='a'),
                             columns=Index(['C1', 'C2', 'C3', 'C4'], name='b'))
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df.pivot('b', 'a', 'c'), expected.T)

        # GH9491
        df = DataFrame({'a': pd.date_range('2014-02-01', periods=6, freq='D'),
                        'c': 100 + np.arange(6)})
        df['b'] = df['a'] - pd.Timestamp('2014-02-02')
        df.loc[1, 'a'] = df.loc[3, 'a'] = nan
        df.loc[1, 'b'] = df.loc[4, 'b'] = nan

        pv = df.pivot('a', 'b', 'c')
        self.assertEqual(pv.notnull().values.sum(), len(df))

        for _, row in df.iterrows():
            self.assertEqual(pv.loc[row['a'], row['b']], row['c'])

        tm.assert_frame_equal(df.pivot('b', 'a', 'c'), pv.T)

    def test_pivot_with_tz(self):
        # GH 5878
        df = DataFrame({'dt1': [datetime(2013, 1, 1, 9, 0),
                                datetime(2013, 1, 2, 9, 0),
                                datetime(2013, 1, 1, 9, 0),
                                datetime(2013, 1, 2, 9, 0)],
                        'dt2': [datetime(2014, 1, 1, 9, 0),
                                datetime(2014, 1, 1, 9, 0),
                                datetime(2014, 1, 2, 9, 0),
                                datetime(2014, 1, 2, 9, 0)],
                        'data1': np.arange(4, dtype='int64'),
                        'data2': np.arange(4, dtype='int64')})

        df['dt1'] = df['dt1'].apply(lambda d: pd.Timestamp(d, tz='US/Pacific'))
        df['dt2'] = df['dt2'].apply(lambda d: pd.Timestamp(d, tz='Asia/Tokyo'))

        exp_col1 = Index(['data1', 'data1', 'data2', 'data2'])
        exp_col2 = pd.DatetimeIndex(['2014/01/01 09:00',
                                     '2014/01/02 09:00'] * 2,
                                    name='dt2', tz='Asia/Tokyo')
        exp_col = pd.MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame([[0, 2, 0, 2], [1, 3, 1, 3]],
                             index=pd.DatetimeIndex(['2013/01/01 09:00',
                                                     '2013/01/02 09:00'],
                                                    name='dt1',
                                                    tz='US/Pacific'),
                             columns=exp_col)

        pv = df.pivot(index='dt1', columns='dt2')
        tm.assert_frame_equal(pv, expected)

        expected = DataFrame([[0, 2], [1, 3]],
                             index=pd.DatetimeIndex(['2013/01/01 09:00',
                                                     '2013/01/02 09:00'],
                                                    name='dt1',
                                                    tz='US/Pacific'),
                             columns=pd.DatetimeIndex(['2014/01/01 09:00',
                                                       '2014/01/02 09:00'],
                                                      name='dt2',
                                                      tz='Asia/Tokyo'))

        pv = df.pivot(index='dt1', columns='dt2', values='data1')
        tm.assert_frame_equal(pv, expected)

    def test_pivot_periods(self):
        df = DataFrame({'p1': [pd.Period('2013-01-01', 'D'),
                               pd.Period('2013-01-02', 'D'),
                               pd.Period('2013-01-01', 'D'),
                               pd.Period('2013-01-02', 'D')],
                        'p2': [pd.Period('2013-01', 'M'),
                               pd.Period('2013-01', 'M'),
                               pd.Period('2013-02', 'M'),
                               pd.Period('2013-02', 'M')],
                        'data1': np.arange(4, dtype='int64'),
                        'data2': np.arange(4, dtype='int64')})

        exp_col1 = Index(['data1', 'data1', 'data2', 'data2'])
        exp_col2 = pd.PeriodIndex(['2013-01', '2013-02'] * 2,
                                  name='p2', freq='M')
        exp_col = pd.MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame([[0, 2, 0, 2], [1, 3, 1, 3]],
                             index=pd.PeriodIndex(['2013-01-01', '2013-01-02'],
                                                  name='p1', freq='D'),
                             columns=exp_col)

        pv = df.pivot(index='p1', columns='p2')
        tm.assert_frame_equal(pv, expected)

        expected = DataFrame([[0, 2], [1, 3]],
                             index=pd.PeriodIndex(['2013-01-01', '2013-01-02'],
                                                  name='p1', freq='D'),
                             columns=pd.PeriodIndex(['2013-01', '2013-02'],
                                                    name='p2', freq='M'))

        pv = df.pivot(index='p1', columns='p2', values='data1')
        tm.assert_frame_equal(pv, expected)

    def test_margins(self):
        def _check_output(result, values_col, index=['A', 'B'],
                          columns=['C'],
                          margins_col='All'):
            col_margins = result.ix[:-1, margins_col]
            expected_col_margins = self.data.groupby(index)[values_col].mean()
            tm.assert_series_equal(col_margins, expected_col_margins,
                                   check_names=False)
            self.assertEqual(col_margins.name, margins_col)

            result = result.sortlevel()
            index_margins = result.ix[(margins_col, '')].iloc[:-1]
            expected_ix_margins = self.data.groupby(columns)[values_col].mean()
            tm.assert_series_equal(index_margins, expected_ix_margins,
                                   check_names=False)
            self.assertEqual(index_margins.name, (margins_col, ''))

            grand_total_margins = result.loc[(margins_col, ''), margins_col]
            expected_total_margins = self.data[values_col].mean()
            self.assertEqual(grand_total_margins, expected_total_margins)

        # column specified
        result = self.data.pivot_table(values='D', index=['A', 'B'],
                                       columns='C',
                                       margins=True, aggfunc=np.mean)
        _check_output(result, 'D')

        # Set a different margins_name (not 'All')
        result = self.data.pivot_table(values='D', index=['A', 'B'],
                                       columns='C',
                                       margins=True, aggfunc=np.mean,
                                       margins_name='Totals')
        _check_output(result, 'D', margins_col='Totals')

        # no column specified
        table = self.data.pivot_table(index=['A', 'B'], columns='C',
                                      margins=True, aggfunc=np.mean)
        for value_col in table.columns.levels[0]:
            _check_output(table[value_col], value_col)

        # no col

        # to help with a buglet
        self.data.columns = [k * 2 for k in self.data.columns]
        table = self.data.pivot_table(index=['AA', 'BB'], margins=True,
                                      aggfunc=np.mean)
        for value_col in table.columns:
            totals = table.loc[('All', ''), value_col]
            self.assertEqual(totals, self.data[value_col].mean())

        # no rows
        rtable = self.data.pivot_table(columns=['AA', 'BB'], margins=True,
                                       aggfunc=np.mean)
        tm.assertIsInstance(rtable, Series)

        table = self.data.pivot_table(index=['AA', 'BB'], margins=True,
                                      aggfunc='mean')
        for item in ['DD', 'EE', 'FF']:
            totals = table.loc[('All', ''), item]
            self.assertEqual(totals, self.data[item].mean())

        # issue number #8349: pivot_table with margins and dictionary aggfunc
        data = [
            {'JOB': 'Worker', 'NAME': 'Bob', 'YEAR': 2013,
             'MONTH': 12, 'DAYS': 3, 'SALARY': 17},
            {'JOB': 'Employ', 'NAME':
             'Mary', 'YEAR': 2013, 'MONTH': 12, 'DAYS': 5, 'SALARY': 23},
            {'JOB': 'Worker', 'NAME': 'Bob', 'YEAR': 2014,
             'MONTH': 1, 'DAYS': 10, 'SALARY': 100},
            {'JOB': 'Worker', 'NAME': 'Bob', 'YEAR': 2014,
             'MONTH': 1, 'DAYS': 11, 'SALARY': 110},
            {'JOB': 'Employ', 'NAME': 'Mary', 'YEAR': 2014,
             'MONTH': 1, 'DAYS': 15, 'SALARY': 200},
            {'JOB': 'Worker', 'NAME': 'Bob', 'YEAR': 2014,
             'MONTH': 2, 'DAYS': 8, 'SALARY': 80},
            {'JOB': 'Employ', 'NAME': 'Mary', 'YEAR': 2014,
             'MONTH': 2, 'DAYS': 5, 'SALARY': 190},
        ]

        df = DataFrame(data)

        df = df.set_index(['JOB', 'NAME', 'YEAR', 'MONTH'], drop=False,
                          append=False)

        result = df.pivot_table(index=['JOB', 'NAME'],
                                columns=['YEAR', 'MONTH'],
                                values=['DAYS', 'SALARY'],
                                aggfunc={'DAYS': 'mean', 'SALARY': 'sum'},
                                margins=True)

        expected = df.pivot_table(index=['JOB', 'NAME'],
                                  columns=['YEAR', 'MONTH'], values=['DAYS'],
                                  aggfunc='mean', margins=True)

        tm.assert_frame_equal(result['DAYS'], expected['DAYS'])

        expected = df.pivot_table(index=['JOB', 'NAME'],
                                  columns=['YEAR', 'MONTH'], values=['SALARY'],
                                  aggfunc='sum', margins=True)

        tm.assert_frame_equal(result['SALARY'], expected['SALARY'])

    def test_pivot_integer_columns(self):
        # caused by upstream bug in unstack

        d = date.min
        data = list(product(['foo', 'bar'], ['A', 'B', 'C'], ['x1', 'x2'],
                            [d + timedelta(i)
                             for i in range(20)], [1.0]))
        df = DataFrame(data)
        table = df.pivot_table(values=4, index=[0, 1, 3], columns=[2])

        df2 = df.rename(columns=str)
        table2 = df2.pivot_table(
            values='4', index=['0', '1', '3'], columns=['2'])

        tm.assert_frame_equal(table, table2, check_names=False)

    def test_pivot_no_level_overlap(self):
        # GH #1181

        data = DataFrame({'a': ['a', 'a', 'a', 'a', 'b', 'b', 'b', 'b'] * 2,
                          'b': [0, 0, 0, 0, 1, 1, 1, 1] * 2,
                          'c': (['foo'] * 4 + ['bar'] * 4) * 2,
                          'value': np.random.randn(16)})

        table = data.pivot_table('value', index='a', columns=['b', 'c'])

        grouped = data.groupby(['a', 'b', 'c'])['value'].mean()
        expected = grouped.unstack('b').unstack('c').dropna(axis=1, how='all')
        tm.assert_frame_equal(table, expected)

    def test_pivot_columns_lexsorted(self):

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
        dr = pd.date_range(date(2000, 1, 1),
                           date(2010, 12, 31))
        dates = dr[np.random.randint(0, len(dr), n)]
        items['Year'] = dates.year
        items['Month'] = dates.month
        items['Day'] = dates.day
        items['Price'] = np.random.lognormal(4.0, 2.0, n)

        df = DataFrame(items)

        pivoted = df.pivot_table('Price', index=['Month', 'Day'],
                                 columns=['Index', 'Symbol', 'Year'],
                                 aggfunc='mean')

        self.assertTrue(pivoted.columns.is_monotonic)

    def test_pivot_complex_aggfunc(self):
        f = {'D': ['std'], 'E': ['sum']}
        expected = self.data.groupby(['A', 'B']).agg(f).unstack('B')
        result = self.data.pivot_table(index='A', columns='B', aggfunc=f)

        tm.assert_frame_equal(result, expected)

    def test_margins_no_values_no_cols(self):
        # Regression test on pivot table: no values or cols passed.
        result = self.data[['A', 'B']].pivot_table(
            index=['A', 'B'], aggfunc=len, margins=True)
        result_list = result.tolist()
        self.assertEqual(sum(result_list[:-1]), result_list[-1])

    def test_margins_no_values_two_rows(self):
        # Regression test on pivot table: no values passed but rows are a
        # multi-index
        result = self.data[['A', 'B', 'C']].pivot_table(
            index=['A', 'B'], columns='C', aggfunc=len, margins=True)
        self.assertEqual(result.All.tolist(), [3.0, 1.0, 4.0, 3.0, 11.0])

    def test_margins_no_values_one_row_one_col(self):
        # Regression test on pivot table: no values passed but row and col
        # defined
        result = self.data[['A', 'B']].pivot_table(
            index='A', columns='B', aggfunc=len, margins=True)
        self.assertEqual(result.All.tolist(), [4.0, 7.0, 11.0])

    def test_margins_no_values_two_row_two_cols(self):
        # Regression test on pivot table: no values passed but rows and cols
        # are multi-indexed
        self.data['D'] = ['a', 'b', 'c', 'd',
                          'e', 'f', 'g', 'h', 'i', 'j', 'k']
        result = self.data[['A', 'B', 'C', 'D']].pivot_table(
            index=['A', 'B'], columns=['C', 'D'], aggfunc=len, margins=True)
        self.assertEqual(result.All.tolist(), [3.0, 1.0, 4.0, 3.0, 11.0])

    def test_pivot_table_with_margins_set_margin_name(self):
        # GH 3335
        for margin_name in ['foo', 'one', 666, None, ['a', 'b']]:
            with self.assertRaises(ValueError):
                # multi-index index
                pivot_table(self.data, values='D', index=['A', 'B'],
                            columns=['C'], margins=True,
                            margins_name=margin_name)
            with self.assertRaises(ValueError):
                # multi-index column
                pivot_table(self.data, values='D', index=['C'],
                            columns=['A', 'B'], margins=True,
                            margins_name=margin_name)
            with self.assertRaises(ValueError):
                # non-multi-index index/column
                pivot_table(self.data, values='D', index=['A'],
                            columns=['B'], margins=True,
                            margins_name=margin_name)

    def test_pivot_timegrouper(self):
        df = DataFrame({
            'Branch': 'A A A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Carl Joe Joe Joe Carl'.split(),
            'Quantity': [1, 3, 5, 1, 8, 1, 9, 3],
            'Date': [datetime(2013, 1, 1),
                     datetime(2013, 1, 1),
                     datetime(2013, 10, 1),
                     datetime(2013, 10, 2),
                     datetime(2013, 10, 1),
                     datetime(2013, 10, 2),
                     datetime(2013, 12, 2),
                     datetime(2013, 12, 2), ]}).set_index('Date')

        expected = DataFrame(np.array([10, 18, 3], dtype='int64')
                             .reshape(1, 3),
                             index=[datetime(2013, 12, 31)],
                             columns='Carl Joe Mark'.split())
        expected.index.name = 'Date'
        expected.columns.name = 'Buyer'

        result = pivot_table(df, index=Grouper(freq='A'), columns='Buyer',
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index='Buyer', columns=Grouper(freq='A'),
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

        expected = DataFrame(np.array([1, np.nan, 3, 9, 18, np.nan])
                             .reshape(2, 3),
                             index=[datetime(2013, 1, 1),
                                    datetime(2013, 7, 1)],
                             columns='Carl Joe Mark'.split())
        expected.index.name = 'Date'
        expected.columns.name = 'Buyer'

        result = pivot_table(df, index=Grouper(freq='6MS'), columns='Buyer',
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index='Buyer', columns=Grouper(freq='6MS'),
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

        # passing the name
        df = df.reset_index()
        result = pivot_table(df, index=Grouper(freq='6MS', key='Date'),
                             columns='Buyer',
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index='Buyer',
                             columns=Grouper(freq='6MS', key='Date'),
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

        self.assertRaises(KeyError, lambda: pivot_table(
            df, index=Grouper(freq='6MS', key='foo'),
            columns='Buyer', values='Quantity', aggfunc=np.sum))
        self.assertRaises(KeyError, lambda: pivot_table(
            df, index='Buyer',
            columns=Grouper(freq='6MS', key='foo'),
            values='Quantity', aggfunc=np.sum))

        # passing the level
        df = df.set_index('Date')
        result = pivot_table(df, index=Grouper(freq='6MS', level='Date'),
                             columns='Buyer', values='Quantity',
                             aggfunc=np.sum)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index='Buyer',
                             columns=Grouper(freq='6MS', level='Date'),
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

        self.assertRaises(ValueError, lambda: pivot_table(
            df, index=Grouper(freq='6MS', level='foo'),
            columns='Buyer', values='Quantity', aggfunc=np.sum))
        self.assertRaises(ValueError, lambda: pivot_table(
            df, index='Buyer',
            columns=Grouper(freq='6MS', level='foo'),
            values='Quantity', aggfunc=np.sum))

        # double grouper
        df = DataFrame({
            'Branch': 'A A A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Carl Joe Joe Joe Carl'.split(),
            'Quantity': [1, 3, 5, 1, 8, 1, 9, 3],
            'Date': [datetime(2013, 11, 1, 13, 0), datetime(2013, 9, 1, 13, 5),
                     datetime(2013, 10, 1, 20, 0),
                     datetime(2013, 10, 2, 10, 0),
                     datetime(2013, 11, 1, 20, 0),
                     datetime(2013, 10, 2, 10, 0),
                     datetime(2013, 10, 2, 12, 0),
                     datetime(2013, 12, 5, 14, 0)],
            'PayDay': [datetime(2013, 10, 4, 0, 0),
                       datetime(2013, 10, 15, 13, 5),
                       datetime(2013, 9, 5, 20, 0),
                       datetime(2013, 11, 2, 10, 0),
                       datetime(2013, 10, 7, 20, 0),
                       datetime(2013, 9, 5, 10, 0),
                       datetime(2013, 12, 30, 12, 0),
                       datetime(2013, 11, 20, 14, 0), ]})

        result = pivot_table(df, index=Grouper(freq='M', key='Date'),
                             columns=Grouper(freq='M', key='PayDay'),
                             values='Quantity', aggfunc=np.sum)
        expected = DataFrame(np.array([np.nan, 3, np.nan, np.nan,
                                       6, np.nan, 1, 9,
                                       np.nan, 9, np.nan, np.nan, np.nan,
                                       np.nan, 3, np.nan]).reshape(4, 4),
                             index=[datetime(2013, 9, 30),
                                    datetime(2013, 10, 31),
                                    datetime(2013, 11, 30),
                                    datetime(2013, 12, 31)],
                             columns=[datetime(2013, 9, 30),
                                      datetime(2013, 10, 31),
                                      datetime(2013, 11, 30),
                                      datetime(2013, 12, 31)])
        expected.index.name = 'Date'
        expected.columns.name = 'PayDay'

        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index=Grouper(freq='M', key='PayDay'),
                             columns=Grouper(freq='M', key='Date'),
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

        tuples = [(datetime(2013, 9, 30), datetime(2013, 10, 31)),
                  (datetime(2013, 10, 31),
                   datetime(2013, 9, 30)),
                  (datetime(2013, 10, 31),
                   datetime(2013, 11, 30)),
                  (datetime(2013, 10, 31),
                   datetime(2013, 12, 31)),
                  (datetime(2013, 11, 30),
                   datetime(2013, 10, 31)),
                  (datetime(2013, 12, 31), datetime(2013, 11, 30)), ]
        idx = MultiIndex.from_tuples(tuples, names=['Date', 'PayDay'])
        expected = DataFrame(np.array([3, np.nan, 6, np.nan, 1, np.nan,
                                       9, np.nan, 9, np.nan,
                                       np.nan, 3]).reshape(6, 2),
                             index=idx, columns=['A', 'B'])
        expected.columns.name = 'Branch'

        result = pivot_table(
            df, index=[Grouper(freq='M', key='Date'),
                       Grouper(freq='M', key='PayDay')], columns=['Branch'],
            values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index=['Branch'],
                             columns=[Grouper(freq='M', key='Date'),
                                      Grouper(freq='M', key='PayDay')],
                             values='Quantity', aggfunc=np.sum)
        tm.assert_frame_equal(result, expected.T)

    def test_pivot_datetime_tz(self):
        dates1 = ['2011-07-19 07:00:00', '2011-07-19 08:00:00',
                  '2011-07-19 09:00:00',
                  '2011-07-19 07:00:00', '2011-07-19 08:00:00',
                  '2011-07-19 09:00:00']
        dates2 = ['2013-01-01 15:00:00', '2013-01-01 15:00:00',
                  '2013-01-01 15:00:00',
                  '2013-02-01 15:00:00', '2013-02-01 15:00:00',
                  '2013-02-01 15:00:00']
        df = DataFrame({'label': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'dt1': dates1, 'dt2': dates2,
                        'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2] * 3})
        df['dt1'] = df['dt1'].apply(lambda d: pd.Timestamp(d, tz='US/Pacific'))
        df['dt2'] = df['dt2'].apply(lambda d: pd.Timestamp(d, tz='Asia/Tokyo'))

        exp_idx = pd.DatetimeIndex(['2011-07-19 07:00:00',
                                    '2011-07-19 08:00:00',
                                    '2011-07-19 09:00:00'],
                                   tz='US/Pacific', name='dt1')
        exp_col1 = Index(['value1', 'value1'])
        exp_col2 = Index(['a', 'b'], name='label')
        exp_col = MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame([[0, 3], [1, 4], [2, 5]],
                             index=exp_idx, columns=exp_col)
        result = pivot_table(df, index=['dt1'], columns=[
                             'label'], values=['value1'])
        tm.assert_frame_equal(result, expected)

        exp_col1 = Index(['sum', 'sum', 'sum', 'sum',
                          'mean', 'mean', 'mean', 'mean'])
        exp_col2 = Index(['value1', 'value1', 'value2', 'value2'] * 2)
        exp_col3 = pd.DatetimeIndex(['2013-01-01 15:00:00',
                                     '2013-02-01 15:00:00'] * 4,
                                    tz='Asia/Tokyo', name='dt2')
        exp_col = MultiIndex.from_arrays([exp_col1, exp_col2, exp_col3])
        expected = DataFrame(np.array([[0, 3, 1, 2, 0, 3, 1, 2],
                                       [1, 4, 2, 1, 1, 4, 2, 1],
                                       [2, 5, 1, 2, 2, 5, 1, 2]],
                                      dtype='int64'),
                             index=exp_idx,
                             columns=exp_col)

        result = pivot_table(df, index=['dt1'], columns=['dt2'],
                             values=['value1', 'value2'],
                             aggfunc=[np.sum, np.mean])
        tm.assert_frame_equal(result, expected)

    def test_pivot_dtaccessor(self):
        # GH 8103
        dates1 = ['2011-07-19 07:00:00', '2011-07-19 08:00:00',
                  '2011-07-19 09:00:00',
                  '2011-07-19 07:00:00', '2011-07-19 08:00:00',
                  '2011-07-19 09:00:00']
        dates2 = ['2013-01-01 15:00:00', '2013-01-01 15:00:00',
                  '2013-01-01 15:00:00',
                  '2013-02-01 15:00:00', '2013-02-01 15:00:00',
                  '2013-02-01 15:00:00']
        df = DataFrame({'label': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'dt1': dates1, 'dt2': dates2,
                        'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2] * 3})
        df['dt1'] = df['dt1'].apply(lambda d: pd.Timestamp(d))
        df['dt2'] = df['dt2'].apply(lambda d: pd.Timestamp(d))

        result = pivot_table(df, index='label', columns=df['dt1'].dt.hour,
                             values='value1')

        exp_idx = Index(['a', 'b'], name='label')
        expected = DataFrame({7: [0, 3], 8: [1, 4], 9: [2, 5]},
                             index=exp_idx,
                             columns=Index([7, 8, 9], name='dt1'))
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index=df['dt2'].dt.month,
                             columns=df['dt1'].dt.hour,
                             values='value1')

        expected = DataFrame({7: [0, 3], 8: [1, 4], 9: [2, 5]},
                             index=Index([1, 2], name='dt2'),
                             columns=Index([7, 8, 9], name='dt1'))
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index=df['dt2'].dt.year.values,
                             columns=[df['dt1'].dt.hour, df['dt2'].dt.month],
                             values='value1')

        exp_col = MultiIndex.from_arrays(
            [[7, 7, 8, 8, 9, 9], [1, 2] * 3], names=['dt1', 'dt2'])
        expected = DataFrame(np.array([[0, 3, 1, 4, 2, 5]], dtype='int64'),
                             index=[2013], columns=exp_col)
        tm.assert_frame_equal(result, expected)

        result = pivot_table(df, index=np.array(['X', 'X', 'X',
                                                 'X', 'Y', 'Y']),
                             columns=[df['dt1'].dt.hour, df['dt2'].dt.month],
                             values='value1')
        expected = DataFrame(np.array([[0, 3, 1, np.nan, 2, np.nan],
                                       [np.nan, np.nan, np.nan,
                                        4, np.nan, 5]]),
                             index=['X', 'Y'], columns=exp_col)
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_with_iterator_values(self):
        # GH 12017
        aggs = {'D': 'sum', 'E': 'mean'}

        pivot_values_list = pd.pivot_table(
            self.data, index=['A'], values=list(aggs.keys()), aggfunc=aggs,
        )

        pivot_values_keys = pd.pivot_table(
            self.data, index=['A'], values=aggs.keys(), aggfunc=aggs,
        )
        tm.assert_frame_equal(pivot_values_keys, pivot_values_list)

        agg_values_gen = (value for value in aggs.keys())
        pivot_values_gen = pd.pivot_table(
            self.data, index=['A'], values=agg_values_gen, aggfunc=aggs,
        )
        tm.assert_frame_equal(pivot_values_gen, pivot_values_list)


class TestCrosstab(tm.TestCase):

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

        self.assertEqual(result.index.names, ('a',))
        self.assertEqual(result.columns.names, ['b', 'c'])

        all_cols = result['All', '']
        exp_cols = df.groupby(['a']).size().astype('i8')
        exp_cols = exp_cols.append(Series([len(df)], index=['All']))
        exp_cols.name = ('All', '')

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.ix['All']
        exp_rows = df.groupby(['b', 'c']).size().astype('i8')
        exp_rows = exp_rows.append(Series([len(df)], index=[('All', '')]))
        exp_rows.name = 'All'

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

        expected = df.pivot_table('values', index=['foo', 'bar'],
                                  columns='baz', aggfunc=np.sum)
        tm.assert_frame_equal(table, expected)

    def test_crosstab_dropna(self):
        # GH 3820
        a = np.array(['foo', 'foo', 'foo', 'bar',
                      'bar', 'foo', 'foo'], dtype=object)
        b = np.array(['one', 'one', 'two', 'one',
                      'two', 'two', 'two'], dtype=object)
        c = np.array(['dull', 'dull', 'dull', 'dull',
                      'dull', 'shiny', 'shiny'], dtype=object)
        res = pd.crosstab(a, [b, c], rownames=['a'],
                          colnames=['b', 'c'], dropna=False)
        m = MultiIndex.from_tuples([('one', 'dull'), ('one', 'shiny'),
                                    ('two', 'dull'), ('two', 'shiny')])
        assert_equal(res.columns.values, m.values)

    def test_categorical_margins(self):
        # GH 10989
        df = pd.DataFrame({'x': np.arange(8),
                           'y': np.arange(8) // 4,
                           'z': np.arange(8) % 2})

        expected = pd.DataFrame([[1.0, 2.0, 1.5], [5, 6, 5.5], [3, 4, 3.5]])
        expected.index = Index([0, 1, 'All'], name='y')
        expected.columns = Index([0, 1, 'All'], name='z')

        data = df.copy()
        table = data.pivot_table('x', 'y', 'z', margins=True)
        tm.assert_frame_equal(table, expected)

        data = df.copy()
        data.y = data.y.astype('category')
        data.z = data.z.astype('category')
        table = data.pivot_table('x', 'y', 'z', margins=True)
        tm.assert_frame_equal(table, expected)

    def test_crosstab_no_overlap(self):
        # GS 10291

        s1 = pd.Series([1, 2, 3], index=[1, 2, 3])
        s2 = pd.Series([4, 5, 6], index=[4, 5, 6])

        actual = crosstab(s1, s2)
        expected = pd.DataFrame()

        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna(self):
        # GH 12577
        # pivot_table counts null into margin ('All')
        # when margins=true and dropna=true

        df = pd.DataFrame({'a': [1, 2, 2, 2, 2, np.nan],
                           'b': [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [2, 3, 5]])
        expected.index = Index([1.0, 2.0, 'All'], name='a')
        expected.columns = Index([3, 4, 'All'], name='b')
        tm.assert_frame_equal(actual, expected)

        df = DataFrame({'a': [1, np.nan, np.nan, np.nan, 2, np.nan],
                        'b': [3, np.nan, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = Index([1.0, 2.0, 'All'], name='a')
        expected.columns = Index([3.0, 4.0, 'All'], name='b')
        tm.assert_frame_equal(actual, expected)

        df = DataFrame({'a': [1, np.nan, np.nan, np.nan, np.nan, 2],
                        'b': [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = Index([1.0, 2.0, 'All'], name='a')
        expected.columns = Index([3, 4, 'All'], name='b')
        tm.assert_frame_equal(actual, expected)

        # GH 12642
        # _add_margins raises KeyError: Level None not found
        # when margins=True and dropna=False
        df = pd.DataFrame({'a': [1, 2, 2, 2, 2, np.nan],
                           'b': [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [2, 4, 6]])
        expected.index = Index([1.0, 2.0, 'All'], name='a')
        expected.columns = Index([3, 4, 'All'], name='b')
        tm.assert_frame_equal(actual, expected)

        df = DataFrame({'a': [1, np.nan, np.nan, np.nan, 2, np.nan],
                        'b': [3, np.nan, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 4, 6]])
        expected.index = Index([1.0, 2.0, 'All'], name='a')
        expected.columns = Index([3.0, 4.0, 'All'], name='b')
        tm.assert_frame_equal(actual, expected)

        a = np.array(['foo', 'foo', 'foo', 'bar',
                      'bar', 'foo', 'foo'], dtype=object)
        b = np.array(['one', 'one', 'two', 'one',
                      'two', np.nan, 'two'], dtype=object)
        c = np.array(['dull', 'dull', 'dull', 'dull',
                      'dull', 'shiny', 'shiny'], dtype=object)

        actual = pd.crosstab(a, [b, c], rownames=['a'],
                             colnames=['b', 'c'], margins=True, dropna=False)
        m = MultiIndex.from_arrays([['one', 'one', 'two', 'two', 'All'],
                                    ['dull', 'shiny', 'dull', 'shiny', '']],
                                   names=['b', 'c'])
        expected = DataFrame([[1, 0, 1, 0, 2], [2, 0, 1, 1, 5],
                              [3, 0, 2, 1, 7]], columns=m)
        expected.index = Index(['bar', 'foo', 'All'], name='a')
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab([a, b], c, rownames=['a', 'b'],
                             colnames=['c'], margins=True, dropna=False)
        m = MultiIndex.from_arrays([['bar', 'bar', 'foo', 'foo', 'All'],
                                    ['one', 'two', 'one', 'two', '']],
                                   names=['a', 'b'])
        expected = DataFrame([[1, 0, 1], [1, 0, 1], [2, 0, 2], [1, 1, 2],
                              [5, 2, 7]], index=m)
        expected.columns = Index(['dull', 'shiny', 'All'], name='c')
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab([a, b], c, rownames=['a', 'b'],
                             colnames=['c'], margins=True, dropna=True)
        m = MultiIndex.from_arrays([['bar', 'bar', 'foo', 'foo', 'All'],
                                    ['one', 'two', 'one', 'two', '']],
                                   names=['a', 'b'])
        expected = DataFrame([[1, 0, 1], [1, 0, 1], [2, 0, 2], [1, 1, 2],
                              [5, 1, 6]], index=m)
        expected.columns = Index(['dull', 'shiny', 'All'], name='c')
        tm.assert_frame_equal(actual, expected)

    def test_crosstab_normalize(self):
        # Issue 12578
        df = pd.DataFrame({'a': [1, 2, 2, 2, 2], 'b': [3, 3, 4, 4, 4],
                           'c': [1, 1, np.nan, 1, 1]})

        rindex = pd.Index([1, 2], name='a')
        cindex = pd.Index([3, 4], name='b')
        full_normal = pd.DataFrame([[0.2, 0], [0.2, 0.6]],
                                   index=rindex, columns=cindex)
        row_normal = pd.DataFrame([[1.0, 0], [0.25, 0.75]],
                                  index=rindex, columns=cindex)
        col_normal = pd.DataFrame([[0.5, 0], [0.5, 1.0]],
                                  index=rindex, columns=cindex)

        # Check all normalize args
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize='all'),
                              full_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize=True),
                              full_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize='index'),
                              row_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize='columns'),
                              col_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize=1),
                              pd.crosstab(df.a, df.b, normalize='columns'))
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize=0),
                              pd.crosstab(df.a, df.b, normalize='index'))

        row_normal_margins = pd.DataFrame([[1.0, 0],
                                          [0.25, 0.75],
                                          [0.4, 0.6]],
                                          index=pd.Index([1, 2, 'All'],
                                                         name='a',
                                                         dtype='object'),
                                          columns=pd.Index([3, 4], name='b'))
        col_normal_margins = pd.DataFrame([[0.5, 0, 0.2], [0.5, 1.0, 0.8]],
                                          index=pd.Index([1, 2], name='a',
                                                         dtype='object'),
                                          columns=pd.Index([3, 4, 'All'],
                                                           name='b'))

        all_normal_margins = pd.DataFrame([[0.2, 0, 0.2],
                                          [0.2, 0.6, 0.8],
                                          [0.4, 0.6, 1]],
                                          index=pd.Index([1, 2, 'All'],
                                                         name='a',
                                                         dtype='object'),
                                          columns=pd.Index([3, 4, 'All'],
                                                           name='b'))

        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize='index',
                                          margins=True), row_normal_margins)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize='columns',
                                          margins=True), col_normal_margins)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize=True,
                                          margins=True), all_normal_margins)

        # Test arrays
        pd.crosstab([np.array([1, 1, 2, 2]), np.array([1, 2, 1, 2])],
                    np.array([1, 2, 1, 2]))

        # Test with aggfunc
        norm_counts = pd.DataFrame([[0.25, 0, 0.25],
                                    [0.25, 0.5, 0.75],
                                    [0.5, 0.5, 1]],
                                   index=pd.Index([1, 2, 'All'],
                                                  name='a',
                                                  dtype='object'),
                                   columns=pd.Index([3, 4, 'All'],
                                                    name='b'))
        test_case = pd.crosstab(df.a, df.b, df.c, aggfunc='count',
                                normalize='all',
                                margins=True)
        tm.assert_frame_equal(test_case, norm_counts)

        df = pd.DataFrame({'a': [1, 2, 2, 2, 2], 'b': [3, 3, 4, 4, 4],
                           'c': [0, 4, np.nan, 3, 3]})

        norm_sum = pd.DataFrame([[0, 0, 0.],
                                 [0.4, 0.6, 1],
                                 [0.4, 0.6, 1]],
                                index=pd.Index([1, 2, 'All'],
                                               name='a',
                                               dtype='object'),
                                columns=pd.Index([3, 4, 'All'],
                                                 name='b',
                                                 dtype='object'))
        test_case = pd.crosstab(df.a, df.b, df.c, aggfunc=np.sum,
                                normalize='all',
                                margins=True)
        tm.assert_frame_equal(test_case, norm_sum)

    def test_crosstab_with_empties(self):
        # Check handling of empties
        df = pd.DataFrame({'a': [1, 2, 2, 2, 2], 'b': [3, 3, 4, 4, 4],
                           'c': [np.nan, np.nan, np.nan, np.nan, np.nan]})

        empty = pd.DataFrame([[0.0, 0.0], [0.0, 0.0]],
                             index=pd.Index([1, 2],
                                            name='a',
                                            dtype='int64'),
                             columns=pd.Index([3, 4], name='b'))

        for i in [True, 'index', 'columns']:
            calculated = pd.crosstab(df.a, df.b, values=df.c, aggfunc='count',
                                     normalize=i)
            tm.assert_frame_equal(empty, calculated)

        nans = pd.DataFrame([[0.0, np.nan], [0.0, 0.0]],
                            index=pd.Index([1, 2],
                                           name='a',
                                           dtype='int64'),
                            columns=pd.Index([3, 4], name='b'))

        calculated = pd.crosstab(df.a, df.b, values=df.c, aggfunc='count',
                                 normalize=False)
        tm.assert_frame_equal(nans, calculated)

    def test_crosstab_errors(self):
        # Issue 12578

        df = pd.DataFrame({'a': [1, 2, 2, 2, 2], 'b': [3, 3, 4, 4, 4],
                           'c': [1, 1, np.nan, 1, 1]})

        error = 'values cannot be used without an aggfunc.'
        with tm.assertRaisesRegexp(ValueError, error):
            pd.crosstab(df.a, df.b, values=df.c)

        error = 'aggfunc cannot be used without values'
        with tm.assertRaisesRegexp(ValueError, error):
            pd.crosstab(df.a, df.b, aggfunc=np.mean)

        error = 'Not a valid normalize argument'
        with tm.assertRaisesRegexp(ValueError, error):
            pd.crosstab(df.a, df.b, normalize='42')

        with tm.assertRaisesRegexp(ValueError, error):
            pd.crosstab(df.a, df.b, normalize=42)

        error = 'Not a valid margins argument'
        with tm.assertRaisesRegexp(ValueError, error):
            pd.crosstab(df.a, df.b, normalize='all', margins=42)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
