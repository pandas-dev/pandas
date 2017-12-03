# -*- coding: utf-8 -*-

""" test function application """

import pytest

from string import ascii_lowercase
from pandas import (date_range, Timestamp,
                    Index, MultiIndex, DataFrame, Series)
from pandas.util.testing import assert_frame_equal, assert_series_equal
from pandas.compat import product as cart_product

import numpy as np

import pandas.util.testing as tm
import pandas as pd
from .common import MixIn


# describe
# --------------------------------

class TestDescribe(MixIn):

    def test_apply_describe_bug(self):
        grouped = self.mframe.groupby(level='first')
        grouped.describe()  # it works!

    def test_series_describe_multikey(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.describe()
        assert_series_equal(result['mean'], grouped.mean(), check_names=False)
        assert_series_equal(result['std'], grouped.std(), check_names=False)
        assert_series_equal(result['min'], grouped.min(), check_names=False)

    def test_series_describe_single(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(lambda x: x.describe())
        expected = grouped.describe().stack()
        assert_series_equal(result, expected)

    def test_series_index_name(self):
        grouped = self.df.loc[:, ['C']].groupby(self.df['A'])
        result = grouped.agg(lambda x: x.mean())
        assert result.index.name == 'A'

    def test_frame_describe_multikey(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.describe()
        desc_groups = []
        for col in self.tsframe:
            group = grouped[col].describe()
            # GH 17464 - Remove duplicate MultiIndex levels
            group_col = pd.MultiIndex(
                levels=[[col], group.columns],
                labels=[[0] * len(group.columns), range(len(group.columns))])
            group = pd.DataFrame(group.values,
                                 columns=group_col,
                                 index=group.index)
            desc_groups.append(group)
        expected = pd.concat(desc_groups, axis=1)
        tm.assert_frame_equal(result, expected)

        groupedT = self.tsframe.groupby({'A': 0, 'B': 0,
                                         'C': 1, 'D': 1}, axis=1)
        result = groupedT.describe()
        expected = self.tsframe.describe().T
        expected.index = pd.MultiIndex(
            levels=[[0, 1], expected.index],
            labels=[[0, 0, 1, 1], range(len(expected.index))])
        tm.assert_frame_equal(result, expected)

    def test_frame_describe_tupleindex(self):

        # GH 14848 - regression from 0.19.0 to 0.19.1
        df1 = DataFrame({'x': [1, 2, 3, 4, 5] * 3,
                         'y': [10, 20, 30, 40, 50] * 3,
                         'z': [100, 200, 300, 400, 500] * 3})
        df1['k'] = [(0, 0, 1), (0, 1, 0), (1, 0, 0)] * 5
        df2 = df1.rename(columns={'k': 'key'})
        pytest.raises(ValueError, lambda: df1.groupby('k').describe())
        pytest.raises(ValueError, lambda: df2.groupby('key').describe())

    def test_frame_describe_unstacked_format(self):
        # GH 4792
        prices = {pd.Timestamp('2011-01-06 10:59:05', tz=None): 24990,
                  pd.Timestamp('2011-01-06 12:43:33', tz=None): 25499,
                  pd.Timestamp('2011-01-06 12:54:09', tz=None): 25499}
        volumes = {pd.Timestamp('2011-01-06 10:59:05', tz=None): 1500000000,
                   pd.Timestamp('2011-01-06 12:43:33', tz=None): 5000000000,
                   pd.Timestamp('2011-01-06 12:54:09', tz=None): 100000000}
        df = pd.DataFrame({'PRICE': prices,
                           'VOLUME': volumes})
        result = df.groupby('PRICE').VOLUME.describe()
        data = [df[df.PRICE == 24990].VOLUME.describe().values.tolist(),
                df[df.PRICE == 25499].VOLUME.describe().values.tolist()]
        expected = pd.DataFrame(data,
                                index=pd.Index([24990, 25499], name='PRICE'),
                                columns=['count', 'mean', 'std', 'min',
                                         '25%', '50%', '75%', 'max'])
        tm.assert_frame_equal(result, expected)


# nunique
# --------------------------------

class TestNUnique(MixIn):

    def test_series_groupby_nunique(self):

        def check_nunique(df, keys, as_index=True):
            for sort, dropna in cart_product((False, True), repeat=2):
                gr = df.groupby(keys, as_index=as_index, sort=sort)
                left = gr['julie'].nunique(dropna=dropna)

                gr = df.groupby(keys, as_index=as_index, sort=sort)
                right = gr['julie'].apply(Series.nunique, dropna=dropna)
                if not as_index:
                    right = right.reset_index(drop=True)

                assert_series_equal(left, right, check_names=False)

        days = date_range('2015-08-23', periods=10)

        for n, m in cart_product(10 ** np.arange(2, 6), (10, 100, 1000)):
            frame = DataFrame({
                'jim': np.random.choice(
                    list(ascii_lowercase), n),
                'joe': np.random.choice(days, n),
                'julie': np.random.randint(0, m, n)
            })

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])

            frame.loc[1::17, 'jim'] = None
            frame.loc[3::37, 'joe'] = None
            frame.loc[7::19, 'julie'] = None
            frame.loc[8::19, 'julie'] = None
            frame.loc[9::19, 'julie'] = None

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])
            check_nunique(frame, ['jim'], as_index=False)
            check_nunique(frame, ['jim', 'joe'], as_index=False)

    def test_nunique(self):
        df = DataFrame({
            'A': list('abbacc'),
            'B': list('abxacc'),
            'C': list('abbacx'),
        })

        expected = DataFrame({'A': [1] * 3, 'B': [1, 2, 1], 'C': [1, 1, 2]})
        result = df.groupby('A', as_index=False).nunique()
        tm.assert_frame_equal(result, expected)

        # as_index
        expected.index = list('abc')
        expected.index.name = 'A'
        result = df.groupby('A').nunique()
        tm.assert_frame_equal(result, expected)

        # with na
        result = df.replace({'x': None}).groupby('A').nunique(dropna=False)
        tm.assert_frame_equal(result, expected)

        # dropna
        expected = DataFrame({'A': [1] * 3, 'B': [1] * 3, 'C': [1] * 3},
                             index=list('abc'))
        expected.index.name = 'A'
        result = df.replace({'x': None}).groupby('A').nunique()
        tm.assert_frame_equal(result, expected)

    def test_nunique_with_object(self):
        # GH 11077
        data = pd.DataFrame(
            [[100, 1, 'Alice'],
             [200, 2, 'Bob'],
             [300, 3, 'Charlie'],
             [-400, 4, 'Dan'],
             [500, 5, 'Edith']],
            columns=['amount', 'id', 'name']
        )

        result = data.groupby(['id', 'amount'])['name'].nunique()
        index = MultiIndex.from_arrays([data.id, data.amount])
        expected = pd.Series([1] * 5, name='name', index=index)
        tm.assert_series_equal(result, expected)

    def test_nunique_with_empty_series(self):
        # GH 12553
        data = pd.Series(name='name')
        result = data.groupby(level=0).nunique()
        expected = pd.Series(name='name', dtype='int64')
        tm.assert_series_equal(result, expected)

    def test_nunique_with_timegrouper(self):
        # GH 13453
        test = pd.DataFrame({
            'time': [Timestamp('2016-06-28 09:35:35'),
                     Timestamp('2016-06-28 16:09:30'),
                     Timestamp('2016-06-28 16:46:28')],
            'data': ['1', '2', '3']}).set_index('time')
        result = test.groupby(pd.Grouper(freq='h'))['data'].nunique()
        expected = test.groupby(
            pd.Grouper(freq='h')
        )['data'].apply(pd.Series.nunique)
        tm.assert_series_equal(result, expected)


# count
# --------------------------------

class TestCount(MixIn):

    def test_groupby_timedelta_cython_count(self):
        df = DataFrame({'g': list('ab' * 2),
                        'delt': np.arange(4).astype('timedelta64[ns]')})
        expected = Series([
            2, 2
        ], index=pd.Index(['a', 'b'], name='g'), name='delt')
        result = df.groupby('g').delt.count()
        tm.assert_series_equal(expected, result)

    def test_count(self):
        n = 1 << 15
        dr = date_range('2015-08-30', periods=n // 10, freq='T')

        df = DataFrame({
            '1st': np.random.choice(
                list(ascii_lowercase), n),
            '2nd': np.random.randint(0, 5, n),
            '3rd': np.random.randn(n).round(3),
            '4th': np.random.randint(-10, 10, n),
            '5th': np.random.choice(dr, n),
            '6th': np.random.randn(n).round(3),
            '7th': np.random.randn(n).round(3),
            '8th': np.random.choice(dr, n) - np.random.choice(dr, 1),
            '9th': np.random.choice(
                list(ascii_lowercase), n)
        })

        for col in df.columns.drop(['1st', '2nd', '4th']):
            df.loc[np.random.choice(n, n // 10), col] = np.nan

        df['9th'] = df['9th'].astype('category')

        for key in '1st', '2nd', ['1st', '2nd']:
            left = df.groupby(key).count()
            right = df.groupby(key).apply(DataFrame.count).drop(key, axis=1)
            assert_frame_equal(left, right)

        # GH5610
        # count counts non-nulls
        df = pd.DataFrame([[1, 2, 'foo'],
                           [1, np.nan, 'bar'],
                           [3, np.nan, np.nan]],
                          columns=['A', 'B', 'C'])

        count_as = df.groupby('A').count()
        count_not_as = df.groupby('A', as_index=False).count()

        expected = DataFrame([[1, 2], [0, 0]], columns=['B', 'C'],
                             index=[1, 3])
        expected.index.name = 'A'
        assert_frame_equal(count_not_as, expected.reset_index())
        assert_frame_equal(count_as, expected)

        count_B = df.groupby('A')['B'].count()
        assert_series_equal(count_B, expected['B'])

    def test_count_object(self):
        df = pd.DataFrame({'a': ['a'] * 3 + ['b'] * 3, 'c': [2] * 3 + [3] * 3})
        result = df.groupby('c').a.count()
        expected = pd.Series([
            3, 3
        ], index=pd.Index([2, 3], name='c'), name='a')
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame({'a': ['a', np.nan, np.nan] + ['b'] * 3,
                           'c': [2] * 3 + [3] * 3})
        result = df.groupby('c').a.count()
        expected = pd.Series([
            1, 3
        ], index=pd.Index([2, 3], name='c'), name='a')
        tm.assert_series_equal(result, expected)

    def test_count_cross_type(self):  # GH8169
        vals = np.hstack((np.random.randint(0, 5, (100, 2)), np.random.randint(
            0, 2, (100, 2))))

        df = pd.DataFrame(vals, columns=['a', 'b', 'c', 'd'])
        df[df == 2] = np.nan
        expected = df.groupby(['c', 'd']).count()

        for t in ['float32', 'object']:
            df['a'] = df['a'].astype(t)
            df['b'] = df['b'].astype(t)
            result = df.groupby(['c', 'd']).count()
            tm.assert_frame_equal(result, expected)

    def test_lower_int_prec_count(self):
        df = DataFrame({'a': np.array(
            [0, 1, 2, 100], np.int8),
            'b': np.array(
            [1, 2, 3, 6], np.uint32),
            'c': np.array(
            [4, 5, 6, 8], np.int16),
            'grp': list('ab' * 2)})
        result = df.groupby('grp').count()
        expected = DataFrame({'a': [2, 2],
                              'b': [2, 2],
                              'c': [2, 2]}, index=pd.Index(list('ab'),
                                                           name='grp'))
        tm.assert_frame_equal(result, expected)

    def test_count_uses_size_on_exception(self):
        class RaisingObjectException(Exception):
            pass

        class RaisingObject(object):

            def __init__(self, msg='I will raise inside Cython'):
                super(RaisingObject, self).__init__()
                self.msg = msg

            def __eq__(self, other):
                # gets called in Cython to check that raising calls the method
                raise RaisingObjectException(self.msg)

        df = DataFrame({'a': [RaisingObject() for _ in range(4)],
                        'grp': list('ab' * 2)})
        result = df.groupby('grp').count()
        expected = DataFrame({'a': [2, 2]}, index=pd.Index(
            list('ab'), name='grp'))
        tm.assert_frame_equal(result, expected)


# size
# --------------------------------

class TestSize(MixIn):

    def test_size(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.size()
        for key, group in grouped:
            assert result[key] == len(group)

        grouped = self.df.groupby('A')
        result = grouped.size()
        for key, group in grouped:
            assert result[key] == len(group)

        grouped = self.df.groupby('B')
        result = grouped.size()
        for key, group in grouped:
            assert result[key] == len(group)

        df = DataFrame(np.random.choice(20, (1000, 3)), columns=list('abc'))
        for sort, key in cart_product((False, True), ('a', 'b', ['a', 'b'])):
            left = df.groupby(key, sort=sort).size()
            right = df.groupby(key, sort=sort)['c'].apply(lambda a: a.shape[0])
            assert_series_equal(left, right, check_names=False)

        # GH11699
        df = DataFrame([], columns=['A', 'B'])
        out = Series([], dtype='int64', index=Index([], name='A'))
        assert_series_equal(df.groupby('A').size(), out)
