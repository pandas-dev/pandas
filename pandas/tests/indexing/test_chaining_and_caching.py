import numpy as np
import pandas as pd
from pandas.core import common as com
from pandas import (compat, DataFrame, option_context,
                    Series, MultiIndex, date_range, Timestamp)
from pandas.util import testing as tm


class TestCaching(tm.TestCase):

    def test_slice_consolidate_invalidate_item_cache(self):

        # this is chained assignment, but will 'work'
        with option_context('chained_assignment', None):

            # #3970
            df = DataFrame({"aa": compat.lrange(5), "bb": [2.2] * 5})

            # Creates a second float block
            df["cc"] = 0.0

            # caches a reference to the 'bb' series
            df["bb"]

            # repr machinery triggers consolidation
            repr(df)

            # Assignment to wrong series
            df['bb'].iloc[0] = 0.17
            df._clear_item_cache()
            self.assertAlmostEqual(df['bb'][0], 0.17)

    def test_setitem_cache_updating(self):
        # GH 5424
        cont = ['one', 'two', 'three', 'four', 'five', 'six', 'seven']

        for do_ref in [False, False]:
            df = DataFrame({'a': cont,
                            "b": cont[3:] + cont[:3],
                            'c': np.arange(7)})

            # ref the cache
            if do_ref:
                df.ix[0, "c"]

            # set it
            df.ix[7, 'c'] = 1

            self.assertEqual(df.ix[0, 'c'], 0.0)
            self.assertEqual(df.ix[7, 'c'], 1.0)

        # GH 7084
        # not updating cache on series setting with slices
        expected = DataFrame({'A': [600, 600, 600]},
                             index=date_range('5/7/2014', '5/9/2014'))
        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        df = DataFrame({'C': ['A', 'A', 'A'], 'D': [100, 200, 300]})

        # loop through df to update out
        six = Timestamp('5/7/2014')
        eix = Timestamp('5/9/2014')
        for ix, row in df.iterrows():
            out.loc[six:eix, row['C']] = out.loc[six:eix, row['C']] + row['D']

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])

        # try via a chain indexing
        # this actually works
        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            v = out[row['C']][six:eix] + row['D']
            out[row['C']][six:eix] = v

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])

        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            out.loc[six:eix, row['C']] += row['D']

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])


class TestChaining(tm.TestCase):

    def test_setitem_chained_setfault(self):

        # GH6026
        # setfaults under numpy 1.7.1 (ok on 1.8)
        data = ['right', 'left', 'left', 'left', 'right', 'left', 'timeout']
        mdata = ['right', 'left', 'left', 'left', 'right', 'left', 'none']

        df = DataFrame({'response': np.array(data)})
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata}))

        recarray = np.rec.fromarrays([data], names=['response'])
        df = DataFrame(recarray)
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata}))

        df = DataFrame({'response': data, 'response1': data})
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata,
                                             'response1': data}))

        # GH 6056
        expected = DataFrame(dict(A=[np.nan, 'bar', 'bah', 'foo', 'bar']))
        df = DataFrame(dict(A=np.array(['foo', 'bar', 'bah', 'foo', 'bar'])))
        df['A'].iloc[0] = np.nan
        result = df.head()
        tm.assert_frame_equal(result, expected)

        df = DataFrame(dict(A=np.array(['foo', 'bar', 'bah', 'foo', 'bar'])))
        df.A.iloc[0] = np.nan
        result = df.head()
        tm.assert_frame_equal(result, expected)

    def test_detect_chained_assignment(self):

        pd.set_option('chained_assignment', 'raise')

        # work with the chain
        expected = DataFrame([[-5, 1], [-6, 3]], columns=list('AB'))
        df = DataFrame(np.arange(4).reshape(2, 2),
                       columns=list('AB'), dtype='int64')
        self.assertIsNone(df.is_copy)
        df['A'][0] = -5
        df['A'][1] = -6
        tm.assert_frame_equal(df, expected)

        # test with the chaining
        df = DataFrame({'A': Series(range(2), dtype='int64'),
                        'B': np.array(np.arange(2, 4), dtype=np.float64)})
        self.assertIsNone(df.is_copy)

        def f():
            df['A'][0] = -5

        self.assertRaises(com.SettingWithCopyError, f)

        def f():
            df['A'][1] = np.nan

        self.assertRaises(com.SettingWithCopyError, f)
        self.assertIsNone(df['A'].is_copy)

        # using a copy (the chain), fails
        df = DataFrame({'A': Series(range(2), dtype='int64'),
                        'B': np.array(np.arange(2, 4), dtype=np.float64)})

        def f():
            df.loc[0]['A'] = -5

        self.assertRaises(com.SettingWithCopyError, f)

        # doc example
        df = DataFrame({'a': ['one', 'one', 'two', 'three',
                              'two', 'one', 'six'],
                        'c': Series(range(7), dtype='int64')})
        self.assertIsNone(df.is_copy)
        expected = DataFrame({'a': ['one', 'one', 'two', 'three',
                                    'two', 'one', 'six'],
                              'c': [42, 42, 2, 3, 4, 42, 6]})

        def f():
            indexer = df.a.str.startswith('o')
            df[indexer]['c'] = 42

        self.assertRaises(com.SettingWithCopyError, f)

        expected = DataFrame({'A': [111, 'bbb', 'ccc'], 'B': [1, 2, 3]})
        df = DataFrame({'A': ['aaa', 'bbb', 'ccc'], 'B': [1, 2, 3]})

        def f():
            df['A'][0] = 111

        self.assertRaises(com.SettingWithCopyError, f)

        def f():
            df.loc[0]['A'] = 111

        self.assertRaises(com.SettingWithCopyError, f)

        df.loc[0, 'A'] = 111
        tm.assert_frame_equal(df, expected)

        # make sure that is_copy is picked up reconstruction
        # GH5475
        df = DataFrame({"A": [1, 2]})
        self.assertIsNone(df.is_copy)
        with tm.ensure_clean('__tmp__pickle') as path:
            df.to_pickle(path)
            df2 = pd.read_pickle(path)
            df2["B"] = df2["A"]
            df2["B"] = df2["A"]

        # a suprious raise as we are setting the entire column here
        # GH5597
        from string import ascii_letters as letters

        def random_text(nobs=100):
            df = []
            for i in range(nobs):
                idx = np.random.randint(len(letters), size=2)
                idx.sort()
                df.append([letters[idx[0]:idx[1]]])

            return DataFrame(df, columns=['letters'])

        df = random_text(100000)

        # always a copy
        x = df.iloc[[0, 1, 2]]
        self.assertIsNotNone(x.is_copy)
        x = df.iloc[[0, 1, 2, 4]]
        self.assertIsNotNone(x.is_copy)

        # explicity copy
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df = df.ix[indexer].copy()
        self.assertIsNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)

        # implicity take
        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df = df.ix[indexer]
        self.assertIsNotNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)

        # implicity take 2
        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df = df.ix[indexer]
        self.assertIsNotNone(df.is_copy)
        df.loc[:, 'letters'] = df['letters'].apply(str.lower)

        # should be ok even though it's a copy!
        self.assertIsNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)
        self.assertIsNone(df.is_copy)

        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df.ix[indexer, 'letters'] = df.ix[indexer, 'letters'].apply(str.lower)

        # an identical take, so no copy
        df = DataFrame({'a': [1]}).dropna()
        self.assertIsNone(df.is_copy)
        df['a'] += 1

        # inplace ops
        # original from:
        # http://stackoverflow.com/questions/20508968/series-fillna-in-a-multiindex-dataframe-does-not-fill-is-this-a-bug
        a = [12, 23]
        b = [123, None]
        c = [1234, 2345]
        d = [12345, 23456]
        tuples = [('eyes', 'left'), ('eyes', 'right'), ('ears', 'left'),
                  ('ears', 'right')]
        events = {('eyes', 'left'): a,
                  ('eyes', 'right'): b,
                  ('ears', 'left'): c,
                  ('ears', 'right'): d}
        multiind = MultiIndex.from_tuples(tuples, names=['part', 'side'])
        zed = DataFrame(events, index=['a', 'b'], columns=multiind)

        def f():
            zed['eyes']['right'].fillna(value=555, inplace=True)

        self.assertRaises(com.SettingWithCopyError, f)

        df = DataFrame(np.random.randn(10, 4))
        s = df.iloc[:, 0].sort_values()
        tm.assert_series_equal(s, df.iloc[:, 0].sort_values())
        tm.assert_series_equal(s, df[0].sort_values())

        # false positives GH6025
        df = DataFrame({'column1': ['a', 'a', 'a'], 'column2': [4, 8, 9]})
        str(df)
        df['column1'] = df['column1'] + 'b'
        str(df)
        df = df[df['column2'] != 8]
        str(df)
        df['column1'] = df['column1'] + 'c'
        str(df)

        # from SO:
        # http://stackoverflow.com/questions/24054495/potential-bug-setting-value-for-undefined-column-using-iloc
        df = DataFrame(np.arange(0, 9), columns=['count'])
        df['group'] = 'b'

        def f():
            df.iloc[0:5]['group'] = 'a'

        self.assertRaises(com.SettingWithCopyError, f)

        # mixed type setting
        # same dtype & changing dtype
        df = DataFrame(dict(A=date_range('20130101', periods=5),
                            B=np.random.randn(5),
                            C=np.arange(5, dtype='int64'),
                            D=list('abcde')))

        def f():
            df.ix[2]['D'] = 'foo'

        self.assertRaises(com.SettingWithCopyError, f)

        def f():
            df.ix[2]['C'] = 'foo'

        self.assertRaises(com.SettingWithCopyError, f)

        def f():
            df['C'][2] = 'foo'

        self.assertRaises(com.SettingWithCopyError, f)

    def test_setting_with_copy_bug(self):

        # operating on a copy
        df = pd.DataFrame({'a': list(range(4)),
                           'b': list('ab..'),
                           'c': ['a', 'b', np.nan, 'd']})
        mask = pd.isnull(df.c)

        def f():
            df[['c']][mask] = df[['b']][mask]

        self.assertRaises(com.SettingWithCopyError, f)

        # invalid warning as we are returning a new object
        # GH 8730
        df1 = DataFrame({'x': Series(['a', 'b', 'c']),
                         'y': Series(['d', 'e', 'f'])})
        df2 = df1[['x']]

        # this should not raise
        df2['y'] = ['g', 'h', 'i']

    def test_detect_chained_assignment_warnings(self):

        # warnings
        with option_context('chained_assignment', 'warn'):
            df = DataFrame({'A': ['aaa', 'bbb', 'ccc'], 'B': [1, 2, 3]})
            with tm.assert_produces_warning(
                    expected_warning=com.SettingWithCopyWarning):
                df.loc[0]['A'] = 111
