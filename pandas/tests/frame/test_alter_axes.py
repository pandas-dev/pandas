# -*- coding: utf-8 -*-

from __future__ import print_function

import inspect
import pytest

from datetime import datetime, timedelta

import numpy as np

from pandas.compat import lrange, PY2
from pandas import (DataFrame, Series, Index, MultiIndex, RangeIndex,
                    IntervalIndex, DatetimeIndex, Categorical, cut,
                    Timestamp, date_range, to_datetime)
from pandas.core.dtypes.common import (
    is_object_dtype,
    is_categorical_dtype,
    is_interval_dtype)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


@pytest.fixture
def frame_of_index_cols():
    df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                    'B': ['one', 'two', 'three', 'one', 'two'],
                    'C': ['a', 'b', 'c', 'd', 'e'],
                    'D': np.random.randn(5),
                    'E': np.random.randn(5)})
    return df


class TestDataFrameAlterAxes(TestData):

    def test_set_index_directly(self):
        df = self.mixed_frame.copy()
        idx = Index(np.arange(len(df))[::-1])

        df.index = idx
        tm.assert_index_equal(df.index, idx)
        with tm.assert_raises_regex(ValueError, 'Length mismatch'):
            df.index = idx[::2]

    def test_set_index(self):
        df = self.mixed_frame.copy()
        idx = Index(np.arange(len(df))[::-1])

        df = df.set_index(idx)
        tm.assert_index_equal(df.index, idx)
        with tm.assert_raises_regex(ValueError, 'Length mismatch'):
            df.set_index(idx[::2])

    def test_set_index_cast(self):
        # issue casting an index then set_index
        df = DataFrame({'A': [1.1, 2.2, 3.3], 'B': [5.0, 6.1, 7.2]},
                       index=[2010, 2011, 2012])
        df2 = df.set_index(df.index.astype(np.int32))
        tm.assert_frame_equal(df, df2)

    # A has duplicate values, C does not
    @pytest.mark.parametrize('keys', ['A', 'C', ['A', 'B']])
    @pytest.mark.parametrize('inplace', [True, False])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_drop_inplace(self, frame_of_index_cols,
                                    drop, inplace, keys):
        df = frame_of_index_cols

        if isinstance(keys, list):
            idx = MultiIndex.from_arrays([df[x] for x in keys], names=keys)
        else:
            idx = Index(df[keys], name=keys)
        expected = df.drop(keys, axis=1) if drop else df
        expected.index = idx

        if inplace:
            result = df.copy()
            result.set_index(keys, drop=drop, inplace=True)
        else:
            result = df.set_index(keys, drop=drop)

        tm.assert_frame_equal(result, expected)

    # A has duplicate values, C does not
    @pytest.mark.parametrize('keys', ['A', 'C', ['A', 'B']])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_append(self, frame_of_index_cols, drop, keys):
        df = frame_of_index_cols

        keys = keys if isinstance(keys, list) else [keys]
        idx = MultiIndex.from_arrays([df.index] + [df[x] for x in keys],
                                     names=[None] + keys)
        expected = df.drop(keys, axis=1) if drop else df.copy()
        expected.index = idx

        result = df.set_index(keys, drop=drop, append=True)

        tm.assert_frame_equal(result, expected)

    # A has duplicate values, C does not
    @pytest.mark.parametrize('keys', ['A', 'C', ['A', 'B']])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_append_to_multiindex(self, frame_of_index_cols,
                                            drop, keys):
        # append to existing multiindex
        df = frame_of_index_cols.set_index(['D'], drop=drop, append=True)

        keys = keys if isinstance(keys, list) else [keys]
        expected = frame_of_index_cols.set_index(['D'] + keys,
                                                 drop=drop, append=True)

        result = df.set_index(keys, drop=drop, append=True)

        tm.assert_frame_equal(result, expected)

    def test_set_index_after_mutation(self):
        # GH1590
        df = DataFrame({'val': [0, 1, 2], 'key': ['a', 'b', 'c']})
        expected = DataFrame({'val': [1, 2]},
                             Index(['b', 'c'], name='key'))

        df2 = df.loc[df.index.map(lambda indx: indx >= 1)]
        result = df2.set_index('key')
        tm.assert_frame_equal(result, expected)

    # MultiIndex constructor does not work directly on Series -> lambda
    # also test index name if append=True (name is duplicate here for B)
    @pytest.mark.parametrize('box', [Series, Index, np.array,
                                     lambda x: MultiIndex.from_arrays([x])])
    @pytest.mark.parametrize('append, index_name', [(True, None),
                             (True, 'B'), (True, 'test'), (False, None)])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_pass_single_array(self, frame_of_index_cols,
                                         drop, append, index_name, box):
        df = frame_of_index_cols
        df.index.name = index_name

        key = box(df['B'])
        # np.array and list "forget" the name of B
        name = [None if box in [np.array, list] else 'B']

        result = df.set_index(key, drop=drop, append=append)

        # only valid column keys are dropped
        # since B is always passed as array above, nothing is dropped
        expected = df.set_index(['B'], drop=False, append=append)
        expected.index.names = [index_name] + name if append else name

        tm.assert_frame_equal(result, expected)

    # MultiIndex constructor does not work directly on Series -> lambda
    # also test index name if append=True (name is duplicate here for A & B)
    @pytest.mark.parametrize('box', [Series, Index, np.array, list,
                                     lambda x: MultiIndex.from_arrays([x])])
    @pytest.mark.parametrize('append, index_name',
                             [(True, None), (True, 'A'), (True, 'B'),
                              (True, 'test'), (False, None)])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_pass_arrays(self, frame_of_index_cols,
                                   drop, append, index_name, box):
        df = frame_of_index_cols
        df.index.name = index_name

        keys = ['A', box(df['B'])]
        # np.array and list "forget" the name of B
        names = ['A', None if box in [np.array, list] else 'B']

        result = df.set_index(keys, drop=drop, append=append)

        # only valid column keys are dropped
        # since B is always passed as array above, only A is dropped, if at all
        expected = df.set_index(['A', 'B'], drop=False, append=append)
        expected = expected.drop('A', axis=1) if drop else expected
        expected.index.names = [index_name] + names if append else names

        tm.assert_frame_equal(result, expected)

    # MultiIndex constructor does not work directly on Series -> lambda
    # We also emulate a "constructor" for the label -> lambda
    # also test index name if append=True (name is duplicate here for A)
    @pytest.mark.parametrize('box2', [Series, Index, np.array, list,
                                      lambda x: MultiIndex.from_arrays([x]),
                                      lambda x: x.name])
    @pytest.mark.parametrize('box1', [Series, Index, np.array, list,
                                      lambda x: MultiIndex.from_arrays([x]),
                                      lambda x: x.name])
    @pytest.mark.parametrize('append, index_name', [(True, None),
                             (True, 'A'), (True, 'test'), (False, None)])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_pass_arrays_duplicate(self, frame_of_index_cols, drop,
                                             append, index_name, box1, box2):
        df = frame_of_index_cols
        df.index.name = index_name

        keys = [box1(df['A']), box2(df['A'])]

        # == gives ambiguous Boolean for Series
        if keys[0] is 'A' and keys[1] is 'A':
            with tm.assert_raises_regex(ValueError,
                                        'Passed duplicate column names.*'):
                df.set_index(keys, drop=drop, append=append)
        else:
            result = df.set_index(keys, drop=drop, append=append)

            # to test against already-tested behavior, we add sequentially,
            # hence second append always True; must wrap in list, otherwise
            # list-box will be illegal
            expected = df.set_index([keys[0]], drop=drop, append=append)
            expected = expected.set_index([keys[1]], drop=drop, append=True)

            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('append', [True, False])
    @pytest.mark.parametrize('drop', [True, False])
    def test_set_index_pass_multiindex(self, frame_of_index_cols,
                                       drop, append):
        df = frame_of_index_cols
        keys = MultiIndex.from_arrays([df['A'], df['B']], names=['A', 'B'])

        result = df.set_index(keys, drop=drop, append=append)

        # setting with a MultiIndex will never drop columns
        expected = df.set_index(['A', 'B'], drop=False, append=append)

        tm.assert_frame_equal(result, expected)

    def test_set_index_verify_integrity(self, frame_of_index_cols):
        df = frame_of_index_cols

        with tm.assert_raises_regex(ValueError,
                                    'Index has duplicate keys'):
            df.set_index('A', verify_integrity=True)
        # with MultiIndex
        with tm.assert_raises_regex(ValueError,
                                    'Index has duplicate keys'):
            df.set_index([df['A'], df['A']], verify_integrity=True)

    def test_set_index_raise(self, frame_of_index_cols):
        df = frame_of_index_cols

        with tm.assert_raises_regex(KeyError, '.*'):  # column names are A-E
            df.set_index(['foo', 'bar', 'baz'], verify_integrity=True)

        # non-existent key in list with arrays
        with tm.assert_raises_regex(KeyError, '.*'):
            df.set_index([df['A'], df['B'], 'X'], verify_integrity=True)

    def test_construction_with_categorical_index(self):
        ci = tm.makeCategoricalIndex(10)
        ci.name = 'B'

        # with Categorical
        df = DataFrame({'A': np.random.randn(10),
                        'B': ci.values})
        idf = df.set_index('B')
        tm.assert_index_equal(idf.index, ci)

        # from a CategoricalIndex
        df = DataFrame({'A': np.random.randn(10),
                        'B': ci})
        idf = df.set_index('B')
        tm.assert_index_equal(idf.index, ci)

        # round-trip
        idf = idf.reset_index().set_index('B')
        tm.assert_index_equal(idf.index, ci)

    def test_set_index_cast_datetimeindex(self):
        df = DataFrame({'A': [datetime(2000, 1, 1) + timedelta(i)
                              for i in range(1000)],
                        'B': np.random.randn(1000)})

        idf = df.set_index('A')
        assert isinstance(idf.index, DatetimeIndex)

    def test_convert_dti_to_series(self):
        # don't cast a DatetimeIndex WITH a tz, leave as object
        # GH 6032
        idx = DatetimeIndex(to_datetime(['2013-1-1 13:00',
                                         '2013-1-2 14:00']),
                            name='B').tz_localize('US/Pacific')
        df = DataFrame(np.random.randn(2, 1), columns=['A'])

        expected = Series(np.array([Timestamp('2013-01-01 13:00:00-0800',
                                              tz='US/Pacific'),
                                    Timestamp('2013-01-02 14:00:00-0800',
                                              tz='US/Pacific')],
                                   dtype="object"), name='B')

        # convert index to series
        result = Series(idx)
        tm.assert_series_equal(result, expected)

        # assign to frame
        df['B'] = idx
        result = df['B']
        tm.assert_series_equal(result, expected)

        # convert to series while keeping the timezone
        result = idx.to_series(keep_tz=True, index=[0, 1])
        tm.assert_series_equal(result, expected)

        # convert to utc
        df['B'] = idx.to_series(index=[0, 1])
        result = df['B']
        comp = Series(DatetimeIndex(expected.values).tz_localize(None),
                      name='B')
        tm.assert_series_equal(result, comp)

        # list of datetimes with a tz
        df['B'] = idx.to_pydatetime()
        result = df['B']
        tm.assert_series_equal(result, expected)

        # GH 6785
        # set the index manually
        import pytz
        df = DataFrame(
            [{'ts': datetime(2014, 4, 1, tzinfo=pytz.utc), 'foo': 1}])
        expected = df.set_index('ts')
        df.index = df['ts']
        df.pop('ts')
        tm.assert_frame_equal(df, expected)

    def test_reset_index_tz(self, tz_aware_fixture):
        # GH 3950
        # reset_index with single level
        tz = tz_aware_fixture
        idx = date_range('1/1/2011', periods=5,
                         freq='D', tz=tz, name='idx')
        df = DataFrame({'a': range(5), 'b': ['A', 'B', 'C', 'D', 'E']},
                       index=idx)

        expected = DataFrame({'idx': [datetime(2011, 1, 1),
                                      datetime(2011, 1, 2),
                                      datetime(2011, 1, 3),
                                      datetime(2011, 1, 4),
                                      datetime(2011, 1, 5)],
                              'a': range(5),
                              'b': ['A', 'B', 'C', 'D', 'E']},
                             columns=['idx', 'a', 'b'])
        expected['idx'] = expected['idx'].apply(lambda d: Timestamp(d, tz=tz))
        tm.assert_frame_equal(df.reset_index(), expected)

    def test_set_index_timezone(self):
        # GH 12358
        # tz-aware Series should retain the tz
        idx = to_datetime(["2014-01-01 10:10:10"],
                          utc=True).tz_convert('Europe/Rome')
        df = DataFrame({'A': idx})
        assert df.set_index(idx).index[0].hour == 11
        assert DatetimeIndex(Series(df.A))[0].hour == 11
        assert df.set_index(df.A).index[0].hour == 11

    def test_set_index_dst(self):
        di = date_range('2006-10-29 00:00:00', periods=3,
                        freq='H', tz='US/Pacific')

        df = DataFrame(data={'a': [0, 1, 2], 'b': [3, 4, 5]},
                       index=di).reset_index()
        # single level
        res = df.set_index('index')
        exp = DataFrame(data={'a': [0, 1, 2], 'b': [3, 4, 5]},
                        index=Index(di, name='index'))
        tm.assert_frame_equal(res, exp)

        # GH 12920
        res = df.set_index(['index', 'a'])
        exp_index = MultiIndex.from_arrays([di, [0, 1, 2]],
                                           names=['index', 'a'])
        exp = DataFrame({'b': [3, 4, 5]}, index=exp_index)
        tm.assert_frame_equal(res, exp)

    def test_reset_index_with_intervals(self):
        idx = IntervalIndex.from_breaks(np.arange(11), name='x')
        original = DataFrame({'x': idx, 'y': np.arange(10)})[['x', 'y']]

        result = original.set_index('x')
        expected = DataFrame({'y': np.arange(10)}, index=idx)
        tm.assert_frame_equal(result, expected)

        result2 = result.reset_index()
        tm.assert_frame_equal(result2, original)

    def test_set_index_multiindexcolumns(self):
        columns = MultiIndex.from_tuples([('foo', 1), ('foo', 2), ('bar', 1)])
        df = DataFrame(np.random.randn(3, 3), columns=columns)
        result = df.set_index(df.columns[0])
        expected = df.iloc[:, 1:]
        expected.index = df.iloc[:, 0].values
        expected.index.names = [df.columns[0]]
        tm.assert_frame_equal(result, expected)

    def test_set_index_empty_column(self):
        # GH 1971
        df = DataFrame([
            {'a': 1, 'p': 0},
            {'a': 2, 'm': 10},
            {'a': 3, 'm': 11, 'p': 20},
            {'a': 4, 'm': 12, 'p': 21}
        ], columns=('a', 'm', 'p', 'x'))

        result = df.set_index(['a', 'x'])
        expected = df[['m', 'p']]
        expected.index = MultiIndex.from_arrays([df['a'], df['x']],
                                                names=['a', 'x'])
        tm.assert_frame_equal(result, expected)

    def test_set_columns(self):
        cols = Index(np.arange(len(self.mixed_frame.columns)))
        self.mixed_frame.columns = cols
        with tm.assert_raises_regex(ValueError, 'Length mismatch'):
            self.mixed_frame.columns = cols[::2]

    def test_dti_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = date_range('2011/01/01', periods=6, freq='M', tz='US/Eastern')
        idx2 = date_range('2013', periods=6, freq='A', tz='Asia/Tokyo')

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.reindex(idx2)
        tm.assert_index_equal(df.index, idx2)

        # GH 11314
        # with tz
        index = date_range(datetime(2015, 10, 1),
                           datetime(2015, 10, 1, 23),
                           freq='H', tz='US/Eastern')
        df = DataFrame(np.random.randn(24, 1), columns=['a'], index=index)
        new_index = date_range(datetime(2015, 10, 2),
                               datetime(2015, 10, 2, 23),
                               freq='H', tz='US/Eastern')

        result = df.set_index(new_index)
        assert result.index.freq == index.freq

    # Renaming

    def test_rename(self):
        mapping = {
            'A': 'a',
            'B': 'b',
            'C': 'c',
            'D': 'd'
        }

        renamed = self.frame.rename(columns=mapping)
        renamed2 = self.frame.rename(columns=str.lower)

        tm.assert_frame_equal(renamed, renamed2)
        tm.assert_frame_equal(renamed2.rename(columns=str.upper),
                              self.frame, check_names=False)

        # index
        data = {
            'A': {'foo': 0, 'bar': 1}
        }

        # gets sorted alphabetical
        df = DataFrame(data)
        renamed = df.rename(index={'foo': 'bar', 'bar': 'foo'})
        tm.assert_index_equal(renamed.index, Index(['foo', 'bar']))

        renamed = df.rename(index=str.upper)
        tm.assert_index_equal(renamed.index, Index(['BAR', 'FOO']))

        # have to pass something
        pytest.raises(TypeError, self.frame.rename)

        # partial columns
        renamed = self.frame.rename(columns={'C': 'foo', 'D': 'bar'})
        tm.assert_index_equal(renamed.columns, Index(['A', 'B', 'foo', 'bar']))

        # other axis
        renamed = self.frame.T.rename(index={'C': 'foo', 'D': 'bar'})
        tm.assert_index_equal(renamed.index, Index(['A', 'B', 'foo', 'bar']))

        # index with name
        index = Index(['foo', 'bar'], name='name')
        renamer = DataFrame(data, index=index)
        renamed = renamer.rename(index={'foo': 'bar', 'bar': 'foo'})
        tm.assert_index_equal(renamed.index,
                              Index(['bar', 'foo'], name='name'))
        assert renamed.index.name == renamer.index.name

    def test_rename_axis_inplace(self):
        # GH 15704
        frame = self.frame.copy()
        expected = frame.rename_axis('foo')
        result = frame.copy()
        no_return = result.rename_axis('foo', inplace=True)

        assert no_return is None
        tm.assert_frame_equal(result, expected)

        expected = frame.rename_axis('bar', axis=1)
        result = frame.copy()
        no_return = result.rename_axis('bar', axis=1, inplace=True)

        assert no_return is None
        tm.assert_frame_equal(result, expected)

    def test_rename_axis_warns(self):
        # https://github.com/pandas-dev/pandas/issues/17833
        df = DataFrame({"A": [1, 2], "B": [1, 2]})
        with tm.assert_produces_warning(FutureWarning) as w:
            df.rename_axis(id, axis=0)
            assert 'rename' in str(w[0].message)

        with tm.assert_produces_warning(FutureWarning) as w:
            df.rename_axis({0: 10, 1: 20}, axis=0)
            assert 'rename' in str(w[0].message)

        with tm.assert_produces_warning(FutureWarning) as w:
            df.rename_axis(id, axis=1)
            assert 'rename' in str(w[0].message)

        with tm.assert_produces_warning(FutureWarning) as w:
            df['A'].rename_axis(id)
            assert 'rename' in str(w[0].message)

    def test_rename_multiindex(self):

        tuples_index = [('foo1', 'bar1'), ('foo2', 'bar2')]
        tuples_columns = [('fizz1', 'buzz1'), ('fizz2', 'buzz2')]
        index = MultiIndex.from_tuples(tuples_index, names=['foo', 'bar'])
        columns = MultiIndex.from_tuples(
            tuples_columns, names=['fizz', 'buzz'])
        df = DataFrame([(0, 0), (1, 1)], index=index, columns=columns)

        #
        # without specifying level -> across all levels

        renamed = df.rename(index={'foo1': 'foo3', 'bar2': 'bar3'},
                            columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'})
        new_index = MultiIndex.from_tuples([('foo3', 'bar1'),
                                            ('foo2', 'bar3')],
                                           names=['foo', 'bar'])
        new_columns = MultiIndex.from_tuples([('fizz3', 'buzz1'),
                                              ('fizz2', 'buzz3')],
                                             names=['fizz', 'buzz'])
        tm.assert_index_equal(renamed.index, new_index)
        tm.assert_index_equal(renamed.columns, new_columns)
        assert renamed.index.names == df.index.names
        assert renamed.columns.names == df.columns.names

        #
        # with specifying a level (GH13766)

        # dict
        new_columns = MultiIndex.from_tuples([('fizz3', 'buzz1'),
                                              ('fizz2', 'buzz2')],
                                             names=['fizz', 'buzz'])
        renamed = df.rename(columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'},
                            level=0)
        tm.assert_index_equal(renamed.columns, new_columns)
        renamed = df.rename(columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'},
                            level='fizz')
        tm.assert_index_equal(renamed.columns, new_columns)

        new_columns = MultiIndex.from_tuples([('fizz1', 'buzz1'),
                                              ('fizz2', 'buzz3')],
                                             names=['fizz', 'buzz'])
        renamed = df.rename(columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'},
                            level=1)
        tm.assert_index_equal(renamed.columns, new_columns)
        renamed = df.rename(columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'},
                            level='buzz')
        tm.assert_index_equal(renamed.columns, new_columns)

        # function
        func = str.upper
        new_columns = MultiIndex.from_tuples([('FIZZ1', 'buzz1'),
                                              ('FIZZ2', 'buzz2')],
                                             names=['fizz', 'buzz'])
        renamed = df.rename(columns=func, level=0)
        tm.assert_index_equal(renamed.columns, new_columns)
        renamed = df.rename(columns=func, level='fizz')
        tm.assert_index_equal(renamed.columns, new_columns)

        new_columns = MultiIndex.from_tuples([('fizz1', 'BUZZ1'),
                                              ('fizz2', 'BUZZ2')],
                                             names=['fizz', 'buzz'])
        renamed = df.rename(columns=func, level=1)
        tm.assert_index_equal(renamed.columns, new_columns)
        renamed = df.rename(columns=func, level='buzz')
        tm.assert_index_equal(renamed.columns, new_columns)

        # index
        new_index = MultiIndex.from_tuples([('foo3', 'bar1'),
                                            ('foo2', 'bar2')],
                                           names=['foo', 'bar'])
        renamed = df.rename(index={'foo1': 'foo3', 'bar2': 'bar3'},
                            level=0)
        tm.assert_index_equal(renamed.index, new_index)

    def test_rename_nocopy(self):
        renamed = self.frame.rename(columns={'C': 'foo'}, copy=False)
        renamed['foo'] = 1.
        assert (self.frame['C'] == 1.).all()

    def test_rename_inplace(self):
        self.frame.rename(columns={'C': 'foo'})
        assert 'C' in self.frame
        assert 'foo' not in self.frame

        c_id = id(self.frame['C'])
        frame = self.frame.copy()
        frame.rename(columns={'C': 'foo'}, inplace=True)

        assert 'C' not in frame
        assert 'foo' in frame
        assert id(frame['foo']) != c_id

    def test_rename_bug(self):
        # GH 5344
        # rename set ref_locs, and set_index was not resetting
        df = DataFrame({0: ['foo', 'bar'], 1: ['bah', 'bas'], 2: [1, 2]})
        df = df.rename(columns={0: 'a'})
        df = df.rename(columns={1: 'b'})
        df = df.set_index(['a', 'b'])
        df.columns = ['2001-01-01']
        expected = DataFrame([[1], [2]],
                             index=MultiIndex.from_tuples(
                                 [('foo', 'bah'), ('bar', 'bas')],
                                 names=['a', 'b']),
                             columns=['2001-01-01'])
        tm.assert_frame_equal(df, expected)

    def test_rename_bug2(self):
        # GH 19497
        # rename was changing Index to MultiIndex if Index contained tuples

        df = DataFrame(data=np.arange(3), index=[(0, 0), (1, 1), (2, 2)],
                       columns=["a"])
        df = df.rename({(1, 1): (5, 4)}, axis="index")
        expected = DataFrame(data=np.arange(3), index=[(0, 0), (5, 4), (2, 2)],
                             columns=["a"])
        tm.assert_frame_equal(df, expected)

    def test_reorder_levels(self):
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]],
                           names=['L0', 'L1', 'L2'])
        df = DataFrame({'A': np.arange(6), 'B': np.arange(6)}, index=index)

        # no change, position
        result = df.reorder_levels([0, 1, 2])
        tm.assert_frame_equal(df, result)

        # no change, labels
        result = df.reorder_levels(['L0', 'L1', 'L2'])
        tm.assert_frame_equal(df, result)

        # rotate, position
        result = df.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(levels=[['one', 'two', 'three'], [0, 1], ['bar']],
                           labels=[[0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L1', 'L2', 'L0'])
        expected = DataFrame({'A': np.arange(6), 'B': np.arange(6)},
                             index=e_idx)
        tm.assert_frame_equal(result, expected)

        result = df.reorder_levels([0, 0, 0])
        e_idx = MultiIndex(levels=[['bar'], ['bar'], ['bar']],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L0', 'L0', 'L0'])
        expected = DataFrame({'A': np.arange(6), 'B': np.arange(6)},
                             index=e_idx)
        tm.assert_frame_equal(result, expected)

        result = df.reorder_levels(['L0', 'L0', 'L0'])
        tm.assert_frame_equal(result, expected)

    def test_reset_index(self):
        stacked = self.frame.stack()[::2]
        stacked = DataFrame({'foo': stacked, 'bar': stacked})

        names = ['first', 'second']
        stacked.index.names = names
        deleveled = stacked.reset_index()
        for i, (lev, lab) in enumerate(zip(stacked.index.levels,
                                           stacked.index.labels)):
            values = lev.take(lab)
            name = names[i]
            tm.assert_index_equal(values, Index(deleveled[name]))

        stacked.index.names = [None, None]
        deleveled2 = stacked.reset_index()
        tm.assert_series_equal(deleveled['first'], deleveled2['level_0'],
                               check_names=False)
        tm.assert_series_equal(deleveled['second'], deleveled2['level_1'],
                               check_names=False)

        # default name assigned
        rdf = self.frame.reset_index()
        exp = Series(self.frame.index.values, name='index')
        tm.assert_series_equal(rdf['index'], exp)

        # default name assigned, corner case
        df = self.frame.copy()
        df['index'] = 'foo'
        rdf = df.reset_index()
        exp = Series(self.frame.index.values, name='level_0')
        tm.assert_series_equal(rdf['level_0'], exp)

        # but this is ok
        self.frame.index.name = 'index'
        deleveled = self.frame.reset_index()
        tm.assert_series_equal(deleveled['index'], Series(self.frame.index))
        tm.assert_index_equal(deleveled.index,
                              Index(np.arange(len(deleveled))))

        # preserve column names
        self.frame.columns.name = 'columns'
        resetted = self.frame.reset_index()
        assert resetted.columns.name == 'columns'

        # only remove certain columns
        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index(['A', 'B'])

        # TODO should reset_index check_names ?
        tm.assert_frame_equal(rs, self.frame, check_names=False)

        rs = frame.reset_index(['index', 'A', 'B'])
        tm.assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index(['index', 'A', 'B'])
        tm.assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index('A')
        xp = self.frame.reset_index().set_index(['index', 'B'])
        tm.assert_frame_equal(rs, xp, check_names=False)

        # test resetting in place
        df = self.frame.copy()
        resetted = self.frame.reset_index()
        df.reset_index(inplace=True)
        tm.assert_frame_equal(df, resetted, check_names=False)

        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index('A', drop=True)
        xp = self.frame.copy()
        del xp['A']
        xp = xp.set_index(['B'], append=True)
        tm.assert_frame_equal(rs, xp, check_names=False)

    def test_reset_index_level(self):
        df = DataFrame([[1, 2, 3, 4], [5, 6, 7, 8]],
                       columns=['A', 'B', 'C', 'D'])

        for levels in ['A', 'B'], [0, 1]:
            # With MultiIndex
            result = df.set_index(['A', 'B']).reset_index(level=levels[0])
            tm.assert_frame_equal(result, df.set_index('B'))

            result = df.set_index(['A', 'B']).reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df.set_index('B'))

            result = df.set_index(['A', 'B']).reset_index(level=levels)
            tm.assert_frame_equal(result, df)

            result = df.set_index(['A', 'B']).reset_index(level=levels,
                                                          drop=True)
            tm.assert_frame_equal(result, df[['C', 'D']])

            # With single-level Index (GH 16263)
            result = df.set_index('A').reset_index(level=levels[0])
            tm.assert_frame_equal(result, df)

            result = df.set_index('A').reset_index(level=levels[:1])
            tm.assert_frame_equal(result, df)

            result = df.set_index(['A']).reset_index(level=levels[0],
                                                     drop=True)
            tm.assert_frame_equal(result, df[['B', 'C', 'D']])

        # Missing levels - for both MultiIndex and single-level Index:
        for idx_lev in ['A', 'B'], ['A']:
            with tm.assert_raises_regex(KeyError, 'Level E '):
                df.set_index(idx_lev).reset_index(level=['A', 'E'])
            with tm.assert_raises_regex(IndexError, 'Too many levels'):
                df.set_index(idx_lev).reset_index(level=[0, 1, 2])

    def test_reset_index_right_dtype(self):
        time = np.arange(0.0, 10, np.sqrt(2) / 2)
        s1 = Series((9.81 * time ** 2) / 2,
                    index=Index(time, name='time'),
                    name='speed')
        df = DataFrame(s1)

        resetted = s1.reset_index()
        assert resetted['time'].dtype == np.float64

        resetted = df.reset_index()
        assert resetted['time'].dtype == np.float64

    def test_reset_index_multiindex_col(self):
        vals = np.random.randn(3, 3).astype(object)
        idx = ['x', 'y', 'z']
        full = np.hstack(([[x] for x in idx], vals))
        df = DataFrame(vals, Index(idx, name='a'),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index()
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index(col_fill=None)
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index(col_level=1, col_fill='blah')
        xp = DataFrame(full, columns=[['blah', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

        df = DataFrame(vals,
                       MultiIndex.from_arrays([[0, 1, 2], ['x', 'y', 'z']],
                                              names=['d', 'a']),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index('a', )
        xp = DataFrame(full, Index([0, 1, 2], name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill=None)
        xp = DataFrame(full, Index(lrange(3), name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill='blah', col_level=1)
        xp = DataFrame(full, Index(lrange(3), name='d'),
                       columns=[['blah', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        tm.assert_frame_equal(rs, xp)

    def test_reset_index_multiindex_nan(self):
        # GH6322, testing reset_index on MultiIndexes
        # when we have a nan or all nan
        df = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [0, 1, np.nan],
                        'C': np.random.rand(3)})
        rs = df.set_index(['A', 'B']).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame({'A': [np.nan, 'b', 'c'],
                        'B': [0, 1, 2],
                        'C': np.random.rand(3)})
        rs = df.set_index(['A', 'B']).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [0, 1, 2],
                        'C': [np.nan, 1.1, 2.2]})
        rs = df.set_index(['A', 'B']).reset_index()
        tm.assert_frame_equal(rs, df)

        df = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [np.nan, np.nan, np.nan],
                        'C': np.random.rand(3)})
        rs = df.set_index(['A', 'B']).reset_index()
        tm.assert_frame_equal(rs, df)

    def test_reset_index_with_datetimeindex_cols(self):
        # GH5818
        #
        df = DataFrame([[1, 2], [3, 4]],
                       columns=date_range('1/1/2013', '1/2/2013'),
                       index=['A', 'B'])

        result = df.reset_index()
        expected = DataFrame([['A', 1, 2], ['B', 3, 4]],
                             columns=['index', datetime(2013, 1, 1),
                                      datetime(2013, 1, 2)])
        tm.assert_frame_equal(result, expected)

    def test_reset_index_range(self):
        # GH 12071
        df = DataFrame([[0, 0], [1, 1]], columns=['A', 'B'],
                       index=RangeIndex(stop=2))
        result = df.reset_index()
        assert isinstance(result.index, RangeIndex)
        expected = DataFrame([[0, 0, 0], [1, 1, 1]],
                             columns=['index', 'A', 'B'],
                             index=RangeIndex(stop=2))
        tm.assert_frame_equal(result, expected)

    def test_set_index_names(self):
        df = tm.makeDataFrame()
        df.index.name = 'name'

        assert df.set_index(df.index).index.names == ['name']

        mi = MultiIndex.from_arrays(df[['A', 'B']].T.values, names=['A', 'B'])
        mi2 = MultiIndex.from_arrays(df[['A', 'B', 'A', 'B']].T.values,
                                     names=['A', 'B', 'C', 'D'])

        df = df.set_index(['A', 'B'])

        assert df.set_index(df.index).index.names == ['A', 'B']

        # Check that set_index isn't converting a MultiIndex into an Index
        assert isinstance(df.set_index(df.index).index, MultiIndex)

        # Check actual equality
        tm.assert_index_equal(df.set_index(df.index).index, mi)

        idx2 = df.index.rename(['C', 'D'])

        # Check that [MultiIndex, MultiIndex] yields a MultiIndex rather
        # than a pair of tuples
        assert isinstance(df.set_index([df.index, idx2]).index, MultiIndex)

        # Check equality
        tm.assert_index_equal(df.set_index([df.index, idx2]).index, mi2)

    def test_rename_objects(self):
        renamed = self.mixed_frame.rename(columns=str.upper)

        assert 'FOO' in renamed
        assert 'foo' not in renamed

    def test_rename_axis_style(self):
        # https://github.com/pandas-dev/pandas/issues/12392
        df = DataFrame({"A": [1, 2], "B": [1, 2]}, index=['X', 'Y'])
        expected = DataFrame({"a": [1, 2], "b": [1, 2]}, index=['X', 'Y'])

        result = df.rename(str.lower, axis=1)
        tm.assert_frame_equal(result, expected)

        result = df.rename(str.lower, axis='columns')
        tm.assert_frame_equal(result, expected)

        result = df.rename({"A": 'a', 'B': 'b'}, axis=1)
        tm.assert_frame_equal(result, expected)

        result = df.rename({"A": 'a', 'B': 'b'}, axis='columns')
        tm.assert_frame_equal(result, expected)

        # Index
        expected = DataFrame({"A": [1, 2], "B": [1, 2]}, index=['x', 'y'])
        result = df.rename(str.lower, axis=0)
        tm.assert_frame_equal(result, expected)

        result = df.rename(str.lower, axis='index')
        tm.assert_frame_equal(result, expected)

        result = df.rename({'X': 'x', 'Y': 'y'}, axis=0)
        tm.assert_frame_equal(result, expected)

        result = df.rename({'X': 'x', 'Y': 'y'}, axis='index')
        tm.assert_frame_equal(result, expected)

        result = df.rename(mapper=str.lower, axis='index')
        tm.assert_frame_equal(result, expected)

    def test_rename_mapper_multi(self):
        df = DataFrame({"A": ['a', 'b'], "B": ['c', 'd'],
                        'C': [1, 2]}).set_index(["A", "B"])
        result = df.rename(str.upper)
        expected = df.rename(index=str.upper)
        tm.assert_frame_equal(result, expected)

    def test_rename_positional_named(self):
        # https://github.com/pandas-dev/pandas/issues/12392
        df = DataFrame({"a": [1, 2], "b": [1, 2]}, index=['X', 'Y'])
        result = df.rename(str.lower, columns=str.upper)
        expected = DataFrame({"A": [1, 2], "B": [1, 2]}, index=['x', 'y'])
        tm.assert_frame_equal(result, expected)

    def test_rename_axis_style_raises(self):
        # https://github.com/pandas-dev/pandas/issues/12392
        df = DataFrame({"A": [1, 2], "B": [1, 2]}, index=['0', '1'])

        # Named target and axis
        with tm.assert_raises_regex(TypeError, None):
            df.rename(index=str.lower, axis=1)

        with tm.assert_raises_regex(TypeError, None):
            df.rename(index=str.lower, axis='columns')

        with tm.assert_raises_regex(TypeError, None):
            df.rename(index=str.lower, axis='columns')

        with tm.assert_raises_regex(TypeError, None):
            df.rename(columns=str.lower, axis='columns')

        with tm.assert_raises_regex(TypeError, None):
            df.rename(index=str.lower, axis=0)

        # Multiple targets and axis
        with tm.assert_raises_regex(TypeError, None):
            df.rename(str.lower, str.lower, axis='columns')

        # Too many targets
        with tm.assert_raises_regex(TypeError, None):
            df.rename(str.lower, str.lower, str.lower)

        # Duplicates
        with tm.assert_raises_regex(TypeError, "multiple values"):
            df.rename(id, mapper=id)

    def test_reindex_api_equivalence(self):
        # equivalence of the labels/axis and index/columns API's
        df = DataFrame([[1, 2, 3], [3, 4, 5], [5, 6, 7]],
                       index=['a', 'b', 'c'],
                       columns=['d', 'e', 'f'])

        res1 = df.reindex(['b', 'a'])
        res2 = df.reindex(index=['b', 'a'])
        res3 = df.reindex(labels=['b', 'a'])
        res4 = df.reindex(labels=['b', 'a'], axis=0)
        res5 = df.reindex(['b', 'a'], axis=0)
        for res in [res2, res3, res4, res5]:
            tm.assert_frame_equal(res1, res)

        res1 = df.reindex(columns=['e', 'd'])
        res2 = df.reindex(['e', 'd'], axis=1)
        res3 = df.reindex(labels=['e', 'd'], axis=1)
        for res in [res2, res3]:
            tm.assert_frame_equal(res1, res)

        res1 = df.reindex(index=['b', 'a'], columns=['e', 'd'])
        res2 = df.reindex(columns=['e', 'd'], index=['b', 'a'])
        res3 = df.reindex(labels=['b', 'a'], axis=0).reindex(labels=['e', 'd'],
                                                             axis=1)
        for res in [res2, res3]:
            tm.assert_frame_equal(res1, res)

    def test_rename_positional(self):
        df = DataFrame(columns=['A', 'B'])
        with tm.assert_produces_warning(FutureWarning) as rec:
            result = df.rename(None, str.lower)
        expected = DataFrame(columns=['a', 'b'])
        tm.assert_frame_equal(result, expected)
        assert len(rec) == 1
        message = str(rec[0].message)
        assert 'rename' in message
        assert 'Use named arguments' in message

    def test_assign_columns(self):
        self.frame['hi'] = 'there'

        frame = self.frame.copy()
        frame.columns = ['foo', 'bar', 'baz', 'quux', 'foo2']
        tm.assert_series_equal(self.frame['C'], frame['baz'],
                               check_names=False)
        tm.assert_series_equal(self.frame['hi'], frame['foo2'],
                               check_names=False)

    def test_set_index_preserve_categorical_dtype(self):
        # GH13743, GH13854
        df = DataFrame({'A': [1, 2, 1, 1, 2],
                        'B': [10, 16, 22, 28, 34],
                        'C1': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=False),
                        'C2': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=True)})
        for cols in ['C1', 'C2', ['A', 'C1'], ['A', 'C2'], ['C1', 'C2']]:
            result = df.set_index(cols).reset_index()
            result = result.reindex(columns=df.columns)
            tm.assert_frame_equal(result, df)

    def test_ambiguous_warns(self):
        df = DataFrame({"A": [1, 2]})
        with tm.assert_produces_warning(FutureWarning):
            df.rename(id, id)

        with tm.assert_produces_warning(FutureWarning):
            df.rename({0: 10}, {"A": "B"})

    @pytest.mark.skipif(PY2, reason="inspect.signature")
    def test_rename_signature(self):
        sig = inspect.signature(DataFrame.rename)
        parameters = set(sig.parameters)
        assert parameters == {"self", "mapper", "index", "columns", "axis",
                              "inplace", "copy", "level"}

    @pytest.mark.skipif(PY2, reason="inspect.signature")
    def test_reindex_signature(self):
        sig = inspect.signature(DataFrame.reindex)
        parameters = set(sig.parameters)
        assert parameters == {"self", "labels", "index", "columns", "axis",
                              "limit", "copy", "level", "method",
                              "fill_value", "tolerance"}

    def test_droplevel(self):
        # GH20342
        df = DataFrame([
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12]
        ])
        df = df.set_index([0, 1]).rename_axis(['a', 'b'])
        df.columns = MultiIndex.from_tuples([('c', 'e'), ('d', 'f')],
                                            names=['level_1', 'level_2'])

        # test that dropping of a level in index works
        expected = df.reset_index('a', drop=True)
        result = df.droplevel('a', axis='index')
        tm.assert_frame_equal(result, expected)

        # test that dropping of a level in columns works
        expected = df.copy()
        expected.columns = Index(['c', 'd'], name='level_1')
        result = df.droplevel('level_2', axis='columns')
        tm.assert_frame_equal(result, expected)


class TestIntervalIndex(object):

    def test_setitem(self):

        df = DataFrame({'A': range(10)})
        s = cut(df.A, 5)
        assert isinstance(s.cat.categories, IntervalIndex)

        # B & D end up as Categoricals
        # the remainer are converted to in-line objects
        # contining an IntervalIndex.values
        df['B'] = s
        df['C'] = np.array(s)
        df['D'] = s.values
        df['E'] = np.array(s.values)

        assert is_categorical_dtype(df['B'])
        assert is_interval_dtype(df['B'].cat.categories)
        assert is_categorical_dtype(df['D'])
        assert is_interval_dtype(df['D'].cat.categories)

        assert is_object_dtype(df['C'])
        assert is_object_dtype(df['E'])

        # they compare equal as Index
        # when converted to numpy objects
        c = lambda x: Index(np.array(x))
        tm.assert_index_equal(c(df.B), c(df.B), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.C), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.D), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.D), check_names=False)

        # B & D are the same Series
        tm.assert_series_equal(df['B'], df['B'], check_names=False)
        tm.assert_series_equal(df['B'], df['D'], check_names=False)

        # C & E are the same Series
        tm.assert_series_equal(df['C'], df['C'], check_names=False)
        tm.assert_series_equal(df['C'], df['E'], check_names=False)

    def test_set_reset_index(self):

        df = DataFrame({'A': range(10)})
        s = cut(df.A, 5)
        df['B'] = s
        df = df.set_index('B')

        df = df.reset_index()

    def test_set_axis_inplace(self):
        # GH14636
        df = DataFrame({'A': [1.1, 2.2, 3.3],
                        'B': [5.0, 6.1, 7.2],
                        'C': [4.4, 5.5, 6.6]},
                       index=[2010, 2011, 2012])

        expected = {0: df.copy(),
                    1: df.copy()}
        expected[0].index = list('abc')
        expected[1].columns = list('abc')
        expected['index'] = expected[0]
        expected['columns'] = expected[1]

        for axis in expected:
            # inplace=True
            # The FutureWarning comes from the fact that we would like to have
            # inplace default to False some day
            for inplace, warn in (None, FutureWarning), (True, None):
                kwargs = {'inplace': inplace}

                result = df.copy()
                with tm.assert_produces_warning(warn):
                    result.set_axis(list('abc'), axis=axis, **kwargs)
                tm.assert_frame_equal(result, expected[axis])

            # inplace=False
            result = df.set_axis(list('abc'), axis=axis, inplace=False)
            tm.assert_frame_equal(expected[axis], result)

        # omitting the "axis" parameter
        with tm.assert_produces_warning(None):
            result = df.set_axis(list('abc'), inplace=False)
        tm.assert_frame_equal(result, expected[0])

        # wrong values for the "axis" parameter
        for axis in 3, 'foo':
            with tm.assert_raises_regex(ValueError, 'No axis named'):
                df.set_axis(list('abc'), axis=axis, inplace=False)

    def test_set_axis_prior_to_deprecation_signature(self):
        df = DataFrame({'A': [1.1, 2.2, 3.3],
                        'B': [5.0, 6.1, 7.2],
                        'C': [4.4, 5.5, 6.6]},
                       index=[2010, 2011, 2012])

        expected = {0: df.copy(),
                    1: df.copy()}
        expected[0].index = list('abc')
        expected[1].columns = list('abc')
        expected['index'] = expected[0]
        expected['columns'] = expected[1]

        # old signature
        for axis in expected:
            with tm.assert_produces_warning(FutureWarning):
                result = df.set_axis(axis, list('abc'), inplace=False)
            tm.assert_frame_equal(result, expected[axis])
