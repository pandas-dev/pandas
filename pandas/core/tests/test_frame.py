# pylint: disable-msg=W0612
from copy import deepcopy
from datetime import datetime, timedelta
from cStringIO import StringIO
import cPickle as pickle
import os
import unittest

from numpy import random
import numpy as np

import pandas.core.datetools as datetools
from pandas.core.api import DataFrame, Index, Series, notnull

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 randn)

import pandas.util.testing as common

#-------------------------------------------------------------------------------
# DataFrame test cases

class TestDataFrame(unittest.TestCase):
    klass = DataFrame

    def setUp(self):
        self.seriesd = common.getSeriesData()
        self.tsd = common.getTimeSeriesData()

        self.frame = self.klass(self.seriesd)
        self.tsframe = self.klass(self.tsd)

        self.mixed_frame = self.frame.copy()
        self.mixed_frame['foo'] = 'bar'

        self.ts1 = common.makeTimeSeries()
        self.ts2 = common.makeTimeSeries()[5:]
        self.ts3 = common.makeTimeSeries()[-5:]
        self.ts4 = common.makeTimeSeries()[1:-1]

        self.ts_dict = {
            'col1' : self.ts1,
            'col2' : self.ts2,
            'col3' : self.ts3,
            'col4' : self.ts4,
        }
        self.empty = self.klass({})

        self.unsortable = self.klass(
            {'foo' : [1] * 1000,
             datetime.today() : [1] * 1000,
             'bar' : ['bar'] * 1000,
             datetime.today() + timedelta(1) : ['bar'] * 1000},
            index=np.arange(1000))

    def test_set_index(self):
        idx = Index(np.arange(len(self.mixed_frame)))
        self.mixed_frame.index = idx

        self.assert_(self.mixed_frame['foo'].index  is idx)

    def test_constructor(self):
        df = self.klass()
        self.assert_(len(df.index) == 0)

        df = self.klass(data={})
        self.assert_(len(df.index) == 0)

    def test_constructor_mixed(self):
        index, data = common.getMixedTypeDict()

        indexed_frame = self.klass(data, index=index)
        unindexed_frame = self.klass(data)

    def test_constructor_dict(self):
        frame = self.klass({'col1' : self.ts1,
                            'col2' : self.ts2})

        common.assert_dict_equal(self.ts1, frame['col1'], compare_keys=False)
        common.assert_dict_equal(self.ts2, frame['col2'], compare_keys=False)

        frame = self.klass({'col1' : self.ts1,
                            'col2' : self.ts2},
                           columns=['col2', 'col3', 'col4'])

        self.assertEqual(len(frame), len(self.ts2))
        self.assert_('col1' not in frame)
        self.assert_(np.isnan(frame['col3']).all())

        # Corner cases
        self.assertEqual(len(self.klass({})), 0)
        self.assertRaises(Exception, lambda x: self.klass([self.ts1, self.ts2]))

        # pass dict and array, nicht nicht
        self.assertRaises(Exception, self.klass,
                          {'A' : {'a' : 'a', 'b' : 'b'},
                           'B' : ['a', 'b']})

        # can I rely on the order?
        self.assertRaises(Exception, self.klass,
                          {'A' : ['a', 'b'],
                           'B' : {'a' : 'a', 'b' : 'b'}})
        self.assertRaises(Exception, self.klass,
                          {'A' : ['a', 'b'],
                           'B' : Series(['a', 'b'], index=['a', 'b'])})

        # Length-one dict micro-optimization
        frame = self.klass({'A' : {'1' : 1, '2' : 2}})
        self.assert_(np.array_equal(frame.index, ['1', '2']))

        # empty dict plus index
        idx = Index([0, 1, 2])
        frame = self.klass({}, index=idx)
        self.assert_(frame.index is idx)

    def test_constructor_dict_cast(self):
        # cast float tests
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        frame = self.klass(test_data, dtype=float)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.float_)
        self.assert_(frame['A'].dtype == np.float_)

        frame = self.klass(test_data)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.object_)
        self.assert_(frame['A'].dtype == np.float_)

        # can't cast to float
        test_data = {
                'A' : dict(zip(range(20), common.makeDateIndex(20))),
                'B' : dict(zip(range(15), randn(15)))
        }
        frame = self.klass(test_data, dtype=float)
        self.assertEqual(len(frame), 20)
        self.assert_(frame['A'].dtype == np.object_)
        self.assert_(frame['B'].dtype == np.float_)

    def test_constructor_ndarray(self):
        mat = np.zeros((2, 3), dtype=float)

        # 2-D input
        frame = self.klass(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.cols()), 3)

        # cast type
        frame = self.klass(mat, columns=['A', 'B', 'C'],
                           index=[1, 2], dtype=int)
        self.assert_(frame.values.dtype == np.int_)

        # 1-D input
        frame = self.klass(np.zeros(3), columns=['A'], index=[1, 2, 3])
        self.assertEqual(len(frame.index), 3)
        self.assertEqual(len(frame.cols()), 1)

        frame = self.klass(['foo', 'bar'], index=[0, 1], columns=['A'])
        self.assertEqual(len(frame), 2)

        # higher dim raise exception
        self.assertRaises(Exception, self.klass, np.zeros((3, 3, 3)),
                          columns=['A', 'B', 'C'], index=[1])

        # wrong size axis labels
        self.assertRaises(Exception, self.klass, mat,
                          columns=['A', 'B', 'C'], index=[1])

        self.assertRaises(Exception, self.klass, mat,
                          columns=['A', 'B'], index=[1, 2])

        # have to pass columns and index
        self.assertRaises(Exception, self.klass, mat, index=[1])
        self.assertRaises(Exception, self.klass, mat, columns=['A', 'B', 'C'])

    def test_array_interface(self):
        result = np.sqrt(self.frame)
        self.assert_(type(result) is type(self.frame))
        self.assert_(result.index is self.frame.index)
        self.assert_(result.cols() == self.frame.cols())

        assert_frame_equal(result, self.frame.apply(np.sqrt))

    def test_pickle(self):
        unpickled = pickle.loads(pickle.dumps(self.mixed_frame))
        assert_frame_equal(self.mixed_frame, unpickled)

    def test_toDict(self):
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        recons_data = self.klass(test_data).toDict()

        for k, v in test_data.iteritems():
            for k2, v2 in v.iteritems():
                self.assertEqual(v2, recons_data[k][k2])

    def test_fromRecords(self):
        # from numpy documentation
        arr = np.zeros((2,),dtype=('i4,f4,a10'))
        arr[:] = [(1,2.,'Hello'),(2,3.,"World")]

        frame = self.klass.fromRecords(arr)
        indexed_frame = self.klass.fromRecords(arr, indexField='f1')

        self.assertRaises(Exception, self.klass.fromRecords, np.zeros((2, 3)))

        # what to do?
        records = indexed_frame.toRecords()
        self.assertEqual(len(records.dtype.names), 3)

    def test_get_agg_axis(self):
        cols = self.frame._get_agg_axis(0)
        self.assert_(list(cols) == list(self.frame.columns))

        idx = self.frame._get_agg_axis(1)
        self.assert_(idx is self.frame.index)

        self.assertRaises(Exception, self.frame._get_agg_axis, 2)

    def test_nonzero(self):
        self.assertFalse(self.empty)

        self.assert_(self.frame)
        self.assert_(self.mixed_frame)

        # corner case
        df = self.klass({'A' : [1., 2., 3.],
                         'B' : ['a', 'b', 'c']},
                        index=np.arange(3))
        del df['A']
        self.assert_(df)

    def test_repr(self):
        # empty
        foo = repr(self.empty)

        # empty with index
        frame = self.klass(index=np.arange(1000))
        foo = repr(frame)

        # small one
        foo = repr(self.frame)

        # big one
        biggie = self.klass(np.zeros((1000, 4)), columns=range(4),
                            index=range(1000))
        foo = repr(biggie)

        # mixed
        foo = repr(self.mixed_frame)

        # big mixed
        biggie = self.klass({'A' : randn(1000),
                             'B' : common.makeStringIndex(1000)},
                            index=range(1000))
        biggie['A'][:20] = np.NaN
        biggie['B'][:20] = np.NaN

        foo = repr(biggie)

        # exhausting cases in DataMatrix.info

        # columns but no index
        no_index = self.klass(columns=[0, 1, 3])
        foo = repr(no_index)

        # no columns or index
        buf = StringIO()
        self.empty.info(buffer=buf)

        # columns are not sortable
        foo = repr(self.unsortable)

    def test_toString(self):
        # big mixed
        biggie = self.klass({'A' : randn(1000),
                             'B' : common.makeStringIndex(1000)},
                            index=range(1000))

        biggie['A'][:20] = np.NaN
        biggie['B'][:20] = np.NaN
        buf = StringIO()
        biggie.toString(buffer=buf)

        biggie.toString(buffer=buf, columns=['B', 'A'], colSpace=17)
        biggie.toString(buffer=buf, columns=['B', 'A'],
                        formatters={'A' : lambda x: '%.1f' % x})

        biggie.toString(buffer=buf, columns=['B', 'A'],
                        float_format=str)

        frame = self.klass(index=np.arange(1000))
        frame.toString(buffer=buf)

    def test_getitem(self):
        # slicing

        sl = self.frame[:20]
        self.assertEqual(20, len(sl.index))

        # column access

        for _, series in sl.iteritems():
            self.assertEqual(20, len(series.index))
            self.assert_(common.equalContents(series.index, sl.index))

        for key, _ in self.frame._series.iteritems():
            self.assert_(self.frame[key] is not None)

        self.assert_('random' not in self.frame)
        self.assertRaises(Exception, self.frame.__getitem__, 'random')

        # boolean indexing
        d = self.tsframe.index[10]
        indexer = self.tsframe.index > d

        subindex = self.tsframe.index[indexer]
        subframe = self.tsframe[indexer]

        self.assert_(np.array_equal(subindex, subframe.index))
        self.assertRaises(Exception, self.tsframe.__getitem__, indexer[:-1])

    def test_setitem(self):
        # not sure what else to do here
        series = self.frame['A'][::2]
        self.frame['col5'] = series
        self.assert_('col5' in self.frame)
        common.assert_dict_equal(series, self.frame['col5'],
                                 compare_keys=False)

        series = self.frame['A']
        self.frame['col6'] = series
        common.assert_dict_equal(series, self.frame['col6'],
                                 compare_keys=False)

        self.assertRaises(Exception, self.frame.__setitem__,
                          randn(len(self.frame) + 1))

        # set ndarray
        arr = randn(len(self.frame))
        self.frame['col9'] = arr
        self.assert_((self.frame['col9'] == arr).all())

        # set value, do out of order for DataMatrix
        self.frame['col7'] = 5
        assert((self.frame['col7'] == 5).all())

        self.frame['col8'] = 'foo'
        assert((self.frame['col8'] == 'foo').all())

        smaller = self.frame[:2]
        smaller['col10'] = ['1', '2']
        self.assertEqual(smaller['col10'].dtype, np.object_)
        self.assert_((smaller['col10'] == ['1', '2']).all())

    def test_delitem(self):
        del self.frame['A']
        self.assert_('A' not in self.frame)

    def test_pop(self):
        A = self.frame.pop('A')
        self.assert_('A' not in self.frame)

        self.frame['foo'] = 'bar'
        foo = self.frame.pop('foo')
        self.assert_('foo' not in self.frame)

    def test_iter(self):
        self.assert_(common.equalContents(list(self.frame), self.frame.cols()))

    def test_len(self):
        self.assertEqual(len(self.frame), len(self.frame.index))

    def test_operators(self):
        garbage = random.random(4)
        colSeries = Series(garbage, index=np.array(self.frame.cols()))

        idSum = self.frame + self.frame
        seriesSum = self.frame + colSeries

        for col, series in idSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] * 2
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(np.isnan(origVal))

        for col, series in seriesSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] + colSeries[col]
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(np.isnan(origVal))

    def test_neg(self):
        # what to do?
        assert_frame_equal(-self.frame, -1 * self.frame)

    def test_firstLastValid(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = np.NaN
        mat[-5:] = np.NaN

        frame = self.klass({'foo' : mat}, index=self.frame.index)
        index = frame._firstTimeWithValue()

        self.assert_(index == frame.index[5])

        index = frame._lastTimeWithValue()
        self.assert_(index == frame.index[-6])

    def test_combineFrame(self):
        frame_copy = self.frame.reindex(self.frame.index[::2])

        del frame_copy['D']
        frame_copy['C'][:5] = np.NaN

        added = self.frame + frame_copy
        common.assert_dict_equal(added['A'].valid(),
                                 self.frame['A'] * 2,
                                 compare_keys=False)

        self.assert_(np.isnan(added['C'][:5]).all())
        self.assert_(np.isnan(added['D']).all())

        self_added = self.frame + self.frame
        self.assert_(self_added.index.equals(self.frame.index))

        added_rev = frame_copy + self.frame
        self.assert_(np.isnan(added['D']).all())

        # corner cases

        # empty
        plus_empty = self.frame + self.empty
        self.assert_(np.isnan(plus_empty.values).all())

        empty_plus = self.empty + self.frame
        self.assert_(np.isnan(empty_plus.values).all())

        empty_empty = self.empty + self.empty
        self.assert_(not empty_empty)

        # out of order
        reverse = self.frame.reindex(columns=self.frame.columns[::-1])

        assert_frame_equal(reverse + self.frame, self.frame * 2)

    def test_combineSeries(self):

        # Series
        series = self.frame.getXS(self.frame.index[0])

        added = self.frame + series

        for key, s in added.iteritems():
            assert_series_equal(s, self.frame[key] + series[key])

        larger_series = series.toDict()
        larger_series['E'] = 1
        larger_series = Series(larger_series)
        larger_added = self.frame + larger_series

        for key, s in self.frame.iteritems():
            assert_series_equal(larger_added[key], s + series[key])
        self.assert_('E' in larger_added)
        self.assert_(np.isnan(larger_added['E']).all())

        # TimeSeries
        ts = self.tsframe['A']
        added = self.tsframe + ts

        for key, col in self.tsframe.iteritems():
            assert_series_equal(added[key], col + ts)

        smaller_frame = self.tsframe[:-5]
        smaller_added = smaller_frame + ts

        self.assert_(smaller_added.index.equals(self.tsframe.index))

        # length 0
        result = self.tsframe + ts[:0]

        # Frame is length 0
        result = self.tsframe[:0] + ts
        self.assertEqual(len(result), 0)

        # empty but with non-empty index
        frame = self.tsframe[:1].reindex(columns=[])
        result = frame * ts
        self.assertEqual(len(result), len(ts))

    def test_combineFunc(self):
        result = self.frame * 2
        self.assert_(np.array_equal(result.values, self.frame.values * 2))

        result = self.empty * 2
        self.assert_(result.index is self.empty.index)
        self.assertEqual(len(result.columns), 0)

    def test_toCSV(self):
        path = '__tmp__'

        self.frame['A'][:5] = np.NaN

        self.frame.toCSV(path)
        self.frame.toCSV(path, cols=['A', 'B'])
        self.frame.toCSV(path, header=False)
        self.frame.toCSV(path, index=False)

        os.remove(path)

    def test_toDataMatrix(self):
        dm = self.frame.toDataMatrix()

    def test_info(self):
        pass

    def test_rows(self):
        self.assert_(self.tsframe.rows() is self.tsframe.index)

    def test_cols(self):
        self.assert_(self.tsframe.cols() == list(self.tsframe.columns))

    def test_columns(self):
        pass

    def test_iteritems(self):
        pass

    def test_append(self):
        begin_index = self.frame.index[:5]
        end_index = self.frame.index[5:]

        begin_frame = self.frame.reindex(begin_index)
        end_frame = self.frame.reindex(end_index)

        appended = begin_frame.append(end_frame)
        assert_almost_equal(appended['A'], self.frame['A'])

        del end_frame['A']
        partial_appended = begin_frame.append(end_frame)
        self.assert_('A' in partial_appended)

        partial_appended = end_frame.append(begin_frame)
        self.assert_('A' in partial_appended)

        # mixed type handling
        appended = self.mixed_frame[:5].append(self.mixed_frame[5:])
        assert_frame_equal(appended, self.mixed_frame)

        # what to test here
        mixed_appended = self.mixed_frame[:5].append(self.frame[5:])
        mixed_appended2 = self.frame[:5].append(self.mixed_frame[5:])

        assert_frame_equal(mixed_appended, mixed_appended2)

        # append empty
        appended = self.frame.append(self.empty)
        assert_frame_equal(self.frame, appended)
        self.assert_(appended is not self.frame)

        appended = self.empty.append(self.frame)
        assert_frame_equal(self.frame, appended)
        self.assert_(appended is not self.frame)

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(datetools.bmonthEnd)
        rule_monthly = self.tsframe.asfreq('EOM')

        assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('WEEKDAY', fillMethod='pad')

    def test_asMatrix(self):
        frame = self.frame
        mat = frame.asMatrix()
        smallerCols = ['C', 'A']

        frameCols = frame.cols()
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])

        # mixed type
        mat = self.mixed_frame.asMatrix(['foo', 'A'])
        self.assertEqual(mat[0, 0], 'bar')

    def test_values(self):
        pass

    def test_deepcopy(self):
        cp = deepcopy(self.frame)
        series = cp['A']
        series[:] = 10
        for idx, value in series.iteritems():
            self.assertNotEqual(self.frame['A'][idx], value)

    def test_copy(self):
        pass

    def test_corr(self):
        self.frame['A'][:5] = np.NaN
        self.frame['B'][:10] = np.NaN

        correls = self.frame.corr()

        assert_almost_equal(correls['A']['C'],
                            self.frame['A'].corr(self.frame['C']))

    def test_dropEmptyRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = np.NaN

        frame = self.klass({'foo' : mat}, index=self.frame.index)

        smaller_frame = frame.dropEmptyRows()
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        smaller_frame = frame.dropEmptyRows(['foo'])
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = np.NaN

        frame = self.klass({'foo' : mat}, index=self.frame.index)
        frame['bar'] = 5

        smaller_frame = frame.dropIncompleteRows()
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        samesize_frame = frame.dropIncompleteRows(specificColumns=['bar'])
        self.assert_(samesize_frame.index.equals(self.frame.index))

    def test_fill(self):
        self.tsframe['A'][:5] = np.NaN
        self.tsframe['A'][-5:] = np.NaN

        zero_filled = self.tsframe.fill(0)
        self.assert_((zero_filled['A'][:5] == 0).all())

        padded = self.tsframe.fill(method='pad')
        self.assert_(np.isnan(padded['A'][:5]).all())
        self.assert_((padded['A'][-5:] == padded['A'][-5]).all())

        # mixed type
        self.mixed_frame['foo'][5:20] = np.NaN
        self.mixed_frame['A'][-10:] = np.NaN

        result = self.mixed_frame.fill(value=0)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.tsframe[::3]

        start, end = self.tsframe.index[3], self.tsframe.index[6]

        start_missing = self.tsframe.index[2]
        end_missing = self.tsframe.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_frame_equal(truncated, expected)

    def test_getXS(self):
        idx = self.frame.index[5]
        xs = self.frame.getXS(idx)
        for item, value in xs.iteritems():
            if np.isnan(value):
                self.assert_(np.isnan(self.frame[item][idx]))
            else:
                self.assertEqual(value, self.frame[item][idx])

        # mixed-type getXS
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        frame = self.klass(test_data)
        xs = frame.getXS('1')
        self.assert_(xs.dtype == np.object_)
        self.assertEqual(xs['A'], 1)
        self.assertEqual(xs['B'], '1')

        self.assertRaises(Exception, self.tsframe.getXS,
                          self.tsframe.index[0] - datetools.bday)

    def test_pivot(self):
        data = {
            'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values' : [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data)
        pivoted = frame.pivot(index='index', columns='columns', values='values')

        expected = DataFrame({
            'One' : {'A' : 1., 'B' : 2., 'C' : 3.},
            'Two' : {'A' : 1., 'B' : 2., 'C' : 3.}
        })

        assert_frame_equal(pivoted, expected)

        # corner cases

    def test_reindex(self):
        newFrame = self.frame.reindex(self.ts1.index)

        for col in newFrame.cols():
            for idx, val in newFrame[col].iteritems():
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assert_(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(np.isnan(val))

        for col, series in newFrame.iteritems():
            self.assert_(common.equalContents(series.index, newFrame.index))
        emptyFrame = self.frame.reindex(Index([]))
        self.assert_(len(emptyFrame.index) == 0)

        # Cython code should be unit-tested directly
        nonContigFrame = self.frame.reindex(self.ts1.index[::2])

        for col in nonContigFrame.cols():
            for idx, val in nonContigFrame[col].iteritems():
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assert_(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(np.isnan(val))

        for col, series in nonContigFrame.iteritems():
            self.assert_(common.equalContents(series.index, nonContigFrame.index))

        # corner cases

        # Same index, copies values
        newFrame = self.frame.reindex(self.frame.index)
        self.assert_(newFrame.index is self.frame.index)

        # length zero
        newFrame = self.frame.reindex([])
        self.assert_(not newFrame)
        self.assertEqual(len(newFrame.cols()), len(self.frame.cols()))

        # pass non-Index
        newFrame = self.frame.reindex(list(self.ts1.index))
        self.assert_(newFrame.index.equals(self.ts1.index))

    def test_reindex_columns(self):
        newFrame = self.frame.reindex(columns=['A', 'B', 'E'])

        assert_series_equal(newFrame['B'], self.frame['B'])
        self.assert_(np.isnan(newFrame['E']).all())
        self.assert_('C' not in newFrame)

        # length zero
        newFrame = self.frame.reindex(columns=[])
        self.assert_(not newFrame)

    def test_reindex_mixed(self):
        pass

    def test_transpose(self):
        frame = self.frame
        dft = frame.T
        for idx, series in dft.iteritems():
            for col, value in series.iteritems():
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

        # mixed type
        index, data = common.getMixedTypeDict()
        mixed = self.klass(data, index=index)

        mixed_T = mixed.T
        for col, s in mixed_T.iteritems():
            self.assert_(s.dtype == np.object_)

    def test_diff(self):
        the_diff = self.tsframe.diff(1)

        assert_series_equal(the_diff['A'],
                            self.tsframe['A'] - self.tsframe['A'].shift(1))

    def test_shift(self):
        # naive shift
        shiftedFrame = self.tsframe.shift(5)
        self.assert_(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        shiftedFrame = self.tsframe.shift(-5)
        self.assert_(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(-5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        # shift by 0
        unshifted = self.tsframe.shift(0)
        assert_frame_equal(unshifted, self.tsframe)

        # shift by DateOffset
        shiftedFrame = self.tsframe.shift(5, offset=datetools.BDay())
        self.assert_(len(shiftedFrame) == len(self.tsframe))

        shiftedFrame2 = self.tsframe.shift(5, timeRule='WEEKDAY')
        assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = self.tsframe.index[0]
        shifted_d = d + datetools.BDay(5)
        assert_series_equal(self.tsframe.getXS(d),
                            shiftedFrame.getXS(shifted_d))

    def test_apply(self):
        # ufunc
        applied = self.frame.apply(np.sqrt)
        assert_series_equal(np.sqrt(self.frame['A']), applied['A'])

        # aggregator
        applied = self.frame.apply(np.mean)
        self.assertEqual(applied['A'], np.mean(self.frame['A']))

        d = self.frame.index[0]
        applied = self.frame.apply(np.mean, axis=1)
        self.assertEqual(applied[d], np.mean(self.frame.getXS(d)))
        self.assert_(applied.index is self.frame.index) # want this

        # empty
        applied = self.empty.apply(np.sqrt)
        self.assert_(not applied)

        applied = self.empty.apply(np.mean)
        self.assert_(not applied)

    def test_tapply(self):
        d = self.frame.index[0]
        tapplied = self.frame.tapply(np.mean)
        self.assertEqual(tapplied[d], np.mean(self.frame.getXS(d)))

    def test_applymap(self):
        applied = self.frame.applymap(lambda x: x * 2)
        assert_frame_equal(applied, self.frame * 2)

        result = self.frame.applymap(type)

    def test_groupby(self):
        grouped = self.tsframe.groupby(lambda x: x.weekday())

        # aggregate
        aggregated = grouped.aggregate(np.mean)
        self.assertEqual(len(aggregated), 5)
        self.assertEqual(len(aggregated.cols()), 4)

        # by string
        tscopy = self.tsframe.copy()
        tscopy['weekday'] = [x.weekday() for x in tscopy.index]
        stragged = tscopy.groupby('weekday').aggregate(np.mean)
        del stragged['weekday']
        assert_frame_equal(stragged, aggregated)

        # transform
        transformed = grouped.transform(lambda x: x - x.mean())
        self.assertEqual(len(aggregated), 5)
        self.assertEqual(len(aggregated.cols()), 4)

        # iterate
        for weekday, group in grouped:
            self.assert_(group.index[0].weekday() == weekday)

    def test_tgroupby(self):
        grouping = {
            'A' : 0,
            'B' : 1,
            'C' : 0,
            'D' : 1
        }

        grouped = self.frame.tgroupby(grouping.get, np.mean)
        self.assertEqual(len(grouped), len(self.frame.index))
        self.assertEqual(len(grouped.cols()), 2)

    def test_filter(self):
        # items

        filtered = self.frame.filter(['A', 'B', 'E'])
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('E' not in filtered)

        # like
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        filtered = fcopy.filter(like='A')
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('AA' in filtered)

        # regex
        filterd = fcopy.filter(regex='[A]+')
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('AA' in filtered)

        # pass in None
        self.assertRaises(Exception, self.frame.filter, items=None)

    def test_sortUp(self):
        # what to do?
        sorted = self.frame.sortUp()

        sorted_A = self.frame.sortUp(column='A')

    def test_sortDown(self):
        sorted = self.frame.sortDown()

        sorted_A = self.frame.sortDown(column='A')

    def test_combineFirst(self):
        # disjoint
        head, tail = self.frame[:5], self.frame[5:]

        combined = head.combineFirst(tail)
        reordered_frame = self.frame.reindex(combined.index)
        assert_frame_equal(combined, reordered_frame)
        self.assert_(common.equalContents(combined.cols(), self.frame.cols()))
        assert_series_equal(combined['A'], reordered_frame['A'])

        # same index
        fcopy = self.frame.copy()
        fcopy['A'] = 1
        del fcopy['C']

        fcopy2 = self.frame.copy()
        fcopy2['B'] = 0
        del fcopy2['D']

        combined = fcopy.combineFirst(fcopy2)

        self.assert_((combined['A'] == 1).all())
        assert_series_equal(combined['B'], fcopy['B'])
        assert_series_equal(combined['C'], fcopy2['C'])
        assert_series_equal(combined['D'], fcopy['D'])

        # overlap
        head, tail = reordered_frame[:10].copy(), reordered_frame
        head['A'] = 1

        combined = head.combineFirst(tail)
        self.assert_((combined['A'][:10] == 1).all())

        # reverse overlap
        tail['A'][:10] = 0
        combined = tail.combineFirst(head)
        self.assert_((combined['A'][:10] == 0).all())

        # no overlap
        f = self.frame[:10]
        g = self.frame[10:]
        combined = f.combineFirst(g)
        assert_series_equal(combined['A'].reindex(f.index), f['A'])
        assert_series_equal(combined['A'].reindex(g.index), g['A'])

        # corner cases
        comb = self.frame.combineFirst(self.empty)
        self.assert_(comb is self.frame)

        comb = self.empty.combineFirst(self.frame)
        self.assert_(comb is self.frame)

    def test_combineAdd(self):
        # trivial
        comb = self.frame.combineAdd(self.frame)

        assert_frame_equal(comb, self.frame * 2)

        # corner cases
        comb = self.frame.combineAdd(self.empty)
        self.assert_(comb is self.frame)

        comb = self.empty.combineAdd(self.frame)
        self.assert_(comb is self.frame)

    def test_combineMult(self):
        # trivial
        comb = self.frame.combineMult(self.frame)

        assert_frame_equal(comb, self.frame ** 2)

        # corner cases
        comb = self.frame.combineMult(self.empty)
        self.assert_(comb is self.frame)

        comb = self.empty.combineMult(self.frame)
        self.assert_(comb is self.frame)

    def test_join_index(self):
        # left / right

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2)
        self.assert_(f.index.equals(joined.index))
        self.assertEqual(len(joined.cols()), 4)

        joined = f.join(f2, how='left')
        self.assert_(joined.index.equals(f.index))
        self.assertEqual(len(joined.cols()), 4)

        joined = f.join(f2, how='right')
        self.assert_(joined.index.equals(f2.index))
        self.assertEqual(len(joined.cols()), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='left')

        # inner

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='inner')
        self.assert_(joined.index.equals(f.index.intersection(f2.index)))
        self.assertEqual(len(joined.cols()), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='inner')

        # outer

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='outer')
        self.assert_(common.equalContents(self.frame.index, joined.index))
        self.assertEqual(len(joined.cols()), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='outer')

        self.assertRaises(Exception, f.join, f2, how='foo')

    def test_join(self):
        index, data = common.getMixedTypeDict()
        target = self.klass(data, index=index)

        # Join on string value
        source = self.klass({'MergedA' : data['A'], 'MergedD' : data['D']},
                            index=data['C'])
        merged = target.join(source, on='C')

        self.assert_(np.array_equal(merged['MergedA'], target['A']))
        self.assert_(np.array_equal(merged['MergedD'], target['D']))

        # Test when some are missing

        # merge column not p resent
        self.assertRaises(Exception, target.join, source, on='E')

        # corner cases

        # nothing to merge
        merged = target.join(source.reindex([]), on='C')

        # overlap
        source_copy = source.copy()
        source_copy['A'] = 0
        self.assertRaises(Exception, target.join, source_copy, on='A')

        # can't specify how
        self.assertRaises(Exception, target.join, source, on='C',
                          how='left')

    def test_statistics(self):
        sumFrame = self.frame.apply(np.sum)
        for col, series in self.frame.iteritems():
            self.assertEqual(sumFrame[col], series.sum())

    def _check_statistic(self, frame, name, alternative):
        f = getattr(frame, name)

        result = f(axis=0)
        assert_series_equal(result, frame.apply(alternative))

        result = f(axis=1)
        comp = frame.apply(alternative, axis=1).reindex(result.index)
        assert_series_equal(result, comp)

        self.assertRaises(Exception, f, axis=2)

    def test_count(self):
        f = lambda s: notnull(s).sum()

        self._check_statistic(self.frame, 'count', f)

    def test_sum(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].sum()

        self._check_statistic(self.frame, 'sum', f)

    def test_sum_object(self):
        values = self.frame.values.astype(int)
        frame = self.klass(values, index=self.frame.index,
                           columns=self.frame.cols())
        deltas = frame * timedelta(1)
        deltas.sum()

    def test_product(self):
        def f(x):
            x = np.asarray(x)
            return np.prod(x[notnull(x)])

        self._check_statistic(self.frame, 'product', f)

    def test_mean(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].mean()

        self._check_statistic(self.frame, 'mean', f)

    def test_median(self):
        def f(x):
            x = np.asarray(x)
            return np.median(x[notnull(x)])

        self._check_statistic(self.frame, 'median', f)

    def test_min(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].min()

        self._check_statistic(self.frame, 'min', f)

    def test_max(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].max()

        self._check_statistic(self.frame, 'max', f)

    def test_mad(self):
        f = lambda x: np.abs(x - x.mean()).mean()

        self._check_statistic(self.frame, 'mad', f)

    def test_var(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].var(ddof=1)

        self._check_statistic(self.frame, 'var', f)

    def test_std(self):
        def f(x):
            x = np.asarray(x)
            return x[notnull(x)].std(ddof=1)

        self._check_statistic(self.frame, 'std', f)

    def test_skew(self):
        try:
            from scipy.stats import skew
        except ImportError:
            return

        def f(x):
            x = np.asarray(x)
            return skew(x[notnull(x)], bias=False)

        self._check_statistic(self.frame, 'skew', f)

    def test_cumsum(self):
        cumsum = self.tsframe.cumsum()

        assert_series_equal(cumsum['A'], np.cumsum(self.tsframe['A'].fill(0)))

        df = DataFrame({'A' : np.arange(20)}, index=np.arange(20))

        # works
        result = df.cumsum()

    def test_cumprod(self):
        cumprod = self.tsframe.cumprod()

        assert_series_equal(cumprod['A'],
                            np.cumprod(self.tsframe['A'].fill(1)))

if __name__ == '__main__':
    unittest.main()
