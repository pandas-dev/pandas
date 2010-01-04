# pylint: disable-msg=W0612

from copy import deepcopy
from cStringIO import StringIO
import os
import unittest

from numpy import random
import numpy as np

from pandas.core.api import DataFrame, Index, Series, notnull
import pandas.core.datetools as datetools

from pandas.core.tests.common import (assert_almost_equal,
                                      assert_series_equal,
                                      assert_frame_equal,
                                      randn)

import pandas.core.tests.common as common

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

    def test_constructor(self):
        df = self.klass()
        self.assert_(len(df.index) == 0)

        df = self.klass(data={})
        self.assert_(len(df.index) == 0)

    def test_constructor_mixed(self):
        index, data = common.getMixedTypeDict()

        indexed_frame = self.klass(data, index=index)
        unindexed_frame = self.klass(data)

    def test_fromDict(self):
        frame = self.klass.fromDict(col1=self.ts1, col2 = self.ts2)

        common.assert_dict_equal(self.ts1, frame['col1'], compare_keys=False)
        common.assert_dict_equal(self.ts2, frame['col2'], compare_keys=False)

        # cast float tests
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        frame = self.klass.fromDict(test_data)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.float_)
        self.assert_(frame['A'].dtype == np.float_)

        frame = self.klass.fromDict(test_data, castFloat=False)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.object_)
        self.assert_(frame['A'].dtype == np.float_)

        # can't cast to float
        test_data = {
                'A' : dict(zip(range(20), common.makeDateIndex(20))),
                'B' : dict(zip(range(15), randn(15)))
        }
        frame = self.klass.fromDict(test_data)
        self.assertEqual(len(frame), 20)
        self.assert_(frame['A'].dtype == np.object_)
        self.assert_(frame['B'].dtype == np.float_)

        # Corner cases
        self.assertEqual(len(self.klass.fromDict({})), 0)
        self.assertEqual(len(self.klass.fromDict()), 0)
        self.assertRaises(Exception, self.klass.fromDict, [self.ts1, self.ts2])

        # Length-one dict micro-optimization
        frame = self.klass.fromDict({'A' : {'1' : 1, '2' : 2}})
        self.assert_(np.array_equal(frame.index, ['1', '2']))

    def test_toDict(self):
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        recons_data = self.klass.fromDict(test_data, castFloat=False).toDict()

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

    def test_fromMatrix(self):
        mat = np.zeros((2, 3), dtype=float)

        frame = self.klass.fromMatrix(mat, ['A', 'B', 'C'], [1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.cols()), 3)

        self.assertRaises(Exception, self.klass.fromMatrix,
                          mat, ['A', 'B', 'C'], [1])
        self.assertRaises(Exception, self.klass.fromMatrix,
                          mat, ['A', 'B'], [1, 2])

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
        biggie = self.klass.fromMatrix(np.zeros((1000, 4)),
                                       range(4), range(1000))
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

        buf = StringIO()
        biggie.toString(buffer=buf)

    def test_toString(self):
        pass

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
        f = -self.frame

    def test_firstTimeWithNValues(self):
        self.frame['A'][:5] = np.NaN

        index = self.frame._firstTimeWithNValues()
        self.assert_(index == self.frame.index[5])

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

    def test_combineSeries(self):

        # Series
        series = self.frame.getXS(self.frame.index[0])

        added = self.frame + series

        for key, s in added.iteritems():
            assert_series_equal(s, self.frame[key] + series[key])

        larger_series = series.toDict()
        larger_series['E'] = 1
        larger_series = Series.fromDict(larger_series)
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

    def test_combineFunc(self):
        pass

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
        pass

    def test_cols(self):
        pass

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

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(datetools.bmonthEnd)
        rule_monthly = self.tsframe.asfreq('EOM')

        assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('WEEKDAY', fillMethod='pad')

    def test_asMatrix(self):
        frame = self.frame
        mat = frame.asMatrix()
        smallerCols = ['C', 'A']
        # smallerMat = frame.asMatrix(smallerCols)
        frameCols = frame.cols()
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])

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

    def test_getTS(self):
        frame = self.tsframe

        tsFrame = frame.getTS(fromDate=frame.index[5], nPeriods=5)
        assert_frame_equal(tsFrame, frame[5:10])

        tsFrame = frame.getTS(fromDate=frame.index[5], toDate=frame.index[9])
        assert_frame_equal(tsFrame, frame[5:10])

        tsFrame = frame.getTS(nPeriods=5, toDate=frame.index[9])
        assert_frame_equal(tsFrame, frame[5:10])

        A = frame.getTS(colName='A', nPeriods=5, toDate=frame.index[9])
        assert_series_equal(A, frame['A'][5:10])

        self.assertRaises(Exception, frame.getTS, nPeriods=5)

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

    def test_pivot(self):
        pass

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

        # pass non-Index
        newFrame = self.frame.reindex(list(self.ts1.index))
        self.assert_(newFrame.index.equals(self.ts1.index))

    def test_reindex_columns(self):
        newFrame = self.frame.reindex(columns=['A', 'B', 'E'])

        assert_series_equal(newFrame['B'], self.frame['B'])
        self.assert_(np.isnan(newFrame['E']).all())
        self.assert_('C' not in newFrame)

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
        pass

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

        # shift by DateOffset
        shiftedFrame = self.tsframe.shift(5, offset=datetools.BDay())
        self.assert_(len(shiftedFrame) == len(self.tsframe))

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

    def test_tapply(self):
        d = self.frame.index[0]
        tapplied = self.frame.tapply(np.mean)
        self.assertEqual(tapplied[d], np.mean(self.frame.getXS(d)))

    def test_applymap(self):
        f = lambda x: x * 2
        applied = self.frame.applymap(f)

        assert_frame_equal(applied, self.frame * 2)

    def test_groupby(self):
        grouped = self.tsframe.groupby(lambda x: x.weekday())

        # aggregate
        aggregated = grouped.aggregate(np.mean)
        self.assertEqual(len(aggregated), 5)
        self.assertEqual(len(aggregated.cols()), 4)

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

        filtered = self.frame.filterItems(['A', 'B', 'E'])
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('E' not in filtered)

        # like
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        filtered = fcopy.filterLike('A')
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('AA' in filtered)

        # regex
        filterd = fcopy.filter(regex='[A]+')
        self.assertEqual(len(filtered.cols()), 2)
        self.assert_('AA' in filtered)

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

    def test_leftJoin(self):
        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.leftJoin(f2)
        self.assert_(f.index.equals(joined.index))
        self.assertEqual(len(joined.cols()), 4)

        # corner case
        self.assertRaises(Exception, self.frame.leftJoin, self.frame)

    def test_outerJoin(self):
        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.outerJoin(f2)
        self.assert_(common.equalContents(self.frame.index, joined.index))
        self.assertEqual(len(joined.cols()), 4)

        # corner case
        self.assertRaises(Exception, self.frame.outerJoin, self.frame)

    def test_merge(self):
        index, data = common.getMixedTypeDict()
        target = self.klass(data, index=index)

        # Merge on string value
        source = self.klass({'MergedA' : data['A'], 'MergedD' : data['D']},
                            index=data['C'])
        merged = target.merge(source, on='C')

        self.assert_(np.array_equal(merged['MergedA'], target['A']))
        self.assert_(np.array_equal(merged['MergedD'], target['D']))

        # Test when some are missing

        # corner case

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
        from scipy.stats import skew
        def f(x):
            x = np.asarray(x)
            return skew(x[notnull(x)], bias=False)

        self._check_statistic(self.frame, 'skew', f)

    def test_cumsum(self):
        cumsum = self.tsframe.cumsum()

        assert_series_equal(cumsum['A'], np.cumsum(self.tsframe['A'].fill(0)))

if __name__ == '__main__':
    unittest.main()
