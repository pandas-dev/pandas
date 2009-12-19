# pylint: disable-msg=W0612

from copy import deepcopy
import os
import unittest

from numpy import random
import numpy as np

from pandas.core.api import DateRange, DataFrame, Index, Series
from pandas.core.datetools import bday
import pandas.core.datetools as datetools

from pandas.core.tests.common import assert_almost_equal, randn
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
        self.empty = DataFrame({})

    def test_constructor(self):
        df = DataFrame()
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

    def test_repr(self):
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
        foo = repr(biggie)

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

        # set value
        self.frame['col7'] = 5
        assert((self.frame['col7'] == 5).all())

        self.frame['col8'] = 'foo'
        assert((self.frame['col8'] == 'foo').all())

        self.assertRaises(Exception, self.frame.__setitem__,
                          randn(len(self.frame) + 1))

        # set ndarray
        arr = randn(len(self.frame))
        self.frame['col9'] = arr
        self.assert_((self.frame['col9'] == arr).all())

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

        frame = DataFrame({'foo' : mat}, index=self.frame.index)
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
        self.assert_(self_added.index is self.frame.index)

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
        pass

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

    def test_toString(self):
        pass

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

        frame = DataFrame({'foo' : mat}, index=self.frame.index)

        smaller_frame = frame.dropEmptyRows()
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = np.NaN

        frame = DataFrame({'foo' : mat}, index=self.frame.index)
        frame['bar'] = 5

        smaller_frame = frame.dropIncompleteRows()
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        samesize_frame = frame.dropIncompleteRows(specificColumns=['bar'])
        self.assert_(samesize_frame.index is self.frame.index)

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
        common.assert_frame_equal(tsFrame, frame[5:10])

        tsFrame = frame.getTS(fromDate=frame.index[5], toDate=frame.index[9])
        common.assert_frame_equal(tsFrame, frame[5:10])

        tsFrame = frame.getTS(nPeriods=5, toDate=frame.index[9])
        common.assert_frame_equal(tsFrame, frame[5:10])

        A = frame.getTS(colName='A', nPeriods=5, toDate=frame.index[9])
        common.assert_series_equal(A, frame['A'][5:10])

        self.assertRaises(Exception, frame.getTS, nPeriods=5)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.tsframe[::3]

        start, end = self.tsframe.index[3], self.tsframe.index[6]

        start_missing = self.tsframe.index[2]
        end_missing = self.tsframe.index[7]

        # neither specified
        truncated = ts.truncate()
        common.assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        common.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        common.assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        common.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        common.assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        common.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        common.assert_frame_equal(truncated, expected)

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

    def test_diff(self):
        pass

    def test_shift(self):
        shiftedFrame = self.tsframe.shift(5)
        for i, idx in enumerate(shiftedFrame.index):
            self.assert_(idx-5*bday == self.tsframe.index[i])
        series = shiftedFrame['A']
        for i, idx in enumerate(series.index):
            self.assert_(idx-5*bday == self.tsframe.index[i])

    def test_apply(self):
        pass

    def test_tapply(self):
        pass

    def test_applymap(self):
        pass

    def test_tgroupby(self):
        pass

    def test_filterItems(self):
        pass

    def test_sortUp(self):
        # what to do?
        sorted = self.frame.sortUp()

        sorted_A = self.frame.sortUp(column='A')

    def test_sortDown(self):
        sorted = self.frame.sortDown()

        sorted_A = self.frame.sortDown(column='A')

    def test_filterLike(self):
        pass

    def test_combineFirst(self):
        pass

    def test_combineAdd(self):
        pass

    def test_combineMult(self):
        pass

    def test_outerJoin(self):
        pass

    def test_leftJoin(self):
        pass

    def test_merge(self):
        index, data = common.getMixedTypeDict()
        target = DataFrame(data, index=index)

        # Merge on string value
        source = DataFrame({'MergedA' : data['A'], 'MergedD' : data['D']},
                           index=data['C'])
        merged = target.merge(source, on='C')

        self.assert_(np.array_equal(merged['MergedA'], target['A']))
        self.assert_(np.array_equal(merged['MergedD'], target['D']))

        # Test when some are missing

    def test_statistics(self):
        sumFrame = self.frame.apply(np.sum)
        for col, series in self.frame.iteritems():
            self.assertEqual(sumFrame[col], series.sum())

    def test_count(self):
        pass

    def test_sum(self):
        pass

    def test_product(self):
        pass

    def test_mean(self):
        pass

    def test_median(self):
        pass

    def test_min(self):
        pass

    def test_max(self):
        pass

    def test_mad(self):
        pass

    def test_var(self):
        pass

    def test_std(self):
        pass

    def test_skew(self):
        pass

    def test_withColumns(self):
        pass

    def testGroupBy(self):
        pass

if __name__ == '__main__':
    unittest.main()
