# pylint: disable-msg=W0612

from copy import deepcopy
from datetime import datetime
import unittest

from numpy import random
import numpy as np

from pandas.core.api import DateRange, DataFrame, Index, Series
from pandas.core.datetools import bday

from pandas.core.tests.common import assert_almost_equal
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

    def test_constructor(self):

        self.assertRaises(Exception, DataFrame)

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
                'B' : dict(zip(range(15), common.randn(15)))
        }
        frame = self.klass.fromDict(test_data)
        self.assertEqual(len(frame), 20)
        self.assert_(frame['A'].dtype == np.object_)
        self.assert_(frame['B'].dtype == np.float_)

        # Corner cases
        self.assertEqual(len(self.klass.fromDict({})), 0)
        self.assertEqual(len(self.klass.fromDict()), 0)
        self.assertRaises(Exception, self.klass.fromDict, [self.ts1, self.ts2])

    def test_toDict(self):
        pass

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
        pass

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
        biggie = self.klass({'A' : common.randn(1000),
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
                          common.randn(len(self.frame) + 1))


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
        pass

    def test_len(self):
        pass

    def test_contains(self):
        pass

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
        pass

    def test_firstTimeWithValue(self):
        pass

    def test_lastTimeWithValue(self):
        pass

    def test_combineFrame(self):
        pass

    def test_combineSeries(self):
        pass

    def test_combineFunc(self):
        pass

    def test_toCSV(self):
        pass

    def test_toDataMatrix(self):
        pass

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
        pass

    def test_asfreq(self):
        pass

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
        pass

    def test_dropEmptyRows(self):
        pass

    def test_dropIncompleteRows(self):
        pass

    def test_fill(self):
        pass

    def test_getTS(self):
        frame = self.tsframe
        tsFrame = frame.getTS(fromDate=frame.index[5], nPeriods=5)
        for i, idx in enumerate(tsFrame.index):
            self.assertEqual(idx, frame.index[5+i])
            for col, series in tsFrame.iteritems():
                self.assertEqual(idx, series.index[i])
        for col, series in frame.iteritems():
            for idx, value in series.iteritems():
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

    def test_truncate(self):
        pass

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
        pass

    def test_sortDown(self):
        pass

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
