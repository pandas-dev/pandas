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
        index1 = DateRange(datetime(2008,4,22), periods=50)
        index2 = DateRange(datetime(2008,4,29), periods=50)
        index3 = DateRange(datetime(2008,4,28), periods=50)
        ts1 = Series(random.random(50), index=index1)
        ts2 = Series(random.random(50), index=index2)
        ts3 = Series(random.random(50), index=index3)
        ts4 = Series(random.random(50), index=index1)
        data = {'col1' : ts1,'col2' : ts2,'col3' : ts3, 'col4' : ts4}
        self.frame = self.klass(data=data, index=index3)
        self.ts1 = ts1
        self.ts2 = ts2
        self.ts3 = ts3
        self.ts4 = ts4

    def test_constructor(self):
        pass

    def test_constructor_mixed(self):
        index, data = common.getMixedTypeDict()

        indexed_frame = self.klass(data, index=index)
        unindexed_frame = self.klass(data)

    def test_fromDict(self):
        newFrame = self.klass.fromDict(col1=self.ts1, col2 = self.ts2)
        for idx in newFrame.index:
            if idx in self.ts1.index:
                self.assertEqual(newFrame['col1'][idx], self.ts1[idx])
            if idx in self.ts2.index:
                self.assertEqual(newFrame['col2'][idx], self.ts2[idx])

    def test_fromRecords(self):
        pass

    def test_toRecords(self):
        pass

    def test_fromMatrix(self):
        pass

    def test_nonzero(self):
        pass

    def test_repr(self):
        pass

    def test_getitem(self):
        """Slicing NOT intended for production code"""
        sl = self.frame[:20]
        self.assertEqual(20, len(sl.index))
        for _, series in sl.iteritems():
            self.assertEqual(20, len(series.index))
            self.assert_(common.equalContents(series.index, sl.index))

        for key, _ in self.frame._series.iteritems():
            self.assert_(self.frame[key] is not None)
        self.assert_('random' not in self.frame)

    def test_setitem(self):
        # not sure what else to do here
        series = self.frame['col1']
        self.frame['col5'] = series
        self.assert_('col5' in self.frame)

    def test_delitem(self):
        del self.frame['col1']
        self.assert_('col1' not in self.frame)
        self.assert_('col1' not in self.frame._series)

    def test_pop(self):
        pass

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
        pass

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

    def test_toDict(self):
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
        smallerCols = ['col3', 'col1']
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
        series = cp['col1']
        series[:] = 10
        for idx, value in series.iteritems():
            self.assertNotEqual(self.frame['col1'][idx], value)

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
        frame = self.frame
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
        shiftedFrame = self.frame.shift(5)
        for i, idx in enumerate(shiftedFrame.index):
            self.assert_(idx-5*bday == self.frame.index[i])
        series = shiftedFrame['col1']
        for i, idx in enumerate(series.index):
            self.assert_(idx-5*bday == self.frame.index[i])

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
