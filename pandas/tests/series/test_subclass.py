# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pandas.util.testing as tm


class TestSeriesSubclassing(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_indexing_sliced(self):
        s = tm.SubclassedSeries([1, 2, 3, 4], index=list('abcd'))
        res = s.loc[['a', 'b']]
        exp = tm.SubclassedSeries([1, 2], index=list('ab'))
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = s.iloc[[2, 3]]
        exp = tm.SubclassedSeries([3, 4], index=list('cd'))
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = s.ix[['a', 'b']]
        exp = tm.SubclassedSeries([1, 2], index=list('ab'))
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

    def test_to_frame(self):
        s = tm.SubclassedSeries([1, 2, 3, 4], index=list('abcd'), name='xxx')
        res = s.to_frame()
        exp = tm.SubclassedDataFrame({'xxx': [1, 2, 3, 4]}, index=list('abcd'))
        tm.assert_frame_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedDataFrame)
