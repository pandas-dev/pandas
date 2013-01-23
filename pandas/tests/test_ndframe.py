import unittest

import numpy as np

from pandas.core.generic import NDFrame
import pandas.util.testing as t


class TestNDFrame(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        tdf = t.makeTimeDataFrame()
        self.ndf = NDFrame(tdf._data)

    def test_constructor(self):
        # with cast
        ndf = NDFrame(self.ndf._data, dtype=np.int64)
        self.assert_(ndf.values.dtype == np.int64)

    def test_ndim(self):
        self.assertEquals(self.ndf.ndim, 2)

    def test_astype(self):
        casted = self.ndf.astype(int)
        self.assert_(casted.values.dtype == np.int_)

        casted = self.ndf.astype(np.int32)
        self.assert_(casted.values.dtype == np.int32)

    def test_squeeze(self):
        # noop
        for s in [ t.makeFloatSeries(), t.makeStringSeries(), t.makeObjectSeries() ]:
            t.assert_series_equal(s.squeeze(),s)
        for df in [ t.makeTimeDataFrame() ]:
            t.assert_frame_equal(df.squeeze(),df)
        for p in [ t.makePanel() ]:
            t.assert_panel_equal(p.squeeze(),p)
        for p4d in [ t.makePanel4D() ]:
            t.assert_panel4d_equal(p4d.squeeze(),p4d)

        # squeezing
        df = t.makeTimeDataFrame().reindex(columns=['A'])
        t.assert_series_equal(df.squeeze(),df['A'])

        p = t.makePanel().reindex(items=['ItemA'])
        t.assert_frame_equal(p.squeeze(),p['ItemA'])

        p = t.makePanel().reindex(items=['ItemA'],minor_axis=['A'])
        t.assert_series_equal(p.squeeze(),p.ix['ItemA',:,'A'])

        p4d = t.makePanel4D().reindex(labels=['label1'])
        t.assert_panel_equal(p4d.squeeze(),p4d['label1'])

        p4d = t.makePanel4D().reindex(labels=['label1'],items=['ItemA'])
        t.assert_frame_equal(p4d.squeeze(),p4d.ix['label1','ItemA'])

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
