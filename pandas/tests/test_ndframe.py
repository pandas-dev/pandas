import unittest

import numpy as np

from pandas.core.generic import NDFrame
import pandas.util.testing as t

class TestNDFrame(unittest.TestCase):

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
        self.assert_(casted.values.dtype == np.int64)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

