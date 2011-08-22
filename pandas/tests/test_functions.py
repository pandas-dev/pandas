from pandas import Index, isnull
import pandas.sandbox.functions as fns
import numpy as np
import unittest
import numpy as np


class TestReductions(unittest.TestCase):

    def setUp(self):
        self.index = Index(np.arange(100))
        pass

    def test_upsample_mean(self):
        pass

    def test_upsample_max(self):
        pass

    def test_upsample_min(self):
        pass

    def test_upsample_generic(self):
        pass

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

