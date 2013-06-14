# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

from datetime import datetime, timedelta
import operator
import pickle
import unittest
import nose
import os

import numpy as np
import pandas.util.testing as tm
import pandas as pd
from pandas import Index
from pandas.sparse.tests import test_sparse
from pandas.util import py3compat
from pandas.util.misc import is_little_endian

class TestPickle(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_pickles import create_data
        self.data = create_data()

    def compare(self, vf):

        # py3 compat when reading py2 pickle
        
        try:
            with open(vf,'rb') as fh:
                data = pickle.load(fh)
        except (ValueError):

            # we are trying to read a py3 pickle in py2.....
            return
        except:
            if not py3compat.PY3:
                raise
            with open(vf,'rb') as fh:
                data = pickle.load(fh, encoding='latin1')

        for typ, dv in data.items():
            for dt, result in dv.items():

                expected = self.data[typ][dt]

                if isinstance(expected,Index):
                    self.assert_(expected.equals(result))
                    continue

                if typ.startswith('sp_'):
                    comparator = getattr(test_sparse,"assert_%s_equal" % typ)
                    comparator(result,expected,exact_indices=False)
                else:
                    comparator = getattr(tm,"assert_%s_equal" % typ)
                    comparator(result,expected)

    def test_read_pickles_0_10_1(self):
        if not is_little_endian():
            raise nose.SkipTest("known failure of test_read_pickles_0_10_1 on non-little endian")

        pth = tm.get_data_path('legacy_pickle/0.10.1')
        for f in os.listdir(pth):
            vf = os.path.join(pth,f)
            self.compare(vf)

    def test_read_pickles_0_11_0(self):
        if not is_little_endian():
            raise nose.SkipTest("known failure of test_read_pickles_0_11_0 on non-little endian")

        pth = tm.get_data_path('legacy_pickle/0.11.0')
        for f in os.listdir(pth):
            vf = os.path.join(pth,f)
            self.compare(vf)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
