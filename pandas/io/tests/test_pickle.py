# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

from datetime import datetime, timedelta
import operator
import pickle as pkl
import unittest
import nose
import os

import numpy as np
import pandas.util.testing as tm
import pandas as pd
from pandas import Index
from pandas.sparse.tests import test_sparse
from pandas import compat
from pandas.util.misc import is_little_endian
import pandas

class TestPickle(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_pickles import create_data
        self.data = create_data()

    def compare(self, vf):

        # py3 compat when reading py2 pickle
        try:
            data = pandas.read_pickle(vf)
        except (ValueError) as detail:
            # trying to read a py3 pickle in py2
            return

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

    def read_pickles(self, version):
        if not is_little_endian():
            raise nose.SkipTest("known failure on non-little endian")

        pth = tm.get_data_path('legacy_pickle/{0}'.format(str(version)))
        for f in os.listdir(pth):
            vf = os.path.join(pth,f)
            self.compare(vf)

    def test_read_pickles_0_10_1(self):
        self.read_pickles('0.10.1')

    def test_read_pickles_0_11_0(self):
        self.read_pickles('0.11.0')

    def test_read_pickles_0_12_0(self):
        self.read_pickles('0.12.0')

    def test_read_pickles_0_13_0(self):
        self.read_pickles('0.13.0')

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
