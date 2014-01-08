# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

from datetime import datetime, timedelta
import operator
import pickle as pkl
import nose
import os

import numpy as np
import pandas.util.testing as tm
import pandas as pd
from pandas import Index
from pandas.sparse.tests import test_sparse
from pandas import compat
from pandas.compat import u
from pandas.util.misc import is_little_endian
import pandas

def _read_pickle(vf, encoding=None, compat=False):
    from pandas.compat import pickle_compat as pc
    with open(vf,'rb') as fh:
        pc.load(fh, encoding=encoding, compat=compat)

class TestPickle(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_pickles import create_data
        self.data = create_data()
        self.path = u('__%s__.pickle' % tm.rands(10))

    def compare_element(self, typ, result, expected):
        if isinstance(expected,Index):
            self.assert_(expected.equals(result))
            return

        if typ.startswith('sp_'):
            comparator = getattr(test_sparse,"assert_%s_equal" % typ)
            comparator(result,expected,exact_indices=False)
        else:
            comparator = getattr(tm,"assert_%s_equal" % typ)
            comparator(result,expected)

    def compare(self, vf):

        # py3 compat when reading py2 pickle
        try:
            data = pandas.read_pickle(vf)
        except (ValueError) as detail:
            # trying to read a py3 pickle in py2
            return

        for typ, dv in data.items():
            for dt, result in dv.items():
                try:
                    expected = self.data[typ][dt]
                except (KeyError):
                    continue

                self.compare_element(typ, result, expected)

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

    def test_round_trip_current(self):

        for typ, dv in self.data.items():

            for dt, expected in dv.items():

                with tm.ensure_clean(self.path) as path:

                    pd.to_pickle(expected,path)

                    result = pd.read_pickle(path)
                    self.compare_element(typ, result, expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
