# -*- coding: utf-8 -*-
# pylint: disable=E1101,E1103,W0232

import os
import pickle

import pandas as pd

from pandas import Period

import pandas.util.testing as tm

class TestPeriod(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.period = Period('2000Q1')

    def test_round_trip_pickle(self):

        period = Period('2000Q1')
        pickle_path = os.path.join(tm.get_data_path(),
                                   'period.pickle')

        with open(pickle_path, 'wb') as f: pickle.dump(period, f)

        self.assertEqual(self.period, pd.read_pickle(pickle_path))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core']
                    exit=False)