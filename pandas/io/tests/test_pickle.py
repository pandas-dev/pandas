# pylint: disable=E1101,E1103,W0232

from datetime import datetime, timedelta
import operator
import pickle
import unittest
import nose
import os

import numpy as np
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_almost_equal
from pandas.util import py3compat
import pandas.core.common as com

import pandas.util.testing as tm
import pandas as pd

class TestPickle(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    #def test_hash_error(self):
    #    self.assertRaises(TypeError, hash, self.strIndex)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
