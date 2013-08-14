import unittest

import numpy as np
import pandas as pd

class TestPeriodField(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_get_period_field_raises_on_out_of_range(self):
        self.assertRaises(ValueError, pd.tslib.get_period_field, -1, 0, 0)

    def test_get_period_field_array_raises_on_out_of_range(self):
        self.assertRaises(ValueError, pd.tslib.get_period_field_arr, -1, np.empty(1), 0)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
