import numpy as np
from numpy.random import randint

import nose

from pandas import DataFrame
from pandas import read_clipboard
from pandas import get_option
from pandas.util import testing as tm
from pandas.util.testing import makeCustomDataframe as mkdf


try:
    import pandas.util.clipboard
except OSError:
    raise nose.SkipTest("no clipboard found")


class TestClipboard(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestClipboard, cls).setUpClass()
        cls.data = {}
        cls.data['string'] = mkdf(5, 3, c_idx_type='s', r_idx_type='i',
                                  c_idx_names=[None], r_idx_names=[None])
        cls.data['int'] = mkdf(5, 3, data_gen_f=lambda *args: randint(2),
                               c_idx_type='s', r_idx_type='i',
                               c_idx_names=[None], r_idx_names=[None])
        cls.data['float'] = mkdf(5, 3,
                                 data_gen_f=lambda r, c: float(r) + 0.01,
                                 c_idx_type='s', r_idx_type='i',
                                 c_idx_names=[None], r_idx_names=[None])
        cls.data['mixed'] = DataFrame({'a': np.arange(1.0, 6.0) + 0.01,
                                       'b': np.arange(1, 6),
                                       'c': list('abcde')})
        # Test GH-5346
        max_rows = get_option('display.max_rows')
        cls.data['longdf'] = mkdf(max_rows+1, 3, data_gen_f=lambda *args: randint(2),
                                  c_idx_type='s', r_idx_type='i',
                                  c_idx_names=[None], r_idx_names=[None])  
        cls.data_types = list(cls.data.keys())

    @classmethod
    def tearDownClass(cls):
        super(TestClipboard, cls).tearDownClass()
        del cls.data_types, cls.data

    def check_round_trip_frame(self, data_type, excel=None, sep=None):
        data = self.data[data_type]
        data.to_clipboard(excel=excel, sep=sep)
        if sep is not None:
            result = read_clipboard(sep=sep,index_col=0)
        else:
            result = read_clipboard()
        tm.assert_frame_equal(data, result, check_dtype=False)

    def test_round_trip_frame_sep(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt,sep=',')

    def test_round_trip_frame_string(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt,excel=False)

    def test_round_trip_frame(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt)
