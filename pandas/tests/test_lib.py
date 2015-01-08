# -*- coding: utf-8 -*-
from datetime import datetime, timedelta, date, time

import numpy as np

import pandas as pd
from pandas.lib import isscalar, item_from_zerodim, max_len_string_array
import pandas.util.testing as tm
from pandas.compat import u

class TestMisc(tm.TestCase):

    def test_max_len_string_array(self):

        arr = np.array(['foo','b',np.nan],dtype='object')
        self.assertTrue(max_len_string_array(arr),3)

        # unicode
        arr = arr.astype('U')
        self.assertTrue(max_len_string_array(arr),3)

class TestIsscalar(tm.TestCase):
    def test_isscalar_builtin_scalars(self):
        self.assertTrue(isscalar(None))
        self.assertTrue(isscalar(True))
        self.assertTrue(isscalar(False))
        self.assertTrue(isscalar(0.))
        self.assertTrue(isscalar(np.nan))
        self.assertTrue(isscalar('foobar'))
        self.assertTrue(isscalar(b'foobar'))
        self.assertTrue(isscalar(u('efoobar')))
        self.assertTrue(isscalar(datetime(2014, 1, 1)))
        self.assertTrue(isscalar(date(2014, 1, 1)))
        self.assertTrue(isscalar(time(12, 0)))
        self.assertTrue(isscalar(timedelta(hours=1)))
        self.assertTrue(isscalar(pd.NaT))

    def test_isscalar_builtin_nonscalars(self):
        self.assertFalse(isscalar({}))
        self.assertFalse(isscalar([]))
        self.assertFalse(isscalar([1]))
        self.assertFalse(isscalar(()))
        self.assertFalse(isscalar((1,)))
        self.assertFalse(isscalar(slice(None)))
        self.assertFalse(isscalar(Ellipsis))

    def test_isscalar_numpy_array_scalars(self):
        self.assertTrue(isscalar(np.int64(1)))
        self.assertTrue(isscalar(np.float64(1.)))
        self.assertTrue(isscalar(np.int32(1)))
        self.assertTrue(isscalar(np.object_('foobar')))
        self.assertTrue(isscalar(np.str_('foobar')))
        self.assertTrue(isscalar(np.unicode_(u('foobar'))))
        self.assertTrue(isscalar(np.bytes_(b'foobar')))
        self.assertTrue(isscalar(np.datetime64('2014-01-01')))
        self.assertTrue(isscalar(np.timedelta64(1, 'h')))

    def test_isscalar_numpy_zerodim_arrays(self):
        for zerodim in [np.array(1),
                        np.array('foobar'),
                        np.array(np.datetime64('2014-01-01')),
                        np.array(np.timedelta64(1, 'h'))]:
            self.assertFalse(isscalar(zerodim))
            self.assertTrue(isscalar(item_from_zerodim(zerodim)))

    def test_isscalar_numpy_arrays(self):
        self.assertFalse(isscalar(np.array([])))
        self.assertFalse(isscalar(np.array([[]])))
        self.assertFalse(isscalar(np.matrix('1; 2')))

    def test_isscalar_pandas_scalars(self):
        self.assertTrue(isscalar(pd.Timestamp('2014-01-01')))
        self.assertTrue(isscalar(pd.Timedelta(hours=1)))
        self.assertTrue(isscalar(pd.Period('2014-01-01')))

    def test_isscalar_pandas_containers(self):
        self.assertFalse(isscalar(pd.Series()))
        self.assertFalse(isscalar(pd.Series([1])))
        self.assertFalse(isscalar(pd.DataFrame()))
        self.assertFalse(isscalar(pd.DataFrame([[1]])))
        self.assertFalse(isscalar(pd.Panel()))
        self.assertFalse(isscalar(pd.Panel([[[1]]])))
        self.assertFalse(isscalar(pd.Index([])))
        self.assertFalse(isscalar(pd.Index([1])))
