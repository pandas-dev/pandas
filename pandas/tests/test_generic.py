# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import operator
import unittest
import nose

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, Panel,
                    isnull, notnull,date_range)
from pandas.core.index import Index, MultiIndex
from pandas.tseries.index import Timestamp, DatetimeIndex

import pandas.core.common as com

from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_panel_equal,
                                 assert_almost_equal,
                                 ensure_clean)
import pandas.util.testing as tm

#------------------------------------------------------------------------------
# Generic types test cases


class Generic(object):

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

    @property
    def _ndim(self):
        return self._typ._AXIS_LEN

    def _axes(self):
        """ return the axes for my object typ """
        return self._typ._AXIS_ORDERS

    def _construct(self, shape, value=None, **kwargs):
        """ construct an object for the given shape
            if value is specified use that if its a scalar
            if value is an array, repeat it as needed """

        if isinstance(shape,int):
            shape = tuple([shape] * self._ndim)
        if value is not None:
            if np.isscalar(value):
                if value == 'empty':
                    arr = None

                    # remove the info axis
                    kwargs.pop(self._typ._info_axis_name,None)
                else:
                    arr = np.empty(shape)
                    arr.fill(value)
            else:
                fshape = np.prod(shape)
                arr = value.ravel()
                new_shape = fshape/arr.shape[0]
                if fshape % arr.shape[0] != 0:
                    raise Exception("invalid value passed in _construct")

                arr = np.repeat(arr,new_shape).reshape(shape)
        else:
            arr = np.random.randn(*shape)
        return self._typ(arr,**kwargs)


    def _compare(self, result, expected):
        self._comparator(result,expected)

    def test_rename(self):

        # single axis
        for axis in self._axes():
            kwargs = { axis : list('ABCD') }
            o = self._construct(4,**kwargs)

            # no values passed
            #self.assertRaises(Exception, o.rename(str.lower))

            # rename a single axis
            result = o.rename(**{ axis : str.lower })
            expected = o.copy()
            setattr(expected,axis,list('abcd'))
            self._compare(result, expected)

        # multiple axes at once

    def test_get_numeric_data(self):

        n = 4
        kwargs = { }
        for i in range(self._ndim):
            kwargs[self._typ._AXIS_NAMES[i]] = list(range(n))

        # get the numeric data
        o = self._construct(n,**kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # non-inclusion
        result = o._get_bool_data()
        expected = self._construct(n,value='empty',**kwargs)
        self._compare(result,expected)

        # get the bool data
        arr = np.array([True,True,False,True])
        o = self._construct(n,value=arr,**kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # _get_numeric_data is includes _get_bool_data, so can't test for non-inclusion

class TestSeries(unittest.TestCase, Generic):
    _typ = Series
    _comparator = lambda self, x, y: assert_series_equal(x,y)

    def test_rename_mi(self):
        s = Series([11,21,31],
                   index=MultiIndex.from_tuples([("A",x) for x in ["a","B","c"]]))
        result = s.rename(str.lower)

class TestDataFrame(unittest.TestCase, Generic):
    _typ = DataFrame
    _comparator = lambda self, x, y: assert_frame_equal(x,y)

    def test_rename_mi(self):
        df = DataFrame([11,21,31],
                       index=MultiIndex.from_tuples([("A",x) for x in ["a","B","c"]]))
        result = df.rename(str.lower)

class TestPanel(unittest.TestCase, Generic):
    _typ = Panel
    _comparator = lambda self, x, y: assert_panel_equal(x,y)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
