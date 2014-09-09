# pylint: disable-msg=W0612,E1101,W0141
import datetime
import nose

from numpy.random import randn
import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas import Panel, DataFrame, Series, notnull, isnull

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)
import pandas.core.common as com
import pandas.util.testing as tm
from pandas.compat import (range, lrange, StringIO, lzip, u,
                                product as cart_product, zip)
import pandas as pd

import pandas.index as _index


class TestNewGroupByAttr(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.frame = DataFrame(np.random.randn(10, 3), index=index,
                               columns=Index(['A', 'B', 'C'], name='exp'))

        self.frame.ix[1, [1, 2]] = np.nan
        self.frame.ix[7, [0, 1]] = np.nan

    AGG_FUNCTIONS = ['skew', 'mad']

    def test_newattr(self) :
        for op, level, axis, skipna in cart_product(self.AGG_FUNCTIONS,
                                                    lrange(2), lrange(2),
                                                    [True,False]) :
            if axis == 0 :
                frame = self.frame
            else :
                frame = self.frame.T

            grouped = frame.groupby(level=level,axis=axis)
            result = getattr(grouped,op)(skipna=skipna)
            expected = getattr(frame,op)(level=level,axis=axis,skipna=skipna)
            assert_frame_equal(result, expected)

if __name__ == '__main__':

    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
