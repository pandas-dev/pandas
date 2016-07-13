# -*- coding: utf-8 -*-

import nose
import numpy as np

from pandas.types.dtypes import DatetimeTZDtype, CategoricalDtype
from pandas.types.common import pandas_dtype

_multiprocess_can_split_ = True


def test_pandas_dtype():

    assert pandas_dtype('datetime64[ns, US/Eastern]') == DatetimeTZDtype(
        'datetime64[ns, US/Eastern]')
    assert pandas_dtype('category') == CategoricalDtype()
    for dtype in ['M8[ns]', 'm8[ns]', 'object', 'float64', 'int64']:
        assert pandas_dtype(dtype) == np.dtype(dtype)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
