#import unittest

from pandas import DataFrame
from pandas.tools.describe import value_range

import numpy as np

def test_value_range():
    df = DataFrame(np.random.randn(5, 5))
    df.ix[0,2] = -5
    df.ix[2,0] = 5

    res = value_range(df)

    assert( res['Minimum'] == -5 )
    assert( res['Maximum'] == 5 )

    df.ix[0,1] = np.NaN

    assert( res['Minimum'] == -5 )
    assert( res['Maximum'] == 5 )
