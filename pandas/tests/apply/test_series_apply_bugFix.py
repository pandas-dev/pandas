import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    concat,
    date_range,
    timedelta_range,
)
import pandas._testing as tm
from pandas.tests.apply.common import series_transform_kernels

def test_series_map_NAinteger():
    s = pd.Series([1,2,None],dtype="Int32")

    def increment(x):
        if x is None:
            return pd.NA
        return x+1

       
    result = s.map(increment)

    expectedResult = pd.Series([2,3,pd.NA],dtype = "Int32")

    pd.testing.assert_series_equal(result,expectedResult)

   

def test_series_apply_NAinteger():
    s = pd.Series([1,2,None],dtype="Int32")

    def increment(x):
        if x is None:
            return pd.NA
        return x+1

       
    result = s.apply(increment)

    expectedResult = pd.Series([2,3,pd.NA],dtype = "Int32")

    pd.testing.assert_series_equal(result,expectedResult)