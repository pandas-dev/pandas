from datetime import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm


def test_dictmap():
    df = DataFrame([[1, 4, 7],
                    [2, 5, 8],
                    [3, 6, 9]],
                   columns=['x', 'y', 'z'])
    result = df.dictmap({'x': lambda x: x ** 2,
                         'z': lambda z: 1 - z/2})
    tm.assert_frame_equal(result, df)
