import numpy as np

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    date_range,
    timedelta_range,
)
import pandas._testing as tm


def _check_mixed_int(df, dtype=None):
    # GH#41672
    result = DataFrame([], columns=['lang', 'name'])
    result = result.agg({'name': lambda y: y.values})
    assert type(result) == Series

    result = DataFrame([['a', 'boof']], columns=['lang', 'name'])
    result = result.agg({'name': lambda y: y.values})
    assert type(result) == Series
