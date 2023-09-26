import numpy as np

from pandas import (
    DataFrame,
    case_when,
)
import pandas._testing as tm


# use fixture and parametrize
def test_case_when_multiple_conditions():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    result = df.assign(new_column=case_when(df.a.eq(1), 1, df.a.gt(1) & df.b.eq(5), 2))
    expected = df.assign(new_column=[1, 2, np.nan])
    tm.assert_frame_equal(result, expected)
