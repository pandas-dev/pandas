import numpy as np

import pandas as pd
import pandas._testing as tm


def test_unsure_about_the_name_and_location():
    # https://github.com/pandas-dev/pandas/issues/33198
    arr = np.random.randn(10 ** 6).reshape(500, 2000).astype(np.float64)
    df = pd.DataFrame(arr)
    df.iloc[:, 1000:] = df.iloc[:, 1000:].astype(np.float32)
    expected = pd.Series(
        {np.dtype("float32"): 1000, np.dtype("float64"): 1000}, dtype="int64"
    )
    result = df.dtypes.value_counts()
    tm.assert_series_equal(result, expected)
