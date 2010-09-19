import numpy as np
import pandas.stats.moments as moments

def test_rolling_median():
    arr = np.random.randn(100)
    arr[20:40] = np.NaN

    result = moments.rolling_median(arr, 50)

    assert(np.isnan(result[20]))

    assert(result[-1] == np.median(arr[-50:]))

    result = moments.rolling_median(arr, 49)

    assert(np.isnan(result[20]))

    assert(result[-1] == np.median(arr[-49:]))
