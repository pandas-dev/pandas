import numpy as np

import pandas as pd
import pandas._testing as tm


def test_rolling_timedelta_inf():
    # GH 57549
    idx = pd.date_range("2024-02-01", "2024-02-21", freq="1d")
    df = pd.DataFrame(index=idx)
    df["A"] = df.index.day.astype(float)
    df.iloc[-2, 0] = np.inf
    df.iloc[-9, 0] = np.nan

    res1 = df.rolling(window=pd.Timedelta("4d"), min_periods=1).mean()
    res2 = df.rolling(window=4, min_periods=1).mean()

    tm.assert_frame_equal(res1, res2)
