from datetime import datetime

import numpy as np

import pandas as pd
import pandas._testing as tm


def test_drop_duplicates_series_vs_dataframe():
    # GH 14192
    df = pd.DataFrame(
        {
            "a": [1, 1, 1, "one", "one"],
            "b": [2, 2, np.nan, np.nan, np.nan],
            "c": [3, 3, np.nan, np.nan, "three"],
            "d": [1, 2, 3, 4, 4],
            "e": [
                datetime(2015, 1, 1),
                datetime(2015, 1, 1),
                datetime(2015, 2, 1),
                pd.NaT,
                pd.NaT,
            ],
        }
    )
    for column in df.columns:
        for keep in ["first", "last", False]:
            dropped_frame = df[[column]].drop_duplicates(keep=keep)
            dropped_series = df[column].drop_duplicates(keep=keep)
            tm.assert_frame_equal(dropped_frame, dropped_series.to_frame())
