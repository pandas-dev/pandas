from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


def test_iterrows_deprecated():
    # GH#43874
    df = pd.DataFrame({"a": [1, 2]})
    with tm.assert_produces_warning(Pandas4Warning, match="iterrows is deprecated"):
        next(df.iterrows())


def test_no_overflow_of_freq_and_time_in_dataframe():
    # GH 35665
    df = pd.DataFrame(
        {
            "some_string": ["2222Y3"],
            "time": [pd.Timedelta("0 days 00:00:00.990000")],
        }
    )
    with tm.assert_produces_warning(Pandas4Warning, match="iterrows is deprecated"):
        for _, row in df.iterrows():
            assert row.dtype == "object"
