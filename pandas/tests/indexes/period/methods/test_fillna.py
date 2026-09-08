from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestFillNA:
    def test_fillna_period(self):
        # GH#11343
        idx = pd.PeriodIndex(["2011-01-01 09:00", pd.NaT, "2011-01-01 11:00"], freq="h")

        exp = pd.PeriodIndex(
            ["2011-01-01 09:00", "2011-01-01 10:00", "2011-01-01 11:00"], freq="h"
        )
        result = idx.fillna(pd.Period("2011-01-01 10:00", freq="h"))
        tm.assert_index_equal(result, exp)

        exp = pd.Index(
            [
                pd.Period("2011-01-01 09:00", freq="h"),
                "x",
                pd.Period("2011-01-01 11:00", freq="h"),
            ],
            dtype=object,
        )
        # GH#45153 filling with incompatible value is deprecated
        with tm.assert_produces_warning(Pandas4Warning, match="fill value"):
            result = idx.fillna("x")
        tm.assert_index_equal(result, exp)

        exp = pd.Index(
            [
                pd.Period("2011-01-01 09:00", freq="h"),
                pd.Period("2011-01-01", freq="D"),
                pd.Period("2011-01-01 11:00", freq="h"),
            ],
            dtype=object,
        )
        # GH#45153 filling with incompatible value is deprecated
        with tm.assert_produces_warning(Pandas4Warning, match="fill value"):
            result = idx.fillna(pd.Period("2011-01-01", freq="D"))
        tm.assert_index_equal(result, exp)
