from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestFillNA:
    def test_fillna_timedelta(self):
        # GH#11343
        idx = pd.TimedeltaIndex(["1 day", pd.NaT, "3 day"])

        exp = pd.TimedeltaIndex(["1 day", "2 day", "3 day"])
        tm.assert_index_equal(idx.fillna(pd.Timedelta("2 day")), exp)

        exp = pd.TimedeltaIndex(["1 day", "3 hour", "3 day"])
        idx.fillna(pd.Timedelta("3 hour"))

        exp = pd.Index(
            [pd.Timedelta("1 day"), "x", pd.Timedelta("3 day")], dtype=object
        )
        # GH#45153 filling with incompatible value is deprecated
        with tm.assert_produces_warning(Pandas4Warning, match="fill value"):
            result = idx.fillna("x")
        tm.assert_index_equal(result, exp)
