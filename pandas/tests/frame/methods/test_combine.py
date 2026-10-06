import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestCombine:
    @pytest.mark.parametrize(
        "data",
        [
            pd.date_range("2000", periods=4),
            pd.date_range("2000", periods=4, tz="US/Central"),
            pd.period_range("2000", periods=4),
            pd.timedelta_range(0, periods=4),
        ],
    )
    def test_combine_datetlike_udf(self, data):
        # GH#23079
        df = pd.DataFrame({"A": data})
        other = df.copy()
        df.iloc[1, 0] = None

        def combiner(a, b):
            return b

        result = df.combine(other, combiner)
        tm.assert_frame_equal(result, other)

    def test_combine_generic(self, float_frame):
        df1 = float_frame
        df2 = float_frame.loc[float_frame.index[:-5], ["A", "B", "C"]]

        combined = df1.combine(df2, np.add)
        combined2 = df2.combine(df1, np.add)
        assert combined["D"].isna().all()
        assert combined2["D"].isna().all()

        chunk = combined.loc[combined.index[:-5], ["A", "B", "C"]]
        chunk2 = combined2.loc[combined2.index[:-5], ["A", "B", "C"]]

        exp = (
            float_frame.loc[float_frame.index[:-5], ["A", "B", "C"]].reindex_like(chunk)
            * 2
        )
        tm.assert_frame_equal(chunk, exp)
        tm.assert_frame_equal(chunk2, exp)

    def test_combine_nonunique_columns(self):
        # GH#51340

        df = pd.DataFrame({"A": range(5), "B": range(5)})
        df.columns = ["A", "A"]

        other = df.copy()
        df.iloc[1, :] = None

        def combiner(a, b):
            return b

        result = df.combine(other, combiner)
        expected = other.astype("float64")
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("value", [2**63, 2**64, -(2**64), 2.0**64, np.inf])
def test_combine_object_result_not_castable_to_int(value):
    # GH#66394 a func result that doesn't fit the common dtype must be left
    #  alone rather than raising OverflowError out of maybe_downcast_numeric
    df1 = pd.DataFrame({"a": [1, 2]})
    df2 = pd.DataFrame({"a": [3, 4]})

    result = df1.combine(
        df2, lambda ser1, ser2: pd.Series([value, value], dtype=object)
    )

    expected = pd.DataFrame({"a": pd.Series([value, value], dtype=object)})
    tm.assert_frame_equal(result, expected)


def test_combine_object_result_castable_to_int():
    # GH#66394 in-range results still downcast to the common dtype
    df1 = pd.DataFrame({"a": [1, 2]})
    df2 = pd.DataFrame({"a": [3, 4]})

    result = df1.combine(df2, lambda ser1, ser2: pd.Series([5, 6], dtype=object))

    expected = pd.DataFrame({"a": [5, 6]})
    tm.assert_frame_equal(result, expected)
