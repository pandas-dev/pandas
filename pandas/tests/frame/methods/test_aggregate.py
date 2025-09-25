import pandas as pd
import pandas._testing as tm


def test_aggregate_empty_dataframe_returns_series():
    df = pd.DataFrame({"A": [], "B": []})
    result_row = df.aggregate(".".join, axis="columns")
    expected_row = pd.Series([None] * len(df.index), index=df.index, dtype=object)
    tm.assert_series_equal(result_row, expected_row)
    result_col = df.aggregate(".".join, axis="index")
    expected_col = pd.Series([None] * len(df.columns), index=df.columns, dtype=object)
    tm.assert_series_equal(result_col, expected_col)
