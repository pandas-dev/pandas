def test_aggregate_empty_dataframe_returns_series():
    import pandas as pd
    df = pd.DataFrame({"A": [], "B": []})
    # Row-wise aggregation (axis=1)
    result_row = df.aggregate(".".join, axis="columns")
    expected_row = pd.Series([], index=df.index, dtype=object)
    pd.testing.assert_series_equal(result_row, expected_row)
    # Column-wise aggregation (axis=0)
    result_col = df.aggregate(".".join, axis="index")
    expected_col = pd.Series([], dtype=object)
    pd.testing.assert_series_equal(result_col, expected_col)
