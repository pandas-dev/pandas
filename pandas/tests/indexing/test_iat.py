import pandas as pd


def test_iat(float_frame):

    for i, row in enumerate(float_frame.index):
        for j, col in enumerate(float_frame.columns):
            result = float_frame.iat[i, j]
            expected = float_frame.at[row, col]
            assert result == expected


def test_iat_duplicate_columns():
    # https://github.com/pandas-dev/pandas/issues/11754
    df = pd.DataFrame([[1, 2]], columns=["x", "x"])
    assert df.iat[0, 0] == 1
