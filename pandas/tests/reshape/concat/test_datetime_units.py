import numpy as np
import pandas as pd


def test_concat_series_datetime64_different_units():
    a = pd.Series(np.array(["2020-01-01T00:00:00"], dtype="datetime64[s]"))
    b = pd.Series(np.array(["2020-01-02T00:00:00"], dtype="datetime64[ms]"))

    result = pd.concat([a, b], ignore_index=True)

    # dtype should still be datetime64 (M)
    assert result.dtype.kind == "M"

    assert result.iloc[0] == pd.Timestamp("2020-01-01")
    assert result.iloc[1] == pd.Timestamp("2020-01-02")


def test_concat_dataframe_datetime64_different_units():
    df1 = pd.DataFrame({
        "ts": np.array(["2020-01-01T00:00:00"], dtype="datetime64[s]")
    })
    df2 = pd.DataFrame({
        "ts": np.array(["2020-01-02T00:00:00"], dtype="datetime64[ms]")
    })

    result = pd.concat([df1, df2], ignore_index=True)

    assert result["ts"].dtype.kind == "M"
    assert result["ts"].iloc[0] == pd.Timestamp("2020-01-01")
    assert result["ts"].iloc[1] == pd.Timestamp("2020-01-02")


def test_concat_datetime64_preserves_unit_order():
    # Ensures that units are properly aligned, and no value shifts occur
    a = pd.Series(np.array(["2020-01-01T00:00:00"], dtype="datetime64[us]"))
    b = pd.Series(np.array(["2020-01-01T00:00:01"], dtype="datetime64[ns]"))

    result = pd.concat([a, b], ignore_index=True)

    assert result.dtype.kind == "M"
    assert result.iloc[0] == pd.Timestamp("2020-01-01 00:00:00")
    assert result.iloc[1] == pd.Timestamp("2020-01-01 00:00:01")
