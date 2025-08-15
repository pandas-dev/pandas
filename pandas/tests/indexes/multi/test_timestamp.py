import pandas as pd
import pytest

pytest.importorskip("pyarrow")

def test_difference_with_pyarrow_timestamp():
    dates = pd.Series(
        ["2024-01-01", "2024-01-02"], dtype="timestamp[ns][pyarrow]"
    )
    ids = [1, 2]

    mi = pd.MultiIndex.from_arrays([ids, dates], names=["id", "date"])
    to_remove = mi[:1]

    result = mi.difference(to_remove)

    expected_dates = pd.Series(
        ["2024-01-02"], dtype="timestamp[ns][pyarrow]"
    )
    expected_ids = [2]
    expected = pd.MultiIndex.from_arrays(
        [expected_ids, expected_dates], names=["id", "date"]
    )

    pd.testing.assert_index_equal(result, expected)
