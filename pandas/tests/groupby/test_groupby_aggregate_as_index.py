import pytest

import pandas as pd
import pandas.testing as pdt


@pytest.fixture
def sample_df():
    return pd.DataFrame(
        {
            "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
            "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
            "C": [1, 2, 3, 4, 5, 6, 7, 8],
            "D": [2.0, 5.0, 8.0, 1.0, 2.0, 9.0, 7.0, 8.0],
        }
    )


def test_groupby_aggregate_as_index_false_with_column_key(sample_df):
    grouped = sample_df.groupby("A", as_index=False)
    result = grouped.aggregate({"C": "sum"})
    expected = pd.DataFrame({"A": ["bar", "foo"], "C": [12, 24]})
    pdt.assert_frame_equal(result, expected)


def test_groupby_aggregate_as_index_false_with_no_grouping_keys(sample_df):
    grouped = sample_df.groupby("A", as_index=False)
    result = grouped.aggregate({"D": "sum"})
    expected = pd.DataFrame({"A": ["bar", "foo"], "D": [15.0, 27.0]})
    pdt.assert_frame_equal(result, expected)
    assert result.index.equals(pd.RangeIndex(start=0, stop=len(result), step=1))
