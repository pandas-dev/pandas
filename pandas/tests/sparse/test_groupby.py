import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
@pytest.mark.filterwarnings("ignore:DataFrame.to_sparse:FutureWarning")
class TestSparseGroupBy:
    def setup_method(self, method):
        self.dense = pd.DataFrame(
            {
                "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
                "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
                "C": np.random.randn(8),
                "D": np.random.randn(8),
                "E": [np.nan, np.nan, 1, 2, np.nan, 1, np.nan, np.nan],
            }
        )
        self.sparse = self.dense.to_sparse()

    def test_first_last_nth(self):
        # tests for first / last / nth
        sparse_grouped = self.sparse.groupby("A")
        dense_grouped = self.dense.groupby("A")

        sparse_grouped_first = sparse_grouped.first()
        sparse_grouped_last = sparse_grouped.last()
        sparse_grouped_nth = sparse_grouped.nth(1)

        dense_grouped_first = pd.DataFrame(dense_grouped.first().to_sparse())
        dense_grouped_last = pd.DataFrame(dense_grouped.last().to_sparse())
        dense_grouped_nth = pd.DataFrame(dense_grouped.nth(1).to_sparse())

        tm.assert_frame_equal(sparse_grouped_first, dense_grouped_first)
        tm.assert_frame_equal(sparse_grouped_last, dense_grouped_last)
        tm.assert_frame_equal(sparse_grouped_nth, dense_grouped_nth)

    def test_aggfuncs(self):
        sparse_grouped = self.sparse.groupby("A")
        dense_grouped = self.dense.groupby("A")

        result = sparse_grouped.mean().to_sparse()
        expected = dense_grouped.mean().to_sparse()

        tm.assert_frame_equal(result, expected)

        # ToDo: sparse sum includes str column
        # tm.assert_frame_equal(sparse_grouped.sum(),
        #                       dense_grouped.sum())

        result = sparse_grouped.count().to_sparse()
        expected = dense_grouped.count().to_sparse()

        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("fill_value", [0, np.nan])
@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
@pytest.mark.filterwarnings("ignore:DataFrame.to_sparse:FutureWarning")
def test_groupby_includes_fill_value(fill_value):
    # https://github.com/pandas-dev/pandas/issues/5078
    df = pd.DataFrame(
        {
            "a": [fill_value, 1, fill_value, fill_value],
            "b": [fill_value, 1, fill_value, fill_value],
        }
    )
    sdf = df.to_sparse(fill_value=fill_value)
    result = sdf.groupby("a").sum()
    expected = pd.DataFrame(df.groupby("a").sum().to_sparse(fill_value=fill_value))
    tm.assert_frame_equal(result, expected, check_index_type=False)
