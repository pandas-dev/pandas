import numpy as np
import pytest

import pandas as pd
from pandas import _np_version_under1p17
import pandas.util.testing as tm


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
@pytest.mark.filterwarnings("ignore:Series.to_sparse:FutureWarning")
@pytest.mark.filterwarnings("ignore:DataFrame.to_sparse:FutureWarning")
class TestPivotTable:
    def setup_method(self, method):
        rs = np.random.RandomState(0)
        self.dense = pd.DataFrame(
            {
                "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
                "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
                "C": rs.randn(8),
                "D": rs.randn(8),
                "E": [np.nan, np.nan, 1, 2, np.nan, 1, np.nan, np.nan],
            }
        )
        self.sparse = self.dense.to_sparse()

    def test_pivot_table(self):
        res_sparse = pd.pivot_table(self.sparse, index="A", columns="B", values="C")
        res_dense = pd.pivot_table(self.dense, index="A", columns="B", values="C")
        tm.assert_frame_equal(res_sparse, res_dense)

        res_sparse = pd.pivot_table(self.sparse, index="A", columns="B", values="E")
        res_dense = pd.pivot_table(self.dense, index="A", columns="B", values="E")
        tm.assert_frame_equal(res_sparse, res_dense)

        res_sparse = pd.pivot_table(
            self.sparse, index="A", columns="B", values="E", aggfunc="mean"
        )
        res_dense = pd.pivot_table(
            self.dense, index="A", columns="B", values="E", aggfunc="mean"
        )
        tm.assert_frame_equal(res_sparse, res_dense)

    def test_pivot_table_with_nans(self):
        res_sparse = pd.pivot_table(
            self.sparse, index="A", columns="B", values="E", aggfunc="sum"
        )
        res_dense = pd.pivot_table(
            self.dense, index="A", columns="B", values="E", aggfunc="sum"
        )
        tm.assert_frame_equal(res_sparse, res_dense)

    @pytest.mark.xfail(
        not _np_version_under1p17,
        reason="failing occasionally on numpy > 1.17",
        strict=False,
    )
    def test_pivot_table_multi(self):
        res_sparse = pd.pivot_table(
            self.sparse, index="A", columns="B", values=["D", "E"]
        )
        res_dense = pd.pivot_table(
            self.dense, index="A", columns="B", values=["D", "E"]
        )
        res_dense = res_dense.apply(lambda x: x.astype("Sparse[float64]"))
        tm.assert_frame_equal(res_sparse, res_dense)
