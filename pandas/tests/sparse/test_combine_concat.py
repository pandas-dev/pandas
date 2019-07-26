import itertools

import numpy as np
import pytest

from pandas.errors import PerformanceWarning

import pandas as pd
import pandas.util.testing as tm


class TestSparseArrayConcat:
    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_basic(self, kind):
        a = pd.SparseArray([1, 0, 0, 2], kind=kind)
        b = pd.SparseArray([1, 0, 2, 2], kind=kind)

        result = pd.SparseArray._concat_same_type([a, b])
        # Can't make any assertions about the sparse index itself
        # since we aren't don't merge sparse blocs across arrays
        # in to_concat
        expected = np.array([1, 2, 1, 2, 2], dtype="int64")
        tm.assert_numpy_array_equal(result.sp_values, expected)
        assert result.kind == kind

    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_uses_first_kind(self, kind):
        other = "integer" if kind == "block" else "block"
        a = pd.SparseArray([1, 0, 0, 2], kind=kind)
        b = pd.SparseArray([1, 0, 2, 2], kind=other)

        result = pd.SparseArray._concat_same_type([a, b])
        expected = np.array([1, 2, 1, 2, 2], dtype="int64")
        tm.assert_numpy_array_equal(result.sp_values, expected)
        assert result.kind == kind


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
class TestSparseSeriesConcat:
    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_concat(self, kind):
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        sparse1 = pd.SparseSeries(val1, name="x", kind=kind)
        sparse2 = pd.SparseSeries(val2, name="y", kind=kind)

        res = pd.concat([sparse1, sparse2])
        exp = pd.concat([pd.Series(val1), pd.Series(val2)])
        exp = pd.SparseSeries(exp, kind=kind)
        tm.assert_sp_series_equal(res, exp, consolidate_block_indices=True)

        sparse1 = pd.SparseSeries(val1, fill_value=0, name="x", kind=kind)
        sparse2 = pd.SparseSeries(val2, fill_value=0, name="y", kind=kind)

        res = pd.concat([sparse1, sparse2])
        exp = pd.concat([pd.Series(val1), pd.Series(val2)])
        exp = pd.SparseSeries(exp, fill_value=0, kind=kind)
        tm.assert_sp_series_equal(res, exp, consolidate_block_indices=True)

    def test_concat_axis1(self):
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        sparse1 = pd.SparseSeries(val1, name="x")
        sparse2 = pd.SparseSeries(val2, name="y")

        res = pd.concat([sparse1, sparse2], axis=1)
        exp = pd.concat([pd.Series(val1, name="x"), pd.Series(val2, name="y")], axis=1)
        exp = pd.SparseDataFrame(exp)
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

    def test_concat_different_fill(self):
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        for kind in ["integer", "block"]:
            sparse1 = pd.SparseSeries(val1, name="x", kind=kind)
            sparse2 = pd.SparseSeries(val2, name="y", kind=kind, fill_value=0)

            with tm.assert_produces_warning(
                PerformanceWarning, raise_on_extra_warnings=False
            ):
                res = pd.concat([sparse1, sparse2])

            exp = pd.concat([pd.Series(val1), pd.Series(val2)])
            exp = pd.SparseSeries(exp, kind=kind)
            tm.assert_sp_series_equal(res, exp)

            with tm.assert_produces_warning(
                PerformanceWarning, raise_on_extra_warnings=False
            ):
                res = pd.concat([sparse2, sparse1])

            exp = pd.concat([pd.Series(val2), pd.Series(val1)])
            exp = pd.SparseSeries(exp, kind=kind, fill_value=0)
            tm.assert_sp_series_equal(res, exp)

    def test_concat_axis1_different_fill(self):
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        sparse1 = pd.SparseSeries(val1, name="x")
        sparse2 = pd.SparseSeries(val2, name="y", fill_value=0)

        res = pd.concat([sparse1, sparse2], axis=1)
        exp = pd.concat([pd.Series(val1, name="x"), pd.Series(val2, name="y")], axis=1)
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), exp)

    def test_concat_different_kind(self):
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        sparse1 = pd.SparseSeries(val1, name="x", kind="integer")
        sparse2 = pd.SparseSeries(val2, name="y", kind="block")

        res = pd.concat([sparse1, sparse2])
        exp = pd.concat([pd.Series(val1), pd.Series(val2)])
        exp = pd.SparseSeries(exp, kind=sparse1.kind)
        tm.assert_sp_series_equal(res, exp)

        res = pd.concat([sparse2, sparse1])
        exp = pd.concat([pd.Series(val2), pd.Series(val1)])
        exp = pd.SparseSeries(exp, kind=sparse2.kind)
        tm.assert_sp_series_equal(res, exp, consolidate_block_indices=True)

    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_concat_sparse_dense(self, kind):
        # use first input's fill_value
        val1 = np.array([1, 2, np.nan, np.nan, 0, np.nan])
        val2 = np.array([3, np.nan, 4, 0, 0])

        sparse = pd.SparseSeries(val1, name="x", kind=kind)
        dense = pd.Series(val2, name="y")

        res = pd.concat([sparse, dense])
        exp = pd.SparseSeries(pd.concat([pd.Series(val1), dense]), kind=kind)
        tm.assert_sp_series_equal(res, exp)

        res = pd.concat([dense, sparse, dense])
        exp = pd.concat([dense, pd.Series(val1), dense])
        # XXX: changed from SparseSeries to Series[sparse]
        exp = pd.Series(pd.SparseArray(exp, kind=kind), index=exp.index, name=exp.name)
        tm.assert_series_equal(res, exp)

        sparse = pd.SparseSeries(val1, name="x", kind=kind, fill_value=0)
        dense = pd.Series(val2, name="y")

        res = pd.concat([sparse, dense])
        # XXX: changed from SparseSeries to Series[sparse]
        exp = pd.concat([pd.Series(val1), dense])
        exp = pd.Series(
            pd.SparseArray(exp, kind=kind, fill_value=0), index=exp.index, name=exp.name
        )
        tm.assert_series_equal(res, exp)

        res = pd.concat([dense, sparse, dense])
        exp = pd.concat([dense, pd.Series(val1), dense])
        # XXX: changed from SparseSeries to Series[sparse]
        exp = pd.Series(
            pd.SparseArray(exp, kind=kind, fill_value=0), index=exp.index, name=exp.name
        )
        tm.assert_series_equal(res, exp)


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
@pytest.mark.filterwarnings("ignore:DataFrame.to_sparse:FutureWarning")
class TestSparseDataFrameConcat:
    def setup_method(self, method):

        self.dense1 = pd.DataFrame(
            {
                "A": [0.0, 1.0, 2.0, np.nan],
                "B": [0.0, 0.0, 0.0, 0.0],
                "C": [np.nan, np.nan, np.nan, np.nan],
                "D": [1.0, 2.0, 3.0, 4.0],
            }
        )

        self.dense2 = pd.DataFrame(
            {
                "A": [5.0, 6.0, 7.0, 8.0],
                "B": [np.nan, 0.0, 7.0, 8.0],
                "C": [5.0, 6.0, np.nan, np.nan],
                "D": [np.nan, np.nan, np.nan, np.nan],
            }
        )

        self.dense3 = pd.DataFrame(
            {
                "E": [5.0, 6.0, 7.0, 8.0],
                "F": [np.nan, 0.0, 7.0, 8.0],
                "G": [5.0, 6.0, np.nan, np.nan],
                "H": [np.nan, np.nan, np.nan, np.nan],
            }
        )

    def test_concat(self):
        # fill_value = np.nan
        sparse = self.dense1.to_sparse()
        sparse2 = self.dense2.to_sparse()

        res = pd.concat([sparse, sparse])
        exp = pd.concat([self.dense1, self.dense1]).to_sparse()
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse2, sparse2])
        exp = pd.concat([self.dense2, self.dense2]).to_sparse()
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse, sparse2])
        exp = pd.concat([self.dense1, self.dense2]).to_sparse()
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse2, sparse])
        exp = pd.concat([self.dense2, self.dense1]).to_sparse()
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        # fill_value = 0
        sparse = self.dense1.to_sparse(fill_value=0)
        sparse2 = self.dense2.to_sparse(fill_value=0)

        res = pd.concat([sparse, sparse])
        exp = pd.concat([self.dense1, self.dense1]).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse2, sparse2])
        exp = pd.concat([self.dense2, self.dense2]).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse, sparse2])
        exp = pd.concat([self.dense1, self.dense2]).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        res = pd.concat([sparse2, sparse])
        exp = pd.concat([self.dense2, self.dense1]).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

    def test_concat_different_fill_value(self):
        # 1st fill_value will be used
        sparse = self.dense1.to_sparse()
        sparse2 = self.dense2.to_sparse(fill_value=0)

        with tm.assert_produces_warning(
            PerformanceWarning, raise_on_extra_warnings=False
        ):
            res = pd.concat([sparse, sparse2])
        exp = pd.concat([self.dense1, self.dense2]).to_sparse()
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

        with tm.assert_produces_warning(
            PerformanceWarning, raise_on_extra_warnings=False
        ):
            res = pd.concat([sparse2, sparse])
        exp = pd.concat([self.dense2, self.dense1]).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, consolidate_block_indices=True)

    def test_concat_different_columns_sort_warns(self):
        sparse = self.dense1.to_sparse()
        sparse3 = self.dense3.to_sparse()

        # stacklevel is wrong since we have two FutureWarnings,
        # one for depr, one for sorting.
        with tm.assert_produces_warning(
            FutureWarning, check_stacklevel=False, raise_on_extra_warnings=False
        ):
            res = pd.concat([sparse, sparse3])
        with tm.assert_produces_warning(
            FutureWarning, check_stacklevel=False, raise_on_extra_warnings=False
        ):
            exp = pd.concat([self.dense1, self.dense3])

        exp = exp.to_sparse()
        tm.assert_sp_frame_equal(res, exp, check_kind=False)

    def test_concat_different_columns(self):
        # fill_value = np.nan
        sparse = self.dense1.to_sparse()
        sparse3 = self.dense3.to_sparse()

        res = pd.concat([sparse, sparse3], sort=True)
        exp = pd.concat([self.dense1, self.dense3], sort=True).to_sparse()
        tm.assert_sp_frame_equal(res, exp, check_kind=False)

        res = pd.concat([sparse3, sparse], sort=True)
        exp = pd.concat([self.dense3, self.dense1], sort=True).to_sparse()
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp, check_kind=False)

    def test_concat_bug(self):
        from pandas.core.sparse.api import SparseDtype

        x = pd.SparseDataFrame({"A": pd.SparseArray([np.nan, np.nan], fill_value=0)})
        y = pd.SparseDataFrame({"B": []})
        res = pd.concat([x, y], sort=False)[["A"]]
        exp = pd.DataFrame(
            {"A": pd.SparseArray([np.nan, np.nan], dtype=SparseDtype(float, 0))}
        )
        tm.assert_frame_equal(res, exp)

    def test_concat_different_columns_buggy(self):
        sparse = self.dense1.to_sparse(fill_value=0)
        sparse3 = self.dense3.to_sparse(fill_value=0)

        res = pd.concat([sparse, sparse3], sort=True)
        exp = pd.concat([self.dense1, self.dense3], sort=True).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan

        tm.assert_sp_frame_equal(
            res, exp, check_kind=False, consolidate_block_indices=True
        )

        res = pd.concat([sparse3, sparse], sort=True)
        exp = pd.concat([self.dense3, self.dense1], sort=True).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(
            res, exp, check_kind=False, consolidate_block_indices=True
        )

        # different fill values
        sparse = self.dense1.to_sparse()
        sparse3 = self.dense3.to_sparse(fill_value=0)
        # each columns keeps its fill_value, thus compare in dense
        res = pd.concat([sparse, sparse3], sort=True)
        exp = pd.concat([self.dense1, self.dense3], sort=True)
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), exp)

        res = pd.concat([sparse3, sparse], sort=True)
        exp = pd.concat([self.dense3, self.dense1], sort=True)
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), exp)

    def test_concat_series(self):
        # fill_value = np.nan
        sparse = self.dense1.to_sparse()
        sparse2 = self.dense2.to_sparse()

        for col in ["A", "D"]:
            res = pd.concat([sparse, sparse2[col]])
            exp = pd.concat([self.dense1, self.dense2[col]]).to_sparse()
            tm.assert_sp_frame_equal(res, exp, check_kind=False)

            res = pd.concat([sparse2[col], sparse])
            exp = pd.concat([self.dense2[col], self.dense1]).to_sparse()
            tm.assert_sp_frame_equal(res, exp, check_kind=False)

        # fill_value = 0
        sparse = self.dense1.to_sparse(fill_value=0)
        sparse2 = self.dense2.to_sparse(fill_value=0)

        for col in ["C", "D"]:
            res = pd.concat([sparse, sparse2[col]])
            exp = pd.concat([self.dense1, self.dense2[col]]).to_sparse(fill_value=0)
            exp._default_fill_value = np.nan
            tm.assert_sp_frame_equal(
                res, exp, check_kind=False, consolidate_block_indices=True
            )

            res = pd.concat([sparse2[col], sparse])
            exp = pd.concat([self.dense2[col], self.dense1]).to_sparse(fill_value=0)
            exp["C"] = res["C"]
            exp._default_fill_value = np.nan
            tm.assert_sp_frame_equal(
                res, exp, consolidate_block_indices=True, check_kind=False
            )

    def test_concat_axis1(self):
        # fill_value = np.nan
        sparse = self.dense1.to_sparse()
        sparse3 = self.dense3.to_sparse()

        res = pd.concat([sparse, sparse3], axis=1)
        exp = pd.concat([self.dense1, self.dense3], axis=1).to_sparse()
        tm.assert_sp_frame_equal(res, exp)

        res = pd.concat([sparse3, sparse], axis=1)
        exp = pd.concat([self.dense3, self.dense1], axis=1).to_sparse()
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # fill_value = 0
        sparse = self.dense1.to_sparse(fill_value=0)
        sparse3 = self.dense3.to_sparse(fill_value=0)

        res = pd.concat([sparse, sparse3], axis=1)
        exp = pd.concat([self.dense1, self.dense3], axis=1).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        res = pd.concat([sparse3, sparse], axis=1)
        exp = pd.concat([self.dense3, self.dense1], axis=1).to_sparse(fill_value=0)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # different fill values
        sparse = self.dense1.to_sparse()
        sparse3 = self.dense3.to_sparse(fill_value=0)
        # each columns keeps its fill_value, thus compare in dense
        res = pd.concat([sparse, sparse3], axis=1)
        exp = pd.concat([self.dense1, self.dense3], axis=1)
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), exp)

        res = pd.concat([sparse3, sparse], axis=1)
        exp = pd.concat([self.dense3, self.dense1], axis=1)
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), exp)

    @pytest.mark.parametrize(
        "fill_value,sparse_idx,dense_idx",
        itertools.product([None, 0, 1, np.nan], [0, 1], [1, 0]),
    )
    def test_concat_sparse_dense_rows(self, fill_value, sparse_idx, dense_idx):
        frames = [self.dense1, self.dense2]
        sparse_frame = [
            frames[dense_idx],
            frames[sparse_idx].to_sparse(fill_value=fill_value),
        ]
        dense_frame = [frames[dense_idx], frames[sparse_idx]]

        # This will try both directions sparse + dense and dense + sparse
        for _ in range(2):
            res = pd.concat(sparse_frame)
            exp = pd.concat(dense_frame)

            assert isinstance(res, pd.SparseDataFrame)
            tm.assert_frame_equal(res.to_dense(), exp)

            sparse_frame = sparse_frame[::-1]
            dense_frame = dense_frame[::-1]

    @pytest.mark.parametrize(
        "fill_value,sparse_idx,dense_idx",
        itertools.product([None, 0, 1, np.nan], [0, 1], [1, 0]),
    )
    @pytest.mark.xfail(reason="The iloc fails and I can't make expected", strict=False)
    def test_concat_sparse_dense_cols(self, fill_value, sparse_idx, dense_idx):
        # See GH16874, GH18914 and #18686 for why this should be a DataFrame
        from pandas.core.dtypes.common import is_sparse

        frames = [self.dense1, self.dense3]

        sparse_frame = [
            frames[dense_idx],
            frames[sparse_idx].to_sparse(fill_value=fill_value),
        ]
        dense_frame = [frames[dense_idx], frames[sparse_idx]]

        # This will try both directions sparse + dense and dense + sparse
        for _ in range(2):
            res = pd.concat(sparse_frame, axis=1)
            exp = pd.concat(dense_frame, axis=1)
            cols = [i for (i, x) in enumerate(res.dtypes) if is_sparse(x)]

            for col in cols:
                exp.iloc[:, col] = exp.iloc[:, col].astype("Sparse")

            for column in frames[dense_idx].columns:
                if dense_idx == sparse_idx:
                    tm.assert_frame_equal(res[column], exp[column])
                else:
                    tm.assert_series_equal(res[column], exp[column])

            tm.assert_frame_equal(res, exp)

            sparse_frame = sparse_frame[::-1]
            dense_frame = dense_frame[::-1]
