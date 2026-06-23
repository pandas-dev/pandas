"""
Tests for GH#61191: DataFrame.duplicated() loses index when DataFrame has no columns.

Covers:
  - DataFrame with rows but zero columns (Bug #1)
  - DataFrame with columns but empty subset= (Bug #2)
  - All three keep= modes for each case
  - Edge cases: empty rows, single row, MultiIndex, non-default index types
  - Regression: existing normal behaviour must be unchanged

Run with:
    pytest tests/frame/methods/test_duplicated_no_columns.py -v
"""

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, Series
import pandas._testing as tm


class TestDuplicatedNoColumns:
    """Bug #1 — DataFrame with rows but no columns."""

    def _make_no_col(self, index=None):
        """Helper: DataFrame with the given index and zero columns."""
        if index is None:
            index = [10, 20, 30]
        return DataFrame(index=index)

    # ── keep='first' (default) ────────────────────────────────────────────

    def test_no_columns_keep_first_int_index(self):
        df = self._make_no_col([10, 20, 30])
        result = df.duplicated()

        expected = Series([False, True, True], index=Index([10, 20, 30]))
        tm.assert_series_equal(result, expected)

    def test_no_columns_keep_first_string_index(self):
        df = self._make_no_col(["alpha", "beta", "gamma"])
        result = df.duplicated()

        expected = Series([False, True, True], index=Index(["alpha", "beta", "gamma"]))
        tm.assert_series_equal(result, expected)

    def test_no_columns_keep_first_float_index(self):
        df = self._make_no_col([0.1, 0.2, 0.3])
        result = df.duplicated()

        expected = Series([False, True, True], index=Index([0.1, 0.2, 0.3]))
        tm.assert_series_equal(result, expected)

    # ── keep='last' ───────────────────────────────────────────────────────

    def test_no_columns_keep_last(self):
        df = self._make_no_col([10, 20, 30])
        result = df.duplicated(keep="last")

        expected = Series([True, True, False], index=Index([10, 20, 30]))
        tm.assert_series_equal(result, expected)

    # ── keep=False ────────────────────────────────────────────────────────

    def test_no_columns_keep_false(self):
        df = self._make_no_col([10, 20, 30])
        result = df.duplicated(keep=False)

        expected = Series([True, True, True], index=Index([10, 20, 30]))
        tm.assert_series_equal(result, expected)

    # ── index is preserved (explicit check) ──────────────────────────────

    def test_no_columns_index_is_preserved(self):
        """Core regression: result.index must match DataFrame.index exactly."""
        custom_index = [100, 200, 300, 400, 500]
        df = self._make_no_col(custom_index)
        result = df.duplicated()

        assert result.index.tolist() == custom_index, (
            "duplicated() dropped/reset the index when DataFrame has no columns"
        )

    # ── result has the correct dtype ──────────────────────────────────────

    def test_no_columns_result_dtype_is_bool(self):
        df = self._make_no_col()
        result = df.duplicated()
        assert result.dtype == np.dtype("bool"), f"Expected bool dtype, got {result.dtype}"

    # ── result length matches DataFrame ──────────────────────────────────

    def test_no_columns_result_length(self):
        for n in (0, 1, 2, 5, 100):
            df = DataFrame(index=range(n))
            result = df.duplicated()
            assert len(result) == n, f"Expected len {n}, got {len(result)}"

    # ── single-row edge case ──────────────────────────────────────────────

    @pytest.mark.parametrize("keep", ["first", "last", False])
    def test_no_columns_single_row_not_duplicate(self, keep):
        """A single row can never be a duplicate of itself."""
        df = DataFrame(index=[99])
        result = df.duplicated(keep=keep)

        expected = Series([False], index=Index([99]))
        tm.assert_series_equal(result, expected)

    # ── empty DataFrame (no rows AND no columns) ──────────────────────────

    @pytest.mark.parametrize("keep", ["first", "last", False])
    def test_fully_empty_dataframe(self, keep):
        """DataFrame with neither rows nor columns should return empty bool Series."""
        df = DataFrame()
        result = df.duplicated(keep=keep)

        expected = Series([], dtype=bool)
        tm.assert_series_equal(result, expected)

    # ── MultiIndex ────────────────────────────────────────────────────────

    def test_no_columns_multiindex(self):
        midx = pd.MultiIndex.from_tuples([(1, "a"), (2, "b"), (3, "c")], names=["x", "y"])
        df = DataFrame(index=midx)
        result = df.duplicated()

        expected = Series([False, True, True], index=midx)
        tm.assert_series_equal(result, expected)

    # ── RangeIndex (default) ─────────────────────────────────────────────

    def test_no_columns_range_index(self):
        df = DataFrame(index=range(4))
        result = df.duplicated()

        expected = Series([False, True, True, True], index=pd.RangeIndex(4))
        tm.assert_series_equal(result, expected)


class TestDuplicatedEmptySubset:
    """Bug #2 — DataFrame.duplicated(subset=[]) with columns present."""

    def _make(self, n=3, index=None):
        if index is None:
            index = range(n)
        return DataFrame({"a": [1] * n, "b": list(range(n))}, index=index)

    # ── keep='first' (default) ────────────────────────────────────────────

    def test_empty_subset_keep_first(self):
        df = self._make(index=[100, 200, 300])
        result = df.duplicated(subset=[])

        expected = Series([False, True, True], index=Index([100, 200, 300]))
        tm.assert_series_equal(result, expected)

    # ── keep='last' ───────────────────────────────────────────────────────

    def test_empty_subset_keep_last(self):
        df = self._make(index=[100, 200, 300])
        result = df.duplicated(subset=[], keep="last")

        expected = Series([True, True, False], index=Index([100, 200, 300]))
        tm.assert_series_equal(result, expected)

    # ── keep=False ────────────────────────────────────────────────────────

    def test_empty_subset_keep_false(self):
        df = self._make(index=[100, 200, 300])
        result = df.duplicated(subset=[], keep=False)

        expected = Series([True, True, True], index=Index([100, 200, 300]))
        tm.assert_series_equal(result, expected)

    # ── index is preserved ────────────────────────────────────────────────

    def test_empty_subset_index_preserved(self):
        idx = [10, 20, 30, 40]
        df = self._make(n=4, index=idx)
        result = df.duplicated(subset=[])

        assert result.index.tolist() == idx

    # ── no ValueError is raised ───────────────────────────────────────────

    def test_empty_subset_no_value_error(self):
        """Regression: previously raised ValueError: not enough values to unpack."""
        df = DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
        try:
            df.duplicated(subset=[])
        except ValueError as exc:
            pytest.fail(f"duplicated(subset=[]) raised unexpected ValueError: {exc}")

    # ── single-row with empty subset ─────────────────────────────────────

    @pytest.mark.parametrize("keep", ["first", "last", False])
    def test_empty_subset_single_row(self, keep):
        df = DataFrame({"a": [42]}, index=[999])
        result = df.duplicated(subset=[], keep=keep)

        expected = Series([False], index=Index([999]))
        tm.assert_series_equal(result, expected)

    # ── explicit empty list vs no subset ─────────────────────────────────

    def test_empty_subset_list_type_accepted(self):
        """Both [] and tuple() should work as empty subset."""
        df = DataFrame({"a": [1, 1, 2]}, index=[1, 2, 3])
        r_list  = df.duplicated(subset=[])
        r_tuple = df.duplicated(subset=())
        tm.assert_series_equal(r_list, r_tuple)


class TestDuplicatedRegressions:
    """Ensure existing behaviour is entirely unchanged."""

    def test_normal_keep_first(self):
        df = pd.DataFrame(
            {"brand": ["A", "A", "B", "B", "B"], "score": [4, 4, 3, 15, 5]}
        )
        result = df.duplicated()
        expected = Series([False, True, False, False, False])
        tm.assert_series_equal(result, expected)

    def test_normal_keep_last(self):
        df = pd.DataFrame({"x": [1, 1, 2]})
        result = df.duplicated(keep="last")
        expected = Series([True, False, False])
        tm.assert_series_equal(result, expected)

    def test_normal_keep_false(self):
        df = pd.DataFrame({"x": [1, 1, 2]})
        result = df.duplicated(keep=False)
        expected = Series([True, True, False])
        tm.assert_series_equal(result, expected)

    def test_subset_selection(self):
        df = pd.DataFrame({"a": [1, 1, 2], "b": [3, 4, 5]})
        result = df.duplicated(subset=["a"])
        expected = Series([False, True, False])
        tm.assert_series_equal(result, expected)

    def test_invalid_subset_raises_keyerror(self):
        df = pd.DataFrame({"a": [1, 2]})
        with pytest.raises(KeyError):
            df.duplicated(subset=["nonexistent_column"])

    def test_empty_rows_dataframe(self):
        df = pd.DataFrame({"a": [], "b": []})
        result = df.duplicated()
        expected = Series([], dtype=bool)
        tm.assert_series_equal(result, expected)
