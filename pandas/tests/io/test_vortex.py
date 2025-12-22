"""
Tests for Vortex file format support.
"""
import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm


@pytest.mark.single_cpu
class TestVortex:
    """Tests for Vortex I/O operations."""

    def test_basic_roundtrip(self, tmp_path):
        """Test basic write and read roundtrip."""
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({
            "a": [1, 2, 3],
            "b": ["x", "y", "z"],
            "c": [1.1, 2.2, 3.3],
        })
        path = tmp_path / "test.vortex"

        df.to_vortex(path)
        result = read_vortex(path)

        tm.assert_frame_equal(df, result)

    def test_read_column_subset(self, tmp_path):
        """Test reading only a subset of columns."""
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({
            "a": [1, 2, 3],
            "b": [4, 5, 6],
            "c": [7, 8, 9],
        })
        path = tmp_path / "test_cols.vortex"

        df.to_vortex(path)
        result = read_vortex(path, columns=["a", "c"])

        expected = df[["a", "c"]]
        tm.assert_frame_equal(result, expected)

    def test_empty_dataframe(self, tmp_path):
        """Test reading and writing an empty DataFrame."""
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({"a": [], "b": []})
        path = tmp_path / "test_empty.vortex"

        df.to_vortex(path)
        result = read_vortex(path)

        tm.assert_frame_equal(df, result)

    def test_various_dtypes(self, tmp_path):
        """Test DataFrame with various data types."""
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({
            "int": np.array([1, 2, 3], dtype="int64"),
            "float": np.array([1.0, 2.0, 3.0], dtype="float64"),
            "str": ["a", "b", "c"],
            "bool": [True, False, True],
        })
        path = tmp_path / "test_dtypes.vortex"

        df.to_vortex(path)
        result = read_vortex(path)

        tm.assert_frame_equal(df, result)

    def test_with_index(self, tmp_path):
        """Test DataFrame with custom index is preserved as column."""
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame(
            {"a": [1, 2, 3], "b": [4, 5, 6]},
            index=["x", "y", "z"],
        )
        path = tmp_path / "test_index.vortex"

        df.to_vortex(path)
        result = read_vortex(path)

        # Vortex saves index as a column, verify it exists and has correct values
        assert "__index_level_0__" in result.columns
        assert list(result["__index_level_0__"]) == list(df.index)

        # Verify data columns match (ignoring the index type difference)
        tm.assert_frame_equal(
            df.reset_index(drop=True),
            result[["a", "b"]].reset_index(drop=True),
        )
