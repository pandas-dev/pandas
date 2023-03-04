import pytest

from pandas import (
    NA,
    DataFrame,
    Series,
)
from pandas.core.arrays.arrow import ArrowDtype

pa = pytest.importorskip("pyarrow", minversion="7.0.0")
import pandas._testing as tm


@pytest.fixture(scope="module")
def sample_dataframe_numpy_backend():
    return DataFrame(
        {
            "u8": Series([1, 2, 3, NA], dtype="UInt8"),
            "f64": Series([float("NaN"), 1.0, 2.0, 3.0], dtype="float64"),
            "s": Series(["foo", "bar", None, "foobar"], dtype="object"),
        }
    )


@pytest.fixture(scope="module")
def sample_dataframe_pyarrow_backend():
    return DataFrame(
        {
            "u8": Series([1, 2, 3, NA], dtype="uint8[pyarrow]"),
            "f64": Series([NA, 1.0, 2.0, 3.0], dtype="float64[pyarrow]"),
            "s": Series(["foo", "bar", NA, "foobar"], dtype="string[pyarrow]"),
        }
    )


@pytest.fixture(scope="module")
def sample_pyarrow_table():
    return pa.table(
        [
            pa.array([1, 2, 3, None], type=pa.uint8()),
            pa.array([None, 1.0, 2.0, 3.0], type=pa.float64()),
            pa.array(["foo", "bar", None, "foobar"], type=pa.string()),
        ],
        names=["u8", "f64", "s"],
    )


class TestPyArrow:
    @pytest.mark.parametrize(
        "column,dtype", [("u8", pa.uint8()), ("f64", pa.float64()), ("s", pa.string())]
    )
    def test_from_pyarrow_uses_right_pandas_types(
        self, sample_pyarrow_table, column, dtype
    ):
        result = DataFrame.from_pyarrow(sample_pyarrow_table)
        assert result[column].dtype == ArrowDtype(dtype)

    @pytest.mark.parametrize("column", ["u8", "f64", "s"])
    def test_from_pyarrow_keeps_correct_data(self, sample_pyarrow_table, column):
        result = DataFrame.from_pyarrow(sample_pyarrow_table)
        assert result[column]._data.array._data == sample_pyarrow_table[column]

    @pytest.mark.parametrize("column", ["u8", "f64", "s"])
    def test_from_pyarrow_does_not_copy_memory(self, sample_pyarrow_table, column):
        result = DataFrame.from_pyarrow(sample_pyarrow_table)

        result_buffers = result[column]._data.array._data.chunks[0].buffers()
        expected_buffers = sample_pyarrow_table[column].chunks[0].buffers()

        for result_buffer, expected_buffer in zip(result_buffers, expected_buffers):
            if result_buffer is None and expected_buffer is None:
                continue
            assert result_buffer.address == expected_buffer.address
            assert result_buffer.size == expected_buffer.size

    def test_to_pyarrow_numpy_backend(
        self, sample_dataframe_numpy_backend, sample_pyarrow_table
    ):
        result = sample_dataframe_numpy_backend.to_pyarrow()
        assert result == sample_pyarrow_table

    def test_to_pyarrow_pyarrow_backend(
        self, sample_dataframe_pyarrow_backend, sample_pyarrow_table
    ):
        result = sample_dataframe_pyarrow_backend.to_pyarrow()
        assert result == sample_pyarrow_table

    def test_pyarrow_roundtrip(
        self, sample_dataframe_numpy_backend, sample_dataframe_pyarrow_backend
    ):
        result = DataFrame.from_pyarrow(sample_dataframe_numpy_backend.to_pyarrow())
        tm.assert_frame_equal(result, sample_dataframe_pyarrow_backend)
