import numpy as np
import pytest

from pandas.errors import DtypeWarning

import pandas._testing as tm
from pandas.core.arrays import ArrowExtensionArray

from pandas.io.parsers.c_parser_wrapper import _concatenate_chunks


def test_concatenate_chunks_pyarrow():
    # GH#51876
    pa = pytest.importorskip("pyarrow")
    chunks = [
        {0: ArrowExtensionArray(pa.array([1.5, 2.5]))},
        {0: ArrowExtensionArray(pa.array([1, 2]))},
    ]
    result = _concatenate_chunks(chunks, ["column_0", "column_1"])
    expected = ArrowExtensionArray(pa.array([1.5, 2.5, 1.0, 2.0]))
    tm.assert_extension_array_equal(result[0], expected)


def test_concatenate_chunks_pyarrow_strings():
    # GH#51876
    pa = pytest.importorskip("pyarrow")
    chunks = [
        {0: ArrowExtensionArray(pa.array([1.5, 2.5]))},
        {0: ArrowExtensionArray(pa.array(["a", "b"]))},
    ]
    with tm.assert_produces_warning(
        DtypeWarning, match="Columns \\(0: column_0\\) have mixed types"
    ):
        result = _concatenate_chunks(chunks, ["column_0", "column_1"])
    expected = np.concatenate(
        [np.array([1.5, 2.5], dtype=object), np.array(["a", "b"])]
    )
    tm.assert_numpy_array_equal(result[0], expected)


def test_concatenate_chunks_usecols_positions():
    # GH#67375 the chunks are keyed by the column's position in the file, which
    # is not an index into column_names once usecols has narrowed it
    pa = pytest.importorskip("pyarrow")
    chunks = [
        {50: ArrowExtensionArray(pa.array([1.5, 2.5]))},
        {50: ArrowExtensionArray(pa.array(["a", "b"]))},
    ]
    with tm.assert_produces_warning(
        DtypeWarning, match="Columns \\(50: column_50\\) have mixed types"
    ):
        result = _concatenate_chunks(chunks, ["column_50"])
    expected = np.concatenate(
        [np.array([1.5, 2.5], dtype=object), np.array(["a", "b"])]
    )
    tm.assert_numpy_array_equal(result[50], expected)
