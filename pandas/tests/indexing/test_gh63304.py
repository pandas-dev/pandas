import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow", minversion="13.0.0")


def test_drop_na_arrow_index():
    # GH#63304
    # Test that dropping pd.NA from PyArrow-backed Index does not raise ArrowInvalid

    # integer
    df = pd.DataFrame(
        {"A": [1, 2, 3]}, index=pd.Index([1, 2, 3], dtype="int64[pyarrow]")
    )
    # pd.NA is not in index, should raise KeyError, but NOT ArrowInvalid
    with pytest.raises(KeyError, match="not found in axis"):
        df.drop(index=[pd.NA])

    # string
    df = pd.DataFrame(
        {"A": [1, 2, 3]}, index=pd.Index(["a", "b", "c"], dtype="string[pyarrow]")
    )
    with pytest.raises(KeyError, match="not found in axis"):
        df.drop(index=[pd.NA])

    # binary
    df = pd.DataFrame(
        {"A": [1, 2, 3]},
        index=pd.Index([b"a", b"b", b"c"], dtype="binary[pyarrow]"),
    )
    with pytest.raises(KeyError, match="not found in axis"):
        df.drop(index=[pd.NA])

    # Case where NA IS in the index (should verify it drops correctly)
    df = pd.DataFrame(
        {"A": [1, 2]}, index=pd.Index([1, pd.NA], dtype="int64[pyarrow]")
    )
    result = df.drop(index=[pd.NA])
    expected = pd.DataFrame({"A": [1]}, index=pd.Index([1], dtype="int64[pyarrow]"))
    tm.assert_frame_equal(result, expected)

    df = pd.DataFrame(
        {"A": [1, 2]}, index=pd.Index(["a", pd.NA], dtype="string[pyarrow]")
    )
    result = df.drop(index=[pd.NA])
    expected = pd.DataFrame(
        {"A": [1]}, index=pd.Index(["a"], dtype="string[pyarrow]")
    )
    tm.assert_frame_equal(result, expected)
