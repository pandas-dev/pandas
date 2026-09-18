"""
Tests dtype specification during parsing
for all of the parsers defined in parsers.py
"""

from io import StringIO

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_dtype_all_columns_empty(all_parsers):
    # see gh-12048
    parser = all_parsers
    result = parser.read_csv(StringIO("A,B"), dtype=str)

    expected = pd.DataFrame({"A": [], "B": []}, dtype=str)
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two"
    result = parser.read_csv(StringIO(data), dtype={"one": "u1"})

    expected = pd.DataFrame(
        {"one": np.empty(0, dtype="u1"), "two": np.empty(0, dtype=object)},
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_with_index_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two"
    result = parser.read_csv(
        StringIO(data), index_col=["one"], dtype={"one": "u1", 1: "f"}
    )

    expected = pd.DataFrame(
        {"two": np.empty(0, dtype="f")}, index=pd.Index([], dtype="u1", name="one")
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_with_multi_index_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two,three"
    result = parser.read_csv(
        StringIO(data), index_col=["one", "two"], dtype={"one": "u1", 1: "f8"}
    )

    exp_idx = pd.MultiIndex.from_arrays(
        [np.empty(0, dtype="u1"), np.empty(0, dtype=np.float64)],
        names=["one", "two"],
    )
    expected = pd.DataFrame({"three": np.empty(0, dtype=object)}, index=exp_idx)
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_with_mangled_column_pass_dtype_by_names(all_parsers):
    parser = all_parsers

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={"one": "u1", "one.1": "f"})

    expected = pd.DataFrame(
        {"one": np.empty(0, dtype="u1"), "one.1": np.empty(0, dtype="f")},
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_with_mangled_column_pass_dtype_by_indexes(all_parsers):
    parser = all_parsers

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={0: "u1", 1: "f"})

    expected = pd.DataFrame(
        {"one": np.empty(0, dtype="u1"), "one.1": np.empty(0, dtype="f")},
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_with_dup_column_pass_dtype_by_indexes(all_parsers):
    # see gh-9424
    parser = all_parsers
    expected = pd.concat(
        [pd.Series([], name="one", dtype="u1"), pd.Series([], name="one.1", dtype="f")],
        axis=1,
    )

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={0: "u1", 1: "f"})
    tm.assert_frame_equal(result, expected)


def test_empty_with_dup_column_pass_dtype_by_indexes_raises(all_parsers):
    # see gh-9424
    parser = all_parsers
    expected = pd.concat(
        [pd.Series([], name="one", dtype="u1"), pd.Series([], name="one.1", dtype="f")],
        axis=1,
    )
    expected.index = expected.index.astype(object)

    with pytest.raises(ValueError, match="Duplicate names"):
        parser.read_csv(StringIO(""), names=["one", "one"], dtype={0: "u1", 1: "f"})


@pytest.mark.parametrize(
    "dtype,expected",
    [
        (np.float64, pd.DataFrame(columns=["a", "b"], dtype=np.float64)),
        (
            "category",
            pd.DataFrame({"a": pd.Categorical([]), "b": pd.Categorical([])}),
        ),
        (
            {"a": "category", "b": "category"},
            pd.DataFrame({"a": pd.Categorical([]), "b": pd.Categorical([])}),
        ),
        ("datetime64[ns]", pd.DataFrame(columns=["a", "b"], dtype="datetime64[ns]")),
        (
            "timedelta64[ns]",
            pd.DataFrame(
                {
                    "a": pd.Series([], dtype="timedelta64[ns]"),
                    "b": pd.Series([], dtype="timedelta64[ns]"),
                },
            ),
        ),
        (
            {"a": np.int64, "b": np.int32},
            pd.DataFrame(
                {
                    "a": pd.Series([], dtype=np.int64),
                    "b": pd.Series([], dtype=np.int32),
                },
            ),
        ),
        (
            {0: np.int64, 1: np.int32},
            pd.DataFrame(
                {
                    "a": pd.Series([], dtype=np.int64),
                    "b": pd.Series([], dtype=np.int32),
                },
            ),
        ),
        (
            {"a": np.int64, 1: np.int32},
            pd.DataFrame(
                {
                    "a": pd.Series([], dtype=np.int64),
                    "b": pd.Series([], dtype=np.int32),
                },
            ),
        ),
    ],
)
@skip_pyarrow  # CSV parse error: Empty CSV file or block
def test_empty_dtype(all_parsers, dtype, expected):
    # see gh-14712
    parser = all_parsers
    data = "a,b"

    result = parser.read_csv(StringIO(data), header=0, dtype=dtype)
    tm.assert_frame_equal(result, expected)
