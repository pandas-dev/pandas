import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.mark.parametrize(
    "df",
    [
        # default RangeIndex
        pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}),
        # named non-range index
        pd.DataFrame({"a": [1, 2, 3]}, index=pd.Index([10, 20, 30], name="idx")),
        # unnamed non-range index
        pd.DataFrame({"a": [1, 2, 3]}, index=[7, 8, 9]),
        # MultiIndex rows
        pd.DataFrame(
            {"a": [1, 2, 3, 4]},
            index=pd.MultiIndex.from_tuples(
                [("x", 1), ("x", 2), ("y", 1), ("y", 2)], names=["l1", "l2"]
            ),
        ),
        # mixed dtypes
        pd.DataFrame(
            {
                "i": [1, 2, 3],
                "f": [1.5, 2.5, 3.5],
                "b": [True, False, True],
                "dt": pd.to_datetime(["2024-01-01", "2024-06-15", "2024-12-31"]),
            }
        ),
    ],
)
def test_from_arrow_roundtrip(df):
    # GH#59780 pandas-side arrow -> DataFrame conversion round-trips through the
    # pandas metadata that pa.Table.from_pandas writes.
    table = pa.Table.from_pandas(df)
    result = pd.DataFrame._from_pyarrow(table)
    tm.assert_frame_equal(result, df)


def test_from_arrow_no_pandas_metadata():
    # a table without the b"pandas" schema metadata converts with a default index
    table = pa.table({"a": [1, 2, 3], "b": ["x", "y", "z"]})
    assert table.schema.pandas_metadata is None
    result = pd.DataFrame._from_pyarrow(table)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})
    tm.assert_frame_equal(result, expected)
