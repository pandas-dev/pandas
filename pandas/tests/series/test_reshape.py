import pytest

import pandas as pd
import pandas._testing as tm


def test_unstack_tuplename_in_multiindex():
    # GH 19966
    idx = pd.MultiIndex.from_product(
        [["a", "b", "c"], [1, 2, 3]], names=[("A", "a"), ("B", "b")]
    )
    ser = pd.Series(1, index=idx)
    result = ser.unstack(("A", "a"))

    expected = pd.DataFrame(
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
        columns=pd.MultiIndex.from_tuples(
            [("a",), ("b",), ("c",)], names=[("A", "a")],
        ),
        index=pd.Index([1, 2, 3], name=("B", "b")),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "unstack_idx, expected_values, expected_index, expected_columns",
    [
        (
            ("A", "a"),
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            pd.MultiIndex.from_tuples(
                [(1, 3), (1, 4), (2, 3), (2, 4)], names=["B", "C"]
            ),
            pd.MultiIndex.from_tuples([("a",), ("b",)], names=[("A", "a")]),
        ),
        (
            (("A", "a"), "B"),
            [[1, 1, 1, 1], [1, 1, 1, 1]],
            pd.Index([3, 4], name="C"),
            pd.MultiIndex.from_tuples(
                [("a", 1), ("a", 2), ("b", 1), ("b", 2)], names=[("A", "a"), "B"]
            ),
        ),
    ],
)
def test_unstack_mixed_type_name_in_multiindex(
    unstack_idx, expected_values, expected_index, expected_columns
):
    # GH 19966
    idx = pd.MultiIndex.from_product(
        [["a", "b"], [1, 2], [3, 4]], names=[("A", "a"), "B", "C"]
    )
    ser = pd.Series(1, index=idx)
    result = ser.unstack(unstack_idx)

    expected = pd.DataFrame(
        expected_values, columns=expected_columns, index=expected_index,
    )
    tm.assert_frame_equal(result, expected)
