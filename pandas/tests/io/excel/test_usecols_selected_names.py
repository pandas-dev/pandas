"""Regression coverage for Excel positional usecols with duplicate names."""

from io import BytesIO

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "headers, usecols, positions, expected_names",
    [
        (["A", "B", "C", "B"], "C:D", [2, 3], ["C", "B"]),
        (["A", "B", "C", "B"], "D:D", [3], ["B"]),
        (["A", "B", "C", "B"], "B:D", [1, 2, 3], ["B", "C", "B.1"]),
        (["A", "B", "C", "B"], "A:D", [0, 1, 2, 3], ["A", "B", "C", "B.1"]),
        (["a", "a", "a.1"], "B:C", [1, 2], ["a", "a.1"]),
        (["a", "a", "a.1"], "B:B", [1], ["a"]),
    ],
)
def test_read_excel_usecols_mangles_selected_names_only(
    headers, usecols, positions, expected_names
):
    # GH#64198: unselected headers must not affect deduplication.
    pytest.importorskip("openpyxl")

    original = pd.DataFrame(
        [list(range(len(headers))), list(range(10, 10 + len(headers)))],
        columns=headers,
    )
    file = BytesIO()
    original.to_excel(file, index=False, engine="openpyxl")
    file.seek(0)

    result = pd.read_excel(file, usecols=usecols, engine="openpyxl")
    expected = original.iloc[:, positions].copy()
    expected.columns = expected_names
    tm.assert_frame_equal(result, expected)


def test_read_excel_usecols_callable_unchanged():
    # Callable selection uses the original header names before filtering.
    pytest.importorskip("openpyxl")

    original = pd.DataFrame([[1, 2, 3, 4]], columns=["A", "B", "C", "B"])
    file = BytesIO()
    original.to_excel(file, index=False, engine="openpyxl")
    file.seek(0)

    result = pd.read_excel(
        file, usecols=lambda col: col in ("B", "C"), engine="openpyxl"
    )
    expected = pd.DataFrame({"B": [2], "C": [3]})
    tm.assert_frame_equal(result, expected)


def test_read_excel_usecols_label_based_keeps_full_name_resolution():
    # Label-based usecols is deliberately unchanged by GH#64198.
    pytest.importorskip("openpyxl")

    original = pd.DataFrame([[1, 2, 3]], columns=["A", "A", "B"])
    file = BytesIO()
    original.to_excel(file, index=False, engine="openpyxl")
    file.seek(0)

    result = pd.read_excel(file, usecols=["A.1"], engine="openpyxl")
    expected = pd.DataFrame({"A.1": [2]})
    tm.assert_frame_equal(result, expected)
