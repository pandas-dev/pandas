"""Examples covering terminal fitting against explicit column selections."""

import pytest

import pandas as pd
import pandas._testing as tm
from pandas.io.formats.printing import get_adjustment


@pytest.mark.parametrize(
    "example",
    [
        "url",
        "wide_header",
        "unicode",
        "missing",
        "formatter",
        "multiindex",
        "row_truncation",
    ],
)
def test_terminal_width_examples(monkeypatch, example):
    width = 40
    index = example in {"unicode", "multiindex"}
    terminal = (width, 24)
    monkeypatch.setattr("pandas.io.formats.format.get_terminal_size", lambda: terminal)
    monkeypatch.setattr("pandas.io.formats.string.get_terminal_size", lambda: terminal)
    df = pd.DataFrame([["x"] * 7 for _ in range(8)], columns=list("abcdefg"))
    options = {"index": index}
    if example == "url":
        df["e"] = "https://example.com/reports/annual-report.csv"
    elif example == "wide_header":
        df.columns = ["a", "b", "c", "d", "a very long descriptive heading", "f", "g"]
    elif example == "unicode":
        df["e"] = "界" * 20
        df.index = pd.Index(["行" + str(i) for i in range(8)], name="番号")
    elif example == "missing":
        df = pd.DataFrame({c: pd.array([1, None, 3], dtype="Int64") for c in df})
        options["na_rep"] = "missing value with a long description"
    elif example == "formatter":
        options["formatters"] = {"e": lambda value: "formatted value " * 4}
    elif example == "multiindex":
        df["e"] = "long value " * 8
        df.columns = pd.MultiIndex.from_tuples([("group", c) for c in df.columns])
        df.index = pd.MultiIndex.from_product([["A", "B"], range(4)])
    elif example == "row_truncation":
        df["e"] = "long value " * 8
        options["max_rows"] = 4

    before = df.copy(deep=True)
    with pd.option_context("display.unicode.east_asian_width", example == "unicode"):
        adjustment = get_adjustment()
        # Explicit truncation provides independently rendered alternatives.
        candidates = [df.to_string(**options)]
        candidates += [df.to_string(max_cols=n, **options) for n in [6, 4, 2]]
        expected = next(
            (
                text
                for text in candidates
                if max(map(adjustment.len, text.splitlines())) < width
            ),
            candidates[-1],
        )
        result = df.to_string(max_cols=0, line_width=width, **options)

    assert result == expected
    tm.assert_frame_equal(df, before)
