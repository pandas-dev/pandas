from io import StringIO

import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm

pytest.importorskip("tabulate")


def test_keyword_deprecation():
    # GH 57280
    msg = (
        "Starting with pandas version 4.0 all arguments of Series.to_markdown "
        "except for the argument 'buf' will be keyword-only."
    )
    s = pd.Series()
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        s.to_markdown(None, "wt")


def test_simple():
    buf = StringIO()
    df = pd.DataFrame([1, 2, 3])
    df.to_markdown(buf=buf)
    result = buf.getvalue()
    assert (
        result == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
    )


def test_empty_frame():
    buf = StringIO()
    df = pd.DataFrame({"id": [], "first_name": [], "last_name": []}).set_index("id")
    df.to_markdown(buf=buf)
    result = buf.getvalue()
    assert result == (
        "| id   | first_name   | last_name   |\n|------|--------------|-------------|"
    )


def test_other_tablefmt():
    buf = StringIO()
    df = pd.DataFrame([1, 2, 3])
    df.to_markdown(buf=buf, tablefmt="jira")
    result = buf.getvalue()
    assert result == "||    ||   0 ||\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"


def test_other_headers():
    buf = StringIO()
    df = pd.DataFrame([1, 2, 3])
    df.to_markdown(buf=buf, headers=["foo", "bar"])
    result = buf.getvalue()
    assert result == (
        "|   foo |   bar |\n|------:|------:|\n|     0 "
        "|     1 |\n|     1 |     2 |\n|     2 |     3 |"
    )


def test_series():
    buf = StringIO()
    s = pd.Series([1, 2, 3], name="foo")
    s.to_markdown(buf=buf)
    result = buf.getvalue()
    assert result == (
        "|    |   foo |\n|---:|------:|\n|  0 |     1 |\n|  1 |     2 |\n|  2 |     3 |"
    )


def test_no_buf():
    df = pd.DataFrame([1, 2, 3])
    result = df.to_markdown()
    assert (
        result == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
    )


@pytest.mark.parametrize("index", [True, False])
def test_index(index):
    # GH 32667

    df = pd.DataFrame([1, 2, 3])

    result = df.to_markdown(index=index)

    if index:
        expected = (
            "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
        )
    else:
        expected = "|   0 |\n|----:|\n|   1 |\n|   2 |\n|   3 |"
    assert result == expected


def test_showindex_disallowed_in_kwargs():
    # GH 32667; disallowing showindex in kwargs enforced in 2.0
    df = pd.DataFrame([1, 2, 3])
    with pytest.raises(ValueError, match="Pass 'index' instead of 'showindex"):
        df.to_markdown(index=True, showindex=True)


@pytest.mark.parametrize("dtype", ["int64", "uint64"])
@pytest.mark.parametrize("floatfmt", ["g", ".0f"])
def test_large_ints_keep_precision(dtype, floatfmt):
    # GH#49465 a float column alongside makes .values upcast to float64
    df = pd.DataFrame({"a": pd.array([503498111827123021, 2], dtype=dtype)})
    df["b"] = [1.5, 2.5]
    result = df.to_markdown(floatfmt=floatfmt)
    b0, b1 = ("1.5", "2.5") if floatfmt == "g" else ("2", "2")
    expected = (
        "|    |                  a |   b |\n"
        "|---:|-------------------:|----:|\n"
        f"|  0 | 503498111827123021 | {b0:>3} |\n"
        f"|  1 |                  2 | {b1:>3} |"
    )
    assert result == expected


def test_large_ints_keep_precision_single_column():
    # GH#49465
    df = pd.DataFrame({"a": [503498111827123021]})
    result = df.to_markdown(floatfmt=".0f")
    expected = (
        "|    |                  a |\n|---:|-------------------:|\n"
        "|  0 | 503498111827123021 |"
    )
    assert result == expected


@pytest.mark.parametrize("dtype", ["Int64", "int64[pyarrow]"])
def test_nullable_int_missing(dtype):
    # GH#49465 missing values render as tabulate's missingval
    if dtype == "int64[pyarrow]":
        pytest.importorskip("pyarrow")
    df = pd.DataFrame({"a": pd.array([None, 503498111827123021], dtype=dtype)})
    df["b"] = ["x", "y"]
    result = df.to_markdown(missingval="NA")
    expected = (
        "|    |                  a | b   |\n"
        "|---:|-------------------:|:----|\n"
        "|  0 |                 NA | x   |\n"
        "|  1 | 503498111827123021 | y   |"
    )
    assert result == expected
