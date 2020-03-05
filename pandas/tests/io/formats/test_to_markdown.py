from io import StringIO

import pytest

import pandas as pd

pytest.importorskip("tabulate")


def test_simple():
    buf = StringIO()
    df = pd.DataFrame([1, 2, 3])
    df.to_markdown(buf=buf)
    result = buf.getvalue()
    assert (
        result == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
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
        "|    |   foo |\n|---:|------:|\n|  0 |     1 "
        "|\n|  1 |     2 |\n|  2 |     3 |"
    )


def test_no_buf(capsys):
    df = pd.DataFrame([1, 2, 3])
    result = df.to_markdown()
    assert (
        result == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
    )
