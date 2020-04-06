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


@pytest.mark.parametrize("index", [True, False, None])
@pytest.mark.parametrize("showindex", [True, False, None])
def test_index(index, showindex):
    # GH 32667
    kwargs = {}
    if index is not None:
        kwargs["index"] = index
    if showindex is not None:
        kwargs["showindex"] = showindex

    df = pd.DataFrame([1, 2, 3])

    if "showindex" in kwargs:
        # force user to use index instead of tabulate's show index
        msg = (
            "'showindex' is not a valid keyword argument, "
            "please use 'index' to control showing the index"
        )
        with pytest.raises(ValueError, match=msg):
            df.to_markdown(**kwargs)
    else:
        result = df.to_markdown(**kwargs)
        if index in [True, None]:
            assert result == (
                "|    |   0 |\n|---:|----:|\n"
                "|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
            )
        else:
            assert result == "|   0 |\n|----:|\n|   1 |\n|   2 |\n|   3 |"
