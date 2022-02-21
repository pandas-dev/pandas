import pytest

jinja2 = pytest.importorskip("jinja2")

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return DataFrame(
        data=[[0, -0.609], [1, -1.228]],
        columns=["A", "B"],
        index=["x", "y"],
    )


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0)


@pytest.mark.parametrize(
    "kwarg, expected",
    [
        ({"alias": [1, 2, 3]}, "``alias``"),
    ],
)
def test_footer_bad_length(styler, kwarg, expected):
    msg = f"{expected} must have same length as ``func``"
    with pytest.raises(ValueError, match=msg):
        styler.set_footer(func=["mean"], **kwarg)


def test_set_footer_warn():
    df = DataFrame([["a"]])
    styler = df.style.set_footer(["mean"], errors="warn")
    msg = (
        "`Styler.set_footer` raised Exception when calculating method `mean` on "
        "column `0`"
    )
    with tm.assert_produces_warning(Warning, match=msg):
        styler._translate(True, True)


def test_set_footer_raise():
    df = DataFrame([["a"]])
    styler = df.style.set_footer(["mean"], errors="raise")
    msg = (
        "`Styler.set_footer` raised Exception when calculating method `mean` on "
        "column `0`"
    )
    with pytest.raises(Exception, match=msg):
        styler._translate(True, True)
