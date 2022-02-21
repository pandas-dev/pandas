import pytest

jinja2 = pytest.importorskip("jinja2")

from pandas import DataFrame

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
