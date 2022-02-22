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


def test_concat_bad_columns(styler):
    msg = "`other.data` must have same columns as `Styler.data"
    with pytest.raises(ValueError, match=msg):
        styler.concat(DataFrame([[1, 2]]).style)
