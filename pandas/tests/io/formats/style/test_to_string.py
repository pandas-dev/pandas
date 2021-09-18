from textwrap import dedent

import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
)

pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return DataFrame({"A": [0, 1], "B": [-0.61, -1.22], "C": ["ab", "cd"]})


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0, precision=2)


def test_basic_string(styler):
    result = styler.to_string()
    expected = dedent(
        """\
     A B C
    0 0 -0.61 ab
    1 1 -1.22 cd
    """
    )
    assert result == expected


def test_comprehesive_string():
    midx = MultiIndex.from_product(
        [["A", "B"], ["a", "b"], ["X", "Y"]], names=["zero", "one", "two"]
    )
    cidx = MultiIndex.from_product(
        [["a", "b"], ["C", "D"], ["V", "W"]], names=["zero", "one", "two"]
    )
    df = DataFrame(np.arange(64).reshape(8, 8), index=midx, columns=cidx)
    styler = Styler(df)
    styler.hide_index(level=0).hide_columns(level=1)
    styler.hide_index([("A", "a", "X"), ("A", "b", "X")])
    styler.hide_columns([("a", "C", "W"), ("a", "D", "V")])
    result = styler.to_string()
    assert result == 1
