"""
modules collects tests for Styler methods which have been deprecated
"""
import numpy as np
import pytest

jinja2 = pytest.importorskip("jinja2")

from pandas import DataFrame
import pandas._testing as tm


@pytest.fixture
def df():
    return DataFrame({"A": [0, 1], "B": np.random.randn(2)})


def test_render(df):
    with tm.assert_produces_warning(FutureWarning):
        df.style.render()
