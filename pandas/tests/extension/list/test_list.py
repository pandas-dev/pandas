import pytest

import pandas as pd
from pandas.core.arrays.list_ import (
    ListArray,
    ListDtype,
)
from pandas.tests.extension.base.constructors import BaseConstructorsTests


@pytest.fixture
def dtype():
    return ListDtype()


@pytest.fixture
def data():
    """Length-100 ListArray for semantics test."""
    # TODO: make better random data
    data = [list("a"), list("ab"), list("abc")] * 33 + [None]
    return ListArray(data)


class TestListArray(BaseConstructorsTests): ...


def test_to_csv(data):
    # https://github.com/pandas-dev/pandas/issues/28840
    # array with list-likes fail when doing astype(str) on the numpy array
    # which was done in get_values_for_csv
    df = pd.DataFrame({"a": data})
    res = df.to_csv()
    assert str(data[0]) in res
