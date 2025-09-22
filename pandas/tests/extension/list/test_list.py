import pytest

import pandas as pd
from pandas.tests.extension.list.array import (
    ListArray,
    ListDtype,
    make_data,
)


@pytest.fixture
def dtype():
    return ListDtype()


@pytest.fixture
def data():
    """Length-10 ListArray for semantics test."""
    data = make_data(10)

    while len(data[0]) == len(data[1]):
        data = make_data(10)

    return ListArray(data)


def test_to_csv(data):
    # https://github.com/pandas-dev/pandas/issues/28840
    # array with list-likes fail when doing astype(str) on the numpy array
    # which was done in get_values_for_csv
    df = pd.DataFrame({"a": data})
    res = df.to_csv()
    assert str(data[0]) in res
