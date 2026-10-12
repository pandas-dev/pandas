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


def test_append_internal_retains_dtype(data):
    # GH#65431 boxing the values the row pre-cast compares must not read an
    #  array of list-valued scalars as 2-D
    df = pd.DataFrame({"a": data[:2]})
    ser = pd.Series([[5, 6]], index=["a"], name=2)

    result = df._append_internal(ser)

    assert result["a"].dtype == data.dtype
    assert list(result["a"]) == [*list(data[:2]), [5, 6]]
