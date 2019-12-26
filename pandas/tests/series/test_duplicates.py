import numpy as np
import pytest

from pandas import Categorical, Series
from pandas.core.construction import create_series_with_explicit_dtype
import pandas.util.testing as tm


def test_nunique():
    # basics.rst doc example
    series = Series(np.random.randn(500))
    series[20:500] = np.nan
    series[10:20] = 5000
    result = series.nunique()
    assert result == 11

    # GH 18051
    s = Series(Categorical([]))
    assert s.nunique() == 0
    s = Series(Categorical([np.nan]))
    assert s.nunique() == 0


def test_unique():
    # GH714 also, dtype=float
    s = Series([1.2345] * 100)
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    s = Series([1.2345] * 100, dtype="f4")
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    # NAs in object arrays #714
    s = Series(["foo"] * 100, dtype="O")
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    # decision about None
    s = Series([1, 2, 3, None, None, None], dtype=object)
    result = s.unique()
    expected = np.array([1, 2, 3, None], dtype=object)
    tm.assert_numpy_array_equal(result, expected)

    # GH 18051
    s = Series(Categorical([]))
    tm.assert_categorical_equal(s.unique(), Categorical([]), check_dtype=False)
    s = Series(Categorical([np.nan]))
    tm.assert_categorical_equal(s.unique(), Categorical([np.nan]), check_dtype=False)


def test_unique_data_ownership():
    # it works! #1807
    Series(Series(["a", "c", "b"]).unique()).sort_values()


@pytest.mark.parametrize(
    "data, expected",
    [
        (np.random.randint(0, 10, size=1000), False),
        (np.arange(1000), True),
        ([], True),
        ([np.nan], True),
        (["foo", "bar", np.nan], True),
        (["foo", "foo", np.nan], False),
        (["foo", "bar", np.nan, np.nan], False),
    ],
)
def test_is_unique(data, expected):
    # GH11946 / GH25180
    s = create_series_with_explicit_dtype(data, dtype_if_empty=object)
    assert s.is_unique is expected


def test_is_unique_class_ne(capsys):
    # GH 20661
    class Foo:
        def __init__(self, val):
            self._value = val

        def __ne__(self, other):
            raise Exception("NEQ not supported")

    with capsys.disabled():
        li = [Foo(i) for i in range(5)]
        s = Series(li, index=list(range(5)))
    s.is_unique
    captured = capsys.readouterr()
    assert len(captured.err) == 0
