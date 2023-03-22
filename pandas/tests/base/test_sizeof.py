from sys import getsizeof

from pandas import (
    DataFrame,
    Series,
)


def test_sysof():
    assert getsizeof(DataFrame) == 1072
    assert getsizeof(DataFrame()) == 140
    assert getsizeof(DataFrame([])) == 140


def test_sysof_series():
    assert getsizeof(Series) == 1200
    assert getsizeof(Series()) == 140
    getsizeof(Series(str))
    getsizeof(Series(int))
    getsizeof(Series(list))


def test_memory_usage_series():
    Series(str).memory_usage(deep=True)
