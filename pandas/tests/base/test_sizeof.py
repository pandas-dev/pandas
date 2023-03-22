from sys import getsizeof

from pandas import (
    DataFrame,
    Series,
)


def test_sysof():
    getsizeof(DataFrame)
    getsizeof(DataFrame())
    getsizeof(DataFrame([]))


def test_sysof_series():
    getsizeof(Series(str))
    getsizeof(Series(int))
    getsizeof(Series(list))


def test_memory_usage_series():
    Series(str).memory_usage(deep=True)
