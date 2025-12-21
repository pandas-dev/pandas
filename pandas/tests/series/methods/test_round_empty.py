import pytest
import numpy as np
import pandas as pd
import pandas.testing as tm

#GH63444


def test_series_round_empty1():
    result = pd.Series([]).round(2)
    expected = pd.Series([])

def test_series_round_empty2():
    result = pd.Series([], dtype='float64').round()
    expected = pd.Series([], dtype='float64')
    tm.assert_series_equal(result, expected)
def test_series_round_numeric():
    result = pd.Series([1.2345, 2.6789, 3.14159]).round(2)
    expected = pd.Series([1.23, 2.68, 3.14])
    tm.assert_series_equal(result, expected)
def test_series_round_integers():
    result = pd.Series([1, 2, 3]).round(2)
    expected = pd.Series([1, 2, 3])
    tm.assert_series_equal(result, expected)

def test_series_round_str():
    with pytest.raises(TypeError):
        pd.Series(['a', 'b', 'c']).round(2)
def test_series_round_obj():
    with pytest.raises(TypeError):
        pd.Series(['123'], dtype='object').round(2)
