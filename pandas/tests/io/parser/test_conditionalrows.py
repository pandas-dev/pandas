"""
Tests that skipped rows are properly handled during
parsing for all of the parsers defined in parsers.py
"""

from datetime import datetime
from io import StringIO

import numpy as np
import pytest

from pandas.errors import EmptyDataError

from pandas import DataFrame, Index
import pandas._testing as tm


def test_conditional_rows_single_column_less_than(all_parsers):
    # see gh-32072
    parser = all_parsers
    data = """country,capital,area,population
Brazil,Brasilia,8.516,200.4
Russia,Moscow,17.10,143.5
India,New Delhi,3.286,1252
China,Beijing,9.597,1357
South Africa,Pretoria,1.221,52.98
"""
    df = parser.read_csv(StringIO(data), conditionalrows="(area > 8.516)")
    expected = DataFrame(
        data={
            "country": ["Russia", "China"],
            "capital": ["Moscow", "Beijing"],
            "area": [17.10, 9.597],
            "population": [143.5, 1357],
        }
    )
    tm.assert_frame_equal(df, expected)


def test_conditional_rows_single_column_greater_than(all_parsers):
    # see gh-32072
    parser = all_parsers
    data = """country,capital,area,population
Brazil,Brasilia, 8.516, 200.4
Russia,Moscow,17.10,143.5
India,New Delhi,3.286,1252
China,Beijing,9.597,1357
South Africa,Pretoria,1.221,52.98
"""
    df = parser.read_csv(StringIO(data), conditionalrows="(area < 8.516)")
    expected = DataFrame(
        data={
            "country": ["India", "South Africa"],
            "capital": ["New Delhi", "Pretoria"],
            "area": [3.286, 1.221],
            "population": [1252, 52.98],
        }
    )
    tm.assert_frame_equal(df, expected)


def test_conditional_rows_multi_columns_and(all_parsers):
    # see gh-32072
    parser = all_parsers
    data = """country,capital,area,population
Brazil,Brasilia, 8.516, 200.4
Russia,Moscow,17.10,143.5
India,New Delhi,3.286,1252
China,Beijing,9.597,1357
South Africa,Pretoria,1.221,52.98
"""
    df = parser.read_csv(
        StringIO(data), conditionalrows="(area <= 8.516 and population > 1200)"
    )
    expected = DataFrame(
        data={
            "country": ["India"],
            "capital": ["New Delhi"],
            "area": [3.286],
            "population": [1252.0],
        }
    )
    tm.assert_frame_equal(df, expected)


def test_conditional_rows_multi_columns_or(all_parsers):
    # see gh-32072
    parser = all_parsers
    data = """country,capital,area,population
Brazil,Brasilia, 8.516, 200.4
Russia,Moscow,17.10,143.5
India,New Delhi,3.286,1252
China,Beijing,9.597,1357
South Africa,Pretoria,1.221,52.98
"""
    df = parser.read_csv(
        StringIO(data), conditionalrows="(area > 8.516 or area < 3.286)"
    )
    expected = DataFrame(
        data={
            "country": ["Russia", "China", "South Africa"],
            "capital": ["Moscow", "Beijing", "Pretoria"],
            "area": [17.10, 9.597, 1.221],
            "population": [143.5, 1357.00, 52.98],
        }
    )
    tm.assert_frame_equal(df, expected)
