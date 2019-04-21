# -*- coding: utf-8 -*-

"""
Tests that skipped data rows are properly handled for all parser, building ontop of the skiprows argument
"""

from io import StringIO

import numpy as np
import pytest

from pandas.compat import lrange
from pandas.errors import EmptyDataError

from pandas import DataFrame, Index
import pandas.util.testing as tm

@pytest.mark.parametrize("data", [
"""1.,2,3.
A,B,C
row1,row2,row3
!,?,$
4,5,6
7,8,9
"""])
def test_skipdatarows(all_parsers):
    parser = all_parsers
    expected = DataFrame(["row1", "row2", "row3"], ["!", "?", "$"],
                            [4, 5, 6], [7, 8, 9]],
                            columns=["row1", "row2", "row3"])
    result = parser.read_csv(StringIO(data), skipdatarows=2, delimiter=",")
    tm.assert_frame_equal(result, expected)

@pytest.mark.parametrize("data", [
"""X,Y,Z
first,data,row
row1,row2,row3
A,B,C
1,2.,4.
#hello world
#ignore this line
5.,NaN,10.0"""])
def test_comment_skipdatarows(all_parsers):
    parser = all_parsers
    result = parser.read_csv(StringIO(data), comment="#", skipdatarows=2,
                            delimeter=",")

    expected = DataFrame([["A", "B", "C"], [1, 2., 4.], [5., np.nan, 10.]],
                            columns=["X", "Y", "Z"])

    tm.assert_frame_equal(result, expected)

def test_skiprows_skipdatarows(all_parsers):
    parser = all_parsers
    result = parser.read_csv(StringIO(data), skiprows=2, skipdatarows=1,
                            delimeter=",")

    expected = DataFrame([["row1", "row2", "row3"], ["A", "B", "C"],
                            [1, 2., 4.], ["#hello world", np.nan, np.nan],
                            ["#ignore this line", np.nan, np.nan],
                            [5., np.nan, 10.]],
                            columns=["X", "Y", "Z"])

    tm.assert_frame_equal(result, expected)

def test_comment_skiprows_skipdatarows():

    parser = all_parsers
    result = parser.read_csv(StringIO(data), comment="#", skiprows=2,
                            skipdatarows=3, delimeter=",")

    expected = DataFrame([[1, 2., 4.], [5., np.nan, 10.]],
                            columns=["X", "Y", "Z"])

    tm.assert_frame_equal(result, expected)
