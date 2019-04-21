# -*- coding: utf-8 -*-

"""
    Tests that skipped data rows are properly handled for all parser,
    building ontop of the skiprows argument
    """

from io import StringIO

import numpy as np
import pytest

from pandas.compat import lrange
from pandas.errors import EmptyDataError

from pandas import DataFrame, Index
import pandas.util.testing as tm


def test_skipdatarows(all_parsers):
    parser = all_parsers
    data = """1.,2,3.
A,B,C
row1,row2,row3
!,?,$
4,5,6
7,8,9"""
    expected = DataFrame([["!", "?", "$"],
                        ["4", "5", "6"],["7", "8", "9"]],
                        columns=["1.", "2", "3."])
    result = parser.read_csv(StringIO(data), skipdatarows=2, delimiter=",")
    tm.assert_frame_equal(result, expected)

def test_comment_skipdatarows(all_parsers):
    parser = all_parsers
    data = """X,Y,Z
first,data,row
row1,row2,row3
A,B,C
D,E,F
#hello world
#ignore this line
last,data,row"""

    result = parser.read_csv(StringIO(data), comment="#", skipdatarows=2,
                             delimiter=",")

    expected = DataFrame([["A", "B", "C"], ["D", "E", "F"],
                            ["last", "data", "row"]],
                                                  columns=["X", "Y", "Z"])

    tm.assert_frame_equal(result, expected)

def test_skiprows_skipdatarows(all_parsers):
    parser = all_parsers
    data = """X,Y,Z
first,data,row
row1,row2,row3
A,B,C
1,2.,4.0
b,c,d
e,f,g
5,NaN,10.0"""
    result = parser.read_csv(StringIO(data), skiprows=2, skipdatarows=1,
                             delimiter=",")

    expected = DataFrame([["row1", "row2", "row3"], ["A", "B", "C"],
                            ["1", "2.", "4.0"],
                            ["b", "c", "d"],
                            ["e", "f", "g"],
                            ["5", np.nan, "10.0"]],
                            columns=["X", "Y", "Z"])
    tm.assert_frame_equal(result, expected)

def test_comment_skiprows_skipdatarows(all_parsers):

    parser = all_parsers
    data = """X,Y,Z
first,data,row
row1,row2,row3
A,B,C
D,E,F
#hello world
#ignore this line
G,H,I
"""
    result = parser.read_csv(StringIO(data), comment="#", skiprows=2,
                            skipdatarows=3, delimiter=",")

    expected = DataFrame([["D", "E", "F"], ["G", "H", "I"]],
                            columns=["X", "Y", "Z"])
    tm.assert_frame_equal(result, expected)
