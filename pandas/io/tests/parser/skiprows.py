# -*- coding: utf-8 -*-

"""
Tests that skipped rows are properly handled during
parsing for all of the parsers defined in parsers.py
"""

from datetime import datetime

import numpy as np

import pandas.util.testing as tm

from pandas import DataFrame
from pandas.compat import StringIO, range, lrange


class SkipRowsTests(object):

    def test_skiprows_bug(self):
        # see gh-505
        text = """#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
1/1/2000,1.,2.,3.
1/2/2000,4,5,6
1/3/2000,7,8,9
"""
        data = self.read_csv(StringIO(text), skiprows=lrange(6), header=None,
                             index_col=0, parse_dates=True)

        data2 = self.read_csv(StringIO(text), skiprows=6, header=None,
                              index_col=0, parse_dates=True)

        expected = DataFrame(np.arange(1., 10.).reshape((3, 3)),
                             columns=[1, 2, 3],
                             index=[datetime(2000, 1, 1), datetime(2000, 1, 2),
                                    datetime(2000, 1, 3)])
        expected.index.name = 0
        tm.assert_frame_equal(data, expected)
        tm.assert_frame_equal(data, data2)

    def test_deep_skiprows(self):
        # see gh-4382
        text = "a,b,c\n" + \
               "\n".join([",".join([str(i), str(i + 1), str(i + 2)])
                          for i in range(10)])
        condensed_text = "a,b,c\n" + \
                         "\n".join([",".join([str(i), str(i + 1), str(i + 2)])
                                    for i in [0, 1, 2, 3, 4, 6, 8, 9]])
        data = self.read_csv(StringIO(text), skiprows=[6, 8])
        condensed_data = self.read_csv(StringIO(condensed_text))
        tm.assert_frame_equal(data, condensed_data)

    def test_skiprows_blank(self):
        # see gh-9832
        text = """#foo,a,b,c
#foo,a,b,c

#foo,a,b,c
#foo,a,b,c

1/1/2000,1.,2.,3.
1/2/2000,4,5,6
1/3/2000,7,8,9
"""
        data = self.read_csv(StringIO(text), skiprows=6, header=None,
                             index_col=0, parse_dates=True)

        expected = DataFrame(np.arange(1., 10.).reshape((3, 3)),
                             columns=[1, 2, 3],
                             index=[datetime(2000, 1, 1), datetime(2000, 1, 2),
                                    datetime(2000, 1, 3)])
        expected.index.name = 0
        tm.assert_frame_equal(data, expected)
