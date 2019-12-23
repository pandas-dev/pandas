import numpy as np
import pytest

import pandas as pd
from pandas import Series, date_range
import pandas.util.testing as tm


class TestSeriesIsIn:
    def test_isin(self):
        s = Series(["A", "B", "C", "a", "B", "B", "A", "C"])

        result = s.isin(["A", "C"])
        expected = Series([True, False, True, False, False, False, True, True])
        tm.assert_series_equal(result, expected)

        # GH#16012
        # This specific issue has to have a series over 1e6 in len, but the
        # comparison array (in_list) must be large enough so that numpy doesn't
        # do a manual masking trick that will avoid this issue altogether
        s = Series(list("abcdefghijk" * 10 ** 5))
        # If numpy doesn't do the manual comparison/mask, these
        # unorderable mixed types are what cause the exception in numpy
        in_list = [-1, "a", "b", "G", "Y", "Z", "E", "K", "E", "S", "I", "R", "R"] * 6

        assert s.isin(in_list).sum() == 200000

    def test_isin_with_string_scalar(self):
        # GH#4763
        s = Series(["A", "B", "C", "a", "B", "B", "A", "C"])
        msg = (
            r"only list-like objects are allowed to be passed to isin\(\),"
            r" you passed a \[str\]"
        )
        with pytest.raises(TypeError, match=msg):
            s.isin("a")

        s = Series(["aaa", "b", "c"])
        with pytest.raises(TypeError, match=msg):
            s.isin("aaa")

    def test_isin_with_i8(self):
        # GH#5021

        expected = Series([True, True, False, False, False])
        expected2 = Series([False, True, False, False, False])

        # datetime64[ns]
        s = Series(date_range("jan-01-2013", "jan-05-2013"))

        result = s.isin(s[0:2])
        tm.assert_series_equal(result, expected)

        result = s.isin(s[0:2].values)
        tm.assert_series_equal(result, expected)

        # fails on dtype conversion in the first place
        result = s.isin(s[0:2].values.astype("datetime64[D]"))
        tm.assert_series_equal(result, expected)

        result = s.isin([s[1]])
        tm.assert_series_equal(result, expected2)

        result = s.isin([np.datetime64(s[1])])
        tm.assert_series_equal(result, expected2)

        result = s.isin(set(s[0:2]))
        tm.assert_series_equal(result, expected)

        # timedelta64[ns]
        s = Series(pd.to_timedelta(range(5), unit="d"))
        result = s.isin(s[0:2])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("empty", [[], Series(dtype=object), np.array([])])
    def test_isin_empty(self, empty):
        # see GH#16991
        s = Series(["a", "b"])
        expected = Series([False, False])

        result = s.isin(empty)
        tm.assert_series_equal(expected, result)
