import pytest

import pandas as pd
import pandas._testing as tm


class TestStrAccessor:
    def test_str_attribute(self):
        # GH#9068
        methods = ["strip", "rstrip", "lstrip"]
        ser = pd.Series([" jack", "jill ", " jesse ", "frank"])
        for method in methods:
            expected = pd.Series([getattr(str, method)(x) for x in ser.values])
            tm.assert_series_equal(getattr(pd.Series.str, method)(ser.str), expected)

        # str accessor only valid with string values
        ser = pd.Series(range(5))
        msg = "Can only use .str accessor with string values, not integer"
        with pytest.raises(AttributeError, match=msg):
            ser.str.repeat(2)

    def test_str_accessor_updates_on_inplace(self):
        ser = pd.Series(list("abc"))
        ser.replace({"a": "A"}, inplace=True)
        assert ser.str.islower().sum() == 2
