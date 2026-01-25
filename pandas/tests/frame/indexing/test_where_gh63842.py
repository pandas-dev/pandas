import numpy as np
import pytest

from pandas import DataFrame, Series, StringDtype
import pandas._testing as tm


def test_where_stringdtype_preservation_with_list_like():
    # GH#63842
    df = DataFrame({"A": Series(["x", "y"], dtype="string")})
    mask = DataFrame({"A": [True, False]})
    result = df.where(mask, ["z", "w"])

    expected = DataFrame({"A": Series(["x", "w"], dtype="string")})
    tm.assert_frame_equal(result, expected)
