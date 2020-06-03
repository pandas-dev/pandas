""" test get/set & misc """


import pandas as pd
from pandas import MultiIndex, Series


def test_access_none_value_in_multiindex():
    # GH34318: test that you can access a None value using .loc through a Multiindex

    s = Series([None], pd.MultiIndex.from_arrays([["Level1"], ["Level2"]]))
    result = s.loc[("Level1", "Level2")]
    assert result is None

    midx = MultiIndex.from_product([["Level1"], ["Level2_a", "Level2_b"]])
    s = Series([None] * len(midx), dtype=object, index=midx)
    result = s.loc[("Level1", "Level2_a")]
    assert result is None

    s = Series([1] * len(midx), dtype=object, index=midx)
    result = s.loc[("Level1", "Level2_a")]
    assert result == 1
