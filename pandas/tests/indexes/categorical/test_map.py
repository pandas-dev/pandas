import pytest

import pandas as pd
from pandas import CategoricalIndex
import pandas._testing as tm


class TestMap:
    @pytest.mark.parametrize(
        "data, categories",
        [
            (list("abcbca"), list("cab")),
            (pd.interval_range(0, 3).repeat(3), pd.interval_range(0, 3)),
        ],
        ids=["string", "interval"],
    )
    def test_map_str(self, data, categories, ordered_fixture):
        # GH 31202 - override base class since we want to maintain categorical/ordered
        index = CategoricalIndex(data, categories=categories, ordered=ordered_fixture)
        result = index.map(str)
        expected = CategoricalIndex(
            map(str, data), categories=map(str, categories), ordered=ordered_fixture
        )
        tm.assert_index_equal(result, expected)
