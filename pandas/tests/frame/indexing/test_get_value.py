import pytest
from enum import Enum
from pandas import (
    DataFrame,
    MultiIndex,
)


class TestGetValue:
    def test_get_set_value_no_partial_indexing(self):
        # partial w/ MultiIndex raise exception
        index = MultiIndex.from_tuples([(0, 1), (0, 2), (1, 1), (1, 2)])
        df = DataFrame(index=index, columns=range(4))
        with pytest.raises(KeyError, match=r"^0$"):
            df._get_value(0, 1)

    def test_get_value(self, float_frame):
        for idx in float_frame.index:
            for col in float_frame.columns:
                result = float_frame._get_value(idx, col)
                expected = float_frame[col][idx]
                assert result == expected

def test_enum_value(self,):
    Cols = Enum('Cols', 'col1 col2')

    q1 = DataFrame({Cols.col1: [1, 2, 3]})
    q2 = DataFrame({Cols.col1: [1, 2, 3]})

    result=((q1[Cols.col1] == q2[Cols.col1]).all())

    expected=True

    assert result == expected

