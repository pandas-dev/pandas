import pytest

from pandas import DataFrame
import pandas._testing as tm
import pandas.core.reshape.merge as m


# choose chunk sizes such that the merges will happen with a single
# chunk or multiple chunks
@pytest.mark.parametrize("chunk_size", [2, 4, 10])
def test_merge_conditional(monkeypatch, chunk_size):
    # GH#8962

    with monkeypatch.context() as context:
        context.setattr(
            "pandas.core.reshape.merge._DEFAULT_LEFT_CHUNK_SIZE", chunk_size
        )
        context.setattr(
            "pandas.core.reshape.merge._DEFAULT_RIGHT_CHUNK_SIZE", chunk_size
        )

        left = DataFrame({"timestep": range(5)})
        right = DataFrame(
            {
                "mood": ["happy", "jolly", "joy", "cloud9"],
                "timestart": [0, 2, 2, 3],
                "timeend": [1, 3, 4, 4],
            }
        )
        left_copy = left.copy()
        right_copy = right.copy()

        def condition(dfx):
            return (dfx.timestart <= dfx.timestep) & (dfx.timestep <= dfx.timeend)

        result = (
            m.merge(left, right, condition=condition)
            .sort_values(["timestep", "mood", "timestart", "timeend"])
            .reset_index(drop=True)
        )
        expected = (
            m.merge(left, right, how="cross")
            .loc[condition]
            .sort_values(["timestep", "mood", "timestart", "timeend"])
            .reset_index(drop=True)
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(left, left_copy)
        tm.assert_frame_equal(right, right_copy)


def test_merge_conditional_non_cross():
    error_msg = 'Must use `how="cross" | None` if `condition` is specified'
    with pytest.raises(m.MergeError, match=error_msg):
        m.merge(DataFrame(), DataFrame(), condition=lambda dfx: None, how="inner")
