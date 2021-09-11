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

        result = (
            m.merge(
                left,
                right,
                on=lambda l, r: (r.timestart <= l.timestep) & (l.timestep <= r.timeend),
            )
            .sort_values(["timestep", "mood", "timestart", "timeend"])
            .reset_index(drop=True)
        )
        expected = (
            m.merge(left, right, how="cross")
            .loc[lambda l, r: (r.timestart <= l.timestep) & (l.timestep <= r.timeend)]
            .sort_values(["timestep", "mood", "timestart", "timeend"])
            .reset_index(drop=True)
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(left, left_copy)
        tm.assert_frame_equal(right, right_copy)


@pytest.mark.parametrize("how", ["left", "right", "outer"])
def test_merge_conditional_non_cross(how):
    error_msg = (
        '`Conditional merge is currently only available for how="inner". '
        "Other merge types will be available in a future version."
    )
    with pytest.raises(NotImplementedError, match=error_msg):
        m.merge(DataFrame(), DataFrame(), on=lambda dfx: None, how=how)
