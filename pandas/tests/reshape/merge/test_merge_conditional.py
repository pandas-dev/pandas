import pytest

from pandas import DataFrame, option_context
from pandas.core.reshape.merge import merge
import pandas._testing as tm


# choose chunk sizes such that the merges will happen with a single
# chunk or multiple chunks
@pytest.mark.parametrize("chunk_size", [2, 4, 10])
def test_merge_conditional(chunk_size):
    # GH#8962

    with option_context(
        "conditional_merge.left_chunk_size", chunk_size,
        "conditional_merge.right_chunk_size", chunk_size,
    ):
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
            merge(
                left,
                right,
                on=lambda l, r: (r.timestart <= l.timestep) & (l.timestep <= r.timeend),
            )
            .sort_values(["timestep", "mood", "timestart", "timeend"])
            .reset_index(drop=True)
        )
        expected = (
            merge(left, right, how="cross")
            .loc[
                lambda dfx: (dfx.timestart <= dfx.timestep)
                & (dfx.timestep <= dfx.timeend)
            ]
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
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, how=how)


def test_merge_conditional_bad_args():
    error_msg = (
        "Cannot define any of (`left_on`, `right_on`, `left_index`, "
        "`right_index`) in a conditional merge"
    )
    with pytest.raises(ValueError, match=error_msg):
        merge(
            DataFrame(columns=['A']),
            DataFrame(columns=['A']),
            on=lambda l, r: None,
            left_on='A'
        )

    error_msg = "Conditional merge must use `axis=1`"
    with pytest.raises(ValueError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, axis=0)

    error_msg = "Conditional merge must use `copy=True`"
    with pytest.raises(ValueError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, copy=False)

    error_msg = "Conditional merge does not support validation"
    with pytest.raises(NotImplementedError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, validate=True)

    error_msg = "Cannot sort on join keys in a conditional merge"
    with pytest.raises(ValueError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, sort=True)


def test_merge_conditional_bad_on_func():
    # too few on func args
    with pytest.raises(TypeError):
        merge(DataFrame(), DataFrame(), on=lambda l: None)

    # too many on func args
    with pytest.raises(TypeError):
        merge(DataFrame(), DataFrame(), on=lambda l, r, x: None)

    # bad on func return value
    with pytest.raises(Exception):
        merge(DataFrame(), DataFrame(), on=lambda l, r: object())
