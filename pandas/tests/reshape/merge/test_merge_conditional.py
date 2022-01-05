import re

import pytest

from pandas import (
    DataFrame,
    option_context,
    get_option
)
import pandas._testing as tm
from pandas.core.reshape.merge import merge


# choose chunk sizes such that the merges will happen with a single
# chunk or multiple chunks
@pytest.mark.parametrize("chunk_size", [2, 4, 10])
def test_merge_conditional(chunk_size):
    # GH#8962

    with option_context("conditional_merge.chunk_size", chunk_size):
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
def test_non_inner(how):
    error_msg = 'Conditional merge is currently only available for how="inner".'
    with pytest.raises(NotImplementedError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, how=how)


def test_bad_args_clashing_ons():
    error_msg = re.escape(
        "Cannot define any of (`left_on`, `right_on`, `left_index`, "
        "`right_index`) in a conditional merge"
    )
    with pytest.raises(ValueError, match=error_msg):
        merge(
            DataFrame(columns=["A"]),
            DataFrame(columns=["A"]),
            on=lambda l, r: None,
            left_on="A",
        )


def test_bad_args_copy_false():
    error_msg = "Conditional merge must use `copy=True`"
    with pytest.raises(ValueError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, copy=False)


def test_bad_args_validate_true():
    error_msg = "Conditional merge does not support validation"
    with pytest.raises(NotImplementedError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, validate=True)


def test_bad_args_sort_true():
    error_msg = "Cannot sort on join keys in a conditional merge"
    with pytest.raises(ValueError, match=error_msg):
        merge(DataFrame(), DataFrame(), on=lambda l, r: None, sort=True)


def test_on_func_too_few_args():
    with pytest.raises(TypeError, match=".+takes.+argument.+"):
        merge(
            DataFrame([[1, 2]], columns=["A", "B"]),
            DataFrame([[1, 2]], columns=["A", "B"]),
            on=lambda l: None,
        )


def test_on_func_too_many_args():
    with pytest.raises(TypeError, match=".+missing.+argument.+"):
        merge(
            DataFrame([[1, 2]], columns=["A", "B"]),
            DataFrame([[1, 2]], columns=["A", "B"]),
            on=lambda l, r, x: None,
        )


def test_on_func_bad_return_value():
    error_msg = 'Callable `on` condition must return a mask of dtype="bool"'
    with pytest.raises(ValueError, match=error_msg):
        merge(
            DataFrame([[1, 2]], columns=["A", "B"]),
            DataFrame([[1, 2]], columns=["A", "B"]),
            on=lambda l, r: object(),
        )


def test_chunk_size():
    option_name = "conditional_merge.chunk_size"

    with option_context(option_name, 1):
        assert get_option(option_name) == (1, 1)

    with option_context(option_name, (1, 1)):
        assert get_option(option_name) == (1, 1)


@pytest.mark.parametrize(
    "option",
    [
        0
        (1,),
        ('a',),
        (1, 'a'),
        ('a', 'b'),
        (1, 2, 123),
        (-1, -12),
        object(),
        (-1, 1),
        'a',
        -1,
    ]
)
def test_chunk_size_option_bad(option):
    error_msg = (
        "option must be an int greater than 0, or a tuple of 2 ints, "
        "both greater than 0"
    )
    with pytest.raises(ValueError, match=error_msg):
        with option_context("conditional_merge.chunk_size", option):
            pass
