# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import math
from collections.abc import Callable
from sys import float_info
from typing import overload

from hypothesis._native.internal.floats import (
    count_between_floats as count_between_floats,
    float_of as float_of,
    float_to_int as float_to_int,
    int_to_float as int_to_float,
    is_negative as is_negative,
    next_down as next_down,
    next_down_normal as next_down_normal,
    next_up as next_up,
    next_up_normal as next_up_normal,
    width_smallest_normals as width_smallest_normals,
)

assert width_smallest_normals(64) == float_info.min

mantissa_mask = (1 << 52) - 1


def make_float_clamper(
    min_value: float,
    max_value: float,
    *,
    allow_nan: bool,
    smallest_nonzero_magnitude: float,
) -> Callable[[float], float]:
    """
    Return a function that clamps positive floats into the given bounds.
    """
    from hypothesis.internal.conjecture.choice import choice_permitted

    assert sign_aware_lte(min_value, max_value)
    range_size = min(max_value - min_value, float_info.max)

    def float_clamper(f: float) -> float:
        if choice_permitted(
            f,
            {
                "min_value": min_value,
                "max_value": max_value,
                "allow_nan": allow_nan,
                "smallest_nonzero_magnitude": smallest_nonzero_magnitude,
            },
        ):
            return f
        # Outside bounds; pick a new value, sampled from the allowed range,
        # using the mantissa bits.
        mant = float_to_int(abs(f)) & mantissa_mask
        f = min_value + range_size * (mant / mantissa_mask)

        # if we resampled into the space disallowed by smallest_nonzero_magnitude,
        # default to smallest_nonzero_magnitude.
        if 0 < abs(f) < smallest_nonzero_magnitude:
            f = smallest_nonzero_magnitude
            # we must have either -smallest_nonzero_magnitude <= min_value or
            # smallest_nonzero_magnitude >= max_value, or no values would be
            # possible. If smallest_nonzero_magnitude is not valid (because it's
            # larger than max_value), then -smallest_nonzero_magnitude must be valid.
            if smallest_nonzero_magnitude > max_value:
                f *= -1

        # Re-enforce the bounds (just in case of floating point arithmetic error)
        return clamp(min_value, f, max_value)

    return float_clamper


def sign_aware_lte(x: float | int, y: float | int) -> bool:
    """Less-than-or-equals, but strictly orders -0.0 and 0.0"""
    if x == 0.0 == y:
        return math.copysign(1.0, x) <= math.copysign(1.0, y)
    else:
        return x <= y


@overload
def clamp(lower: int, value: int, upper: int) -> int: ...
@overload
def clamp(lower: float, value: float, upper: float) -> float: ...
def clamp(lower: float | int, value: float | int, upper: float | int) -> float | int:
    """Given a value and lower/upper bounds, 'clamp' the value so that
    it satisfies lower <= value <= upper.  NaN is mapped to lower."""
    # this seems pointless (and is for integers), but handles the -0.0/0.0 case.
    if not sign_aware_lte(lower, value):
        return lower
    if not sign_aware_lte(value, upper):
        return upper
    return value


SMALLEST_SUBNORMAL = next_up(0.0)
SIGNALING_NAN = int_to_float(0x7FF8_0000_0000_0001)  # nonzero mantissa
MAX_PRECISE_INTEGER = 2**53
assert math.isnan(SIGNALING_NAN)
assert math.copysign(1, SIGNALING_NAN) == 1
