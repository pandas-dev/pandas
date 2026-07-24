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
import sys

from hypothesis.internal.conjecture.floats import float_to_lex
from hypothesis.internal.conjecture.shrinking.common import Shrinker
from hypothesis.internal.conjecture.shrinking.integer import Integer
from hypothesis.internal.floats import MAX_PRECISE_INTEGER, float_to_int, int_to_float

# Bit pattern of the boundary float, so we can compute float-grid indices
# relative to it without recomputing the constant on every call.
_BOUNDARY_BITS = float_to_int(float(MAX_PRECISE_INTEGER))


def _float_to_position(f: float) -> int:
    """Map a non-negative float to a linear integer position such that adjacent
    representable floats correspond to adjacent integers.

    For ``f <= MAX_PRECISE_INTEGER`` the position is just ``int(f)``. Above the
    boundary, where the gap between adjacent floats exceeds 1, we extend by the
    float's index in the bit-pattern sequence past ``MAX_PRECISE_INTEGER``, so
    that decrementing the position by 1 corresponds to ``next_down(f)``.
    """
    if f <= MAX_PRECISE_INTEGER:
        return int(f)
    return MAX_PRECISE_INTEGER + (float_to_int(f) - _BOUNDARY_BITS)


def _position_to_float(n: int) -> float:
    """Inverse of :func:`_float_to_position` on the integer-valued range. Always
    returns an integer-valued, non-negative float."""
    if n <= MAX_PRECISE_INTEGER:
        return float(n)
    return int_to_float(_BOUNDARY_BITS + (n - MAX_PRECISE_INTEGER))


class Float(Shrinker):
    def setup(self):
        self.debugging_enabled = True

    def make_canonical(self, f):
        if math.isnan(f):
            # Distinguish different NaN bit patterns, while making each equal to itself.
            # Wrap in tuple to avoid potential collision with (huge) finite floats.
            return ("nan", float_to_int(f))
        return f

    def check_invariants(self, value):
        # We only handle positive floats (including NaN) because we encode the sign
        # separately anyway.
        assert not (value < 0)

    def left_is_better(self, left, right):
        lex1 = float_to_lex(left)
        lex2 = float_to_lex(right)
        return lex1 < lex2

    def short_circuit(self):
        # We check for a bunch of standard "large" floats. If we're currently
        # worse than them and the shrink downwards doesn't help, abort early
        # because there's not much useful we can do here.

        for g in [sys.float_info.max, math.inf, math.nan]:
            self.consider(g)

        # If we're stuck at a nasty float don't try to shrink it further.
        if not math.isfinite(self.current):
            return True

    def run_step(self):
        # Above MAX_PRECISE_INTEGER all floats are integers, but the gap between
        # adjacent floats is > 1, so consecutive integers are not all
        # representable. Integer.shrink would step by n - 1, which rounds straight
        # back to n and stalls. We instead shrink on the float grid by delegating
        # to Integer with a bijection that maps each representable float to an
        # adjacent integer position, so n - 1 always corresponds to next_down(n).
        if self.current > MAX_PRECISE_INTEGER:
            self.delegate(
                Integer,
                convert_to=_float_to_position,
                convert_from=_position_to_float,
            )
            return

        # Finally we get to the important bit: Each of these is a small change
        # to the floating point number that corresponds to a large change in
        # the lexical representation. Trying these ensures that our floating
        # point shrink can always move past these obstacles. In particular it
        # ensures we can always move to integer boundaries and shrink past a
        # change that would require shifting the exponent while not changing
        # the float value much.

        # First, try dropping precision bits by rounding the scaled value. We
        # try values ordered from least-precise (integer) to more precise, ie.
        # approximate lexicographical order. Once we find an acceptable shrink,
        # self.consider discards the remaining attempts early and skips test
        # invocation. The loop count sets max fractional bits to keep, and is a
        # compromise between completeness and performance.

        for p in range(10):
            scaled = self.current * 2**p  # note: self.current may change in loop
            for truncate in [math.floor, math.ceil]:
                self.consider(truncate(scaled) / 2**p)

        if self.consider(int(self.current)):
            self.debug("Just an integer now")
            self.delegate(Integer, convert_to=int, convert_from=float)
            return

        # Now try to minimize the top part of the fraction as an integer. This
        # basically splits the float as k + x with 0 <= x < 1 and minimizes
        # k as an integer, but without the precision issues that would have.
        m, n = self.current.as_integer_ratio()
        i, r = divmod(m, n)
        self.call_shrinker(Integer, i, lambda k: self.consider((k * n + r) / n))
