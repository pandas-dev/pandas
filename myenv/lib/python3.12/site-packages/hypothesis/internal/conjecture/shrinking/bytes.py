# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

from hypothesis.internal.compat import int_from_bytes, int_to_bytes
from hypothesis.internal.conjecture.shrinking.integer import Integer


class Bytes(Integer):
    def __init__(self, initial, predicate, **kwargs):
        # shrink by interpreting the bytes as an integer.
        # move to Collection.shrink when we support variable-size bytes,
        # because b'\x00\x02' could shrink to either b'\x00\x01' or b'\x02'.
        super().__init__(
            int_from_bytes(initial),
            lambda n: predicate(int_to_bytes(n, len(initial))),
            **kwargs,
        )
