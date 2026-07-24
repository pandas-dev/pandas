"""
Type aliases for [SPEC 7](https://scientific-python.org/specs/spec-0007/).
"""

from collections.abc import Sequence
from typing import Any

import numpy as np

__all__ = ["RNG", "ToRNG", "ToSeed"]

###

type _Integral = np.integer[Any] | np.timedelta64
type _IntOrSequence = int | _Integral | Sequence[int | _Integral]

type ToSeed = (
    np.random.SeedSequence
    | _IntOrSequence
    | np.ndarray[Any, np.dtype[_Integral | np.flexible | np.object_]]
)

type RNG = np.random.Generator | np.random.RandomState
type ToRNG = RNG | np.random.BitGenerator | ToSeed
