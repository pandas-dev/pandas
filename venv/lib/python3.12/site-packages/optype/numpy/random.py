"""
Type aliases for [SPEC 7](https://scientific-python.org/specs/spec-0007/).
"""

import sys
from collections.abc import Sequence
from typing import Any, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeAliasType
else:
    from typing_extensions import TypeAliasType

import numpy as np

__all__ = ["RNG", "ToRNG", "ToSeed"]

###

_Integral: TypeAlias = np.integer[Any] | np.timedelta64
_IntOrSequence: TypeAlias = int | _Integral | Sequence[int | _Integral]

ToSeed = TypeAliasType(
    "ToSeed",
    np.random.SeedSequence
    | _IntOrSequence
    | np.ndarray[Any, np.dtype[_Integral | np.flexible | np.object_]],
)

RNG = TypeAliasType("RNG", np.random.Generator | np.random.RandomState)
ToRNG = TypeAliasType("ToRNG", RNG | np.random.BitGenerator | ToSeed)
