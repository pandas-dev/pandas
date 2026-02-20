from collections.abc import Generator, Iterable
from types import ModuleType
from typing import Generic
from typing_extensions import TypeVar

import numpy as np

from ._base import NestedFixedRule

_T = TypeVar("_T")
_XPT_co = TypeVar("_XPT_co", default=ModuleType, covariant=True)

###

class GenzMalikCubature(NestedFixedRule[_XPT_co, np.float64], Generic[_XPT_co]):  # undocumented
    ndim: int
    degree: int
    lower_degree: int
    def __init__(self, /, ndim: int, degree: int = 7, lower_degree: int = 5, xp: _XPT_co | None = None) -> None: ...

def _distinct_permutations(iterable: Iterable[_T]) -> Generator[tuple[_T, ...]]: ...  # undocumented
