from typing import Any, Final, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array
from scipy.sparse._base import _spbase

_RealT = TypeVar("_RealT", bound=npc.integer | npc.floating)
_Graph: TypeAlias = onp.CanArrayND[_RealT] | _spbase[_RealT, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

@overload
def minimum_spanning_tree(csgraph: _Graph[_RealT], overwrite: bool = False) -> csr_array[_RealT, tuple[int, int]]: ...
@overload
def minimum_spanning_tree(csgraph: onp.ToJustInt2D, overwrite: bool = False) -> csr_array[np.int_, tuple[int, int]]: ...
@overload
def minimum_spanning_tree(csgraph: onp.ToJustFloat64_2D, overwrite: bool = False) -> csr_array[np.float64, tuple[int, int]]: ...
@overload
def minimum_spanning_tree(csgraph: onp.ToFloat2D, overwrite: bool = False) -> csr_array[Any, tuple[int, int]]: ...
