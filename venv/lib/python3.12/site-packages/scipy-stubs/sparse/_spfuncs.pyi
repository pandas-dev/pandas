from typing import Literal as L, TypeAlias, TypeVar

import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase

__all__ = ["count_blocks", "estimate_blocksize"]

_SizeT = TypeVar("_SizeT", bound=int)
_BlockSize: TypeAlias = tuple[_SizeT, _SizeT]

def estimate_blocksize(A: _spbase | onp.ToComplex2D, efficiency: float | npc.floating = 0.7) -> _BlockSize[L[1, 2, 3, 4, 6]]: ...

#
def count_blocks(A: _spbase | onp.ToComplex2D, blocksize: tuple[onp.ToJustInt, onp.ToJustInt]) -> int: ...
