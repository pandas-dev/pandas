# defined in scipy/linalg/src/_batched_linalg_module.cc

from typing import TypeVar, final

import optype.numpy as onp
import optype.numpy.compat as npc

_ArrayT = TypeVar("_ArrayT", bound=onp.ArrayND[npc.inexact32 | npc.inexact64])

###

@final
class error(Exception): ...  # undocumented

def _inv(Am: _ArrayT, structure: int, overwrite_a: bool, lower: bool, /) -> tuple[_ArrayT, list[dict[str, float]]]: ...
def _solve(
    Am: _ArrayT,
    b: onp.ArrayND[npc.inexact32 | npc.inexact64],
    structure: int,
    overwrite_a: bool,
    transposed: bool,
    lower: bool,
    /,
) -> tuple[_ArrayT, list[dict[str, float]]]: ...
