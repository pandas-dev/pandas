# defined in scipy/linalg/src/_batched_linalg_module.cc

from typing import final

import optype.numpy as onp
import optype.numpy.compat as npc

###

type _InexactND = onp.ArrayND[npc.inexact32 | npc.inexact64]

###

@final
class error(Exception): ...  # undocumented

#
def _inv[ArrayT: _InexactND](
    Am: ArrayT, structure: int, overwrite_a: bool, lower: bool, /
) -> tuple[ArrayT, list[dict[str, float]]]: ...

#
def _solve[ArrayT: _InexactND](
    Am: ArrayT, b: _InexactND, structure: int, overwrite_a: bool, transposed: bool, lower: bool, /
) -> tuple[ArrayT, list[dict[str, float]]]: ...
