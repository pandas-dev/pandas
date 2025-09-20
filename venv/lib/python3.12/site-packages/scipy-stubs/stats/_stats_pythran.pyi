from typing import Final, Literal, TypeAlias

import optype.numpy as onp

from ._stats_mstats_common import SiegelslopesResult

_Method: TypeAlias = Literal["hierarchical", "separate"]

###

__pythran__: Final[tuple[str, str]] = ...

def siegelslopes(y: onp.ToFloatND, x: onp.ToFloatND | None = None, method: _Method = "hierarchical") -> SiegelslopesResult: ...
