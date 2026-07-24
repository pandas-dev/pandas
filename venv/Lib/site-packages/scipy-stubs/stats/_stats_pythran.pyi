from typing import Final, Literal

import optype.numpy as onp

from ._stats_mstats_common import SiegelslopesResult

###

type _Method = Literal["hierarchical", "separate"]

###

__pythran__: Final[tuple[str, str]] = ...  # undocumented

def siegelslopes(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, method: _Method = "hierarchical"
) -> SiegelslopesResult: ...  # undocumented
