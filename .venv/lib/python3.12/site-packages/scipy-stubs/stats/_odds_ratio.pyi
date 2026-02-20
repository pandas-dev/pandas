from typing import Literal, TypeAlias

import optype.numpy as onp
import optype.numpy.compat as npc

from ._common import ConfidenceInterval

_Kind: TypeAlias = Literal["conditional", "sample"]

class OddsRatioResult:
    statistic: float
    def __init__(self, /, _table: onp.Array2D[npc.integer], _kind: _Kind, statistic: float) -> None: ...
    def confidence_interval(self, /, confidence_level: float = 0.95, alternative: str = "two-sided") -> ConfidenceInterval: ...

def odds_ratio(table: onp.ToInt2D, *, kind: _Kind = "conditional") -> OddsRatioResult: ...
