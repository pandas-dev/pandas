from dataclasses import dataclass
from typing import Literal

import numpy as np
import optype as op
import optype.numpy as onp

@dataclass
class PageTrendTestResult:
    statistic: np.float64
    pvalue: np.float64
    method: Literal["asymptotic", "exact"]

def page_trend_test(
    data: onp.ToFloatND,
    ranked: op.CanBool = False,
    predicted_ranks: onp.ToIntND | None = None,
    method: Literal["auto", "asymptotic", "exact"] = "auto",
) -> PageTrendTestResult: ...
