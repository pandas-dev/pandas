from typing import Literal, TypeAlias

from ._common import ConfidenceInterval

_Alternative: TypeAlias = Literal["two-sided", "less", "greater"]

class BinomTestResult:
    k: int
    n: int
    alternative: _Alternative
    statistic: float
    pvalue: float

    def __init__(self, /, k: int, n: int, alternative: _Alternative, statistic: float, pvalue: float) -> None: ...
    def proportion_ci(
        self, /, confidence_level: float = 0.95, method: Literal["exact", "wilson", "wilsoncc"] = "exact"
    ) -> ConfidenceInterval: ...

def binomtest(k: int, n: int, p: float = 0.5, alternative: _Alternative = "two-sided") -> BinomTestResult: ...
