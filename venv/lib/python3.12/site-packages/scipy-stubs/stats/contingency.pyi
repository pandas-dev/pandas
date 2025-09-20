from typing import Literal, Self

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._crosstab import crosstab
from ._odds_ratio import odds_ratio
from ._relative_risk import relative_risk
from ._resampling import ResamplingMethod
from ._typing import BaseBunch, PowerDivergenceStatistic

__all__ = ["association", "chi2_contingency", "crosstab", "expected_freq", "margins", "odds_ratio", "relative_risk"]

class Chi2ContingencyResult(BaseBunch[np.float64, np.float64, int, onp.ArrayND[np.float64]]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    def pvalue(self, /) -> np.float64: ...
    @property
    def dof(self, /) -> int: ...
    @property
    def expected_freq(self, /) -> onp.ArrayND[np.float64]: ...
    def __new__(_cls, statistic: np.float64, pvalue: np.float64, dof: int, expected_freq: onp.ArrayND[np.float64]) -> Self: ...
    def __init__(
        self, /, statistic: np.float64, pvalue: np.float64, dof: int, expected_freq: onp.ArrayND[np.float64]
    ) -> None: ...

#
def margins(a: onp.ArrayND[npc.number | np.bool_ | np.timedelta64]) -> list[onp.ArrayND[npc.number | np.timedelta64]]: ...

#
def expected_freq(observed: onp.ToFloatND) -> np.float64 | onp.ArrayND[np.float64]: ...

#
def chi2_contingency(
    observed: onp.ToFloatND,
    correction: bool = True,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    method: ResamplingMethod | None = None,
) -> Chi2ContingencyResult: ...

#
def association(
    observed: onp.ToFloatND,
    method: Literal["cramer", "tschuprow", "pearson"] = "cramer",
    correction: bool = False,
    lambda_: PowerDivergenceStatistic | float | None = None,
) -> float: ...
