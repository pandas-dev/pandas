from typing import Any, Generic, Literal, Self, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._crosstab import crosstab
from ._odds_ratio import odds_ratio
from ._relative_risk import relative_risk
from ._resampling import ResamplingMethod
from ._typing import BaseBunch, PowerDivergenceStatistic

__all__ = ["association", "chi2_contingency", "crosstab", "expected_freq", "margins", "odds_ratio", "relative_risk"]

###

type _to_floating = npc.floating | npc.integer | np.bool  # noqa: PYI042

_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

###

# NOTE: On numpy<2.1, pyright reports 12 false positive incompatible overload errors here.
# pyright: reportOverlappingOverload=false

class Chi2ContingencyResult(BaseBunch[np.float64, np.float64, int, onp.ArrayND[np.float64]], Generic[_ShapeT_co]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    def pvalue(self, /) -> np.float64: ...
    @property
    def dof(self, /) -> int: ...
    @property
    def expected_freq(self, /) -> onp.ArrayND[np.float64, _ShapeT_co]: ...

    #
    @override
    def __new__(  # pyrefly:ignore[bad-override]
        _cls, statistic: np.float64, pvalue: np.float64, dof: int, expected_freq: onp.ArrayND[np.float64, _ShapeT_co]
    ) -> Self: ...
    @override
    def __init__(  # pyrefly:ignore[bad-override]
        self, /, statistic: np.float64, pvalue: np.float64, dof: int, expected_freq: onp.ArrayND[np.float64, _ShapeT_co]
    ) -> None: ...

#
@overload
def margins[ScalarT: npc.number | np.timedelta64, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[ScalarT, ShapeT],
) -> list[onp.ArrayND[ScalarT, ShapeT]]: ...
@overload
def margins(a: onp.ArrayND[np.bool]) -> list[onp.ArrayND[np.int_]]: ...

#
@overload
def expected_freq[ShapeT: tuple[int, ...]](observed: onp.ArrayND[_to_floating, ShapeT]) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def expected_freq(observed: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

#
@overload
def chi2_contingency[ShapeT: tuple[int, ...]](
    observed: onp.ArrayND[_to_floating, ShapeT],
    correction: bool = True,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    method: ResamplingMethod | None = None,
) -> Chi2ContingencyResult[ShapeT]: ...
@overload
def chi2_contingency(
    observed: onp.ToFloatND,
    correction: bool = True,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    method: ResamplingMethod | None = None,
) -> Chi2ContingencyResult: ...

#
def association(
    observed: onp.ToJustIntND,
    method: Literal["cramer", "tschuprow", "pearson"] = "cramer",
    correction: bool = False,
    lambda_: PowerDivergenceStatistic | float | None = None,
) -> float: ...
