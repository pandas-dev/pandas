from dataclasses import dataclass
from typing import Any, Final, Literal, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

from ._censored_data import CensoredData
from ._common import ConfidenceInterval
from ._typing import Alternative

__all__ = ["ecdf", "logrank"]

_EDFKind: TypeAlias = Literal["cdf", "sf"]
_CIMethod: TypeAlias = Literal["linear", "log-log"]

_Int1D: TypeAlias = onp.Array1D[np.int_]
_Float1D: TypeAlias = onp.Array1D[np.float64]

_KwargsT = TypeVar("_KwargsT")
_KwargsT_contra = TypeVar("_KwargsT_contra", contravariant=True)
_LineT = TypeVar("_LineT")

_SampleData: TypeAlias = onp.ToFloatND | CensoredData

@type_check_only
class _CanStep(Protocol[_KwargsT_contra, _LineT]):
    def step(self, x: _Float1D, y: _Float1D, /, **kwargs: _KwargsT_contra) -> list[_LineT]: ...

###

@dataclass
class EmpiricalDistributionFunction:
    # NOTE: the order of attributes matters
    quantiles: _Float1D
    probabilities: _Float1D
    _n: _Int1D
    _d: _Int1D
    _sf: _Float1D
    _kind: _EDFKind

    def __init__(self, /, q: _Float1D, p: _Float1D, n: _Int1D, d: _Int1D, kind: _EDFKind) -> None: ...
    def evaluate(self, /, x: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...
    @overload
    def plot(self, /, ax: None = None, **kwds: object) -> list[Any]: ...
    @overload
    def plot(self, /, ax: _CanStep[_KwargsT, _LineT], **kwds: _KwargsT) -> list[_LineT]: ...
    def confidence_interval(
        self, /, confidence_level: onp.ToFloat = 0.95, *, method: _CIMethod = "linear"
    ) -> ConfidenceInterval: ...

@dataclass
class ECDFResult:
    cdf: Final[EmpiricalDistributionFunction]
    sf: Final[EmpiricalDistributionFunction]

    def __init__(self, /, q: _Float1D, cdf: _Float1D, sf: _Float1D, n: _Int1D, d: _Int1D) -> None: ...

@dataclass
class LogRankResult:
    statistic: np.float64
    pvalue: np.float64

def ecdf(sample: _SampleData) -> ECDFResult: ...
def logrank(x: _SampleData, y: _SampleData, alternative: Alternative = "two-sided") -> LogRankResult: ...
