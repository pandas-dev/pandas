import abc
from typing import Any, Generic, Literal as L, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._qmc import QMCEngine

_XT_co = TypeVar("_XT_co", bound=npc.number, default=np.float64, covariant=True)

_MedianMethod: TypeAlias = L["formula", "icdf"] | None
_ModeMethod: TypeAlias = L["formula", "optimization"] | None
_SampleMethod: TypeAlias = L["formula", "inverse_transform"] | None
_RMomentMethod: TypeAlias = L["formula", "transform", "quadrature", "cache"] | None
_CMomentMethod: TypeAlias = L["formula", "transform", "quadrature", "cache", "normalize"] | None
_SMomentMethod: TypeAlias = L["formula", "transform", "general", "cache", "normalize"] | None
_EntropyMethod: TypeAlias = L["formula", "logexp", "quadrature"] | None
_PXFMethod: TypeAlias = L["formula", "logexp"] | None
_CDFMethod: TypeAlias = L["formula", "logexp", "complement", "quadrature", "subtraction"] | None
_CCDFMethod: TypeAlias = L["formula", "logexp", "complement", "quadrature", "addition"] | None
_ICDFMethod: TypeAlias = L["formula", "complement", "inversion"] | None

_Float: TypeAlias = np.float64 | np.longdouble
_Complex: TypeAlias = np.complex128 | np.clongdouble
_FloatND: TypeAlias = onp.ArrayND[_Float]
_ToFloat0ND: TypeAlias = onp.ToFloat | onp.ToFloatND

###

class _ProbabilityDistribution(Generic[_XT_co], metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def support(self, /) -> tuple[_XT_co | onp.ArrayND[_XT_co], _XT_co | onp.ArrayND[_XT_co]]: ...
    @abc.abstractmethod
    def median(self, /, *, method: _MedianMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def mode(self, /, *, method: _ModeMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def sample(
        self, /, shape: int | tuple[int, ...], *, method: _SampleMethod, rng: QMCEngine | onp.random.ToRNG | None
    ) -> _XT_co | onp.Array[Any, _XT_co]: ...  # `Any` shape is needed on `numpy<2.1`

    #
    @abc.abstractmethod
    def mean(self, /, *, method: _RMomentMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def variance(self, /, *, method: _CMomentMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def standard_deviation(self, /, *, method: _CMomentMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def skewness(self, /, *, method: _SMomentMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def kurtosis(self, /, *, method: _SMomentMethod) -> _XT_co | onp.ArrayND[_XT_co]: ...

    #
    @overload
    @abc.abstractmethod
    def moment(self, /, order: onp.ToInt, kind: L["raw"], *, method: _RMomentMethod) -> _Float | _FloatND: ...
    @overload
    @abc.abstractmethod
    def moment(self, /, order: onp.ToInt, kind: L["central"], *, method: _CMomentMethod) -> _Float | _FloatND: ...
    @overload
    @abc.abstractmethod
    def moment(self, /, order: onp.ToInt, kind: L["standardized"], *, method: _SMomentMethod) -> _Float | _FloatND: ...

    #
    @abc.abstractmethod
    def entropy(self, /, *, method: _EntropyMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def logentropy(self, /, *, method: _EntropyMethod) -> _Complex | onp.ArrayND[_Complex]: ...

    #
    @abc.abstractmethod
    def pdf(self, x: _ToFloat0ND, /, *, method: _PXFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def logpdf(self, x: _ToFloat0ND, /, *, method: _PXFMethod) -> _Float | _FloatND: ...

    #
    def pmf(self, x: _ToFloat0ND, /, *, method: _PXFMethod = None) -> _Float | _FloatND: ...
    def logpmf(self, x: _ToFloat0ND, /, *, method: _PXFMethod = None) -> _Float | _FloatND: ...

    #
    @abc.abstractmethod
    def cdf(self, x: _ToFloat0ND, y: _ToFloat0ND | None, /, *, method: _CDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def icdf(self, p: _ToFloat0ND, /, *, method: _ICDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def ccdf(self, x: _ToFloat0ND, y: _ToFloat0ND | None, /, *, method: _CCDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def iccdf(self, p: _ToFloat0ND, /, *, method: _ICDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def logcdf(self, x: _ToFloat0ND, y: _ToFloat0ND | None, /, *, method: _CDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def ilogcdf(self, logp: _ToFloat0ND, /, *, method: _ICDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def logccdf(self, x: _ToFloat0ND, y: _ToFloat0ND | None, /, *, method: _CCDFMethod) -> _Float | _FloatND: ...
    @abc.abstractmethod
    def ilogccdf(self, logp: _ToFloat0ND, /, *, method: _ICDFMethod) -> _Float | _FloatND: ...
