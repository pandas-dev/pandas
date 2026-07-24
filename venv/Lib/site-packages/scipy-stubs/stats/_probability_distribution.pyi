import abc
import types
from typing import Any, Generic, Literal as L, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._qmc import QMCEngine

_XT_co = TypeVar("_XT_co", bound=npc.number, default=np.float64, covariant=True)

type _MedianMethod = L["formula", "icdf"] | None
type _ModeMethod = L["formula", "optimization"] | None
type _SampleMethod = L["formula", "inverse_transform"] | None
type _RMomentMethod = L["formula", "transform", "quadrature", "cache"] | None
type _CMomentMethod = L["formula", "transform", "quadrature", "cache", "normalize"] | None
type _SMomentMethod = L["formula", "transform", "general", "cache", "normalize"] | None
type _LMomentMethod = L["formula", "general", "order_statistics", "quadrature_icdf", "cache"]
type _EntropyMethod = L["formula", "logexp", "quadrature"] | None
type _PXFMethod = L["formula", "logexp"] | None
type _CDFMethod = L["formula", "logexp", "complement", "quadrature", "subtraction"] | None
type _CCDFMethod = L["formula", "logexp", "complement", "quadrature", "addition"] | None
type _ICDFMethod = L["formula", "complement", "inversion"] | None

type _Float = np.float64 | np.longdouble
type _Complex = np.complex128 | np.clongdouble
type _FloatND = onp.ArrayND[_Float]
type _ToFloat0ND = onp.ToFloat | onp.ToFloatND

###

class _ProbabilityDistribution(Generic[_XT_co], metaclass=abc.ABCMeta):
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    #
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

    # The `kind` param here is wrong: https://github.com/scipy/scipy/pull/25243
    @abc.abstractmethod
    def lmoment(self, /, order: int, *, standardize: bool, method: _LMomentMethod | None) -> _Float | _FloatND: ...

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
