from collections.abc import Callable
from typing import NamedTuple, Protocol, Self, overload, type_check_only

import numpy as np
import optype.numpy as onp

import scipy.stats as stats

__all__ = ["DiscreteAliasUrn", "NumericalInversePolynomial", "TransformedDensityRejection", "UNURANError"]

@type_check_only
class _HasSupport(Protocol):
    @property
    def support(self, /) -> tuple[float, float]: ...

@type_check_only
class _HasPMF(_HasSupport, Protocol):
    @property
    def pmf(self, /) -> Callable[..., float]: ...

@type_check_only
class _HasPDF(_HasSupport, Protocol):
    @property
    def pdf(self, /) -> Callable[..., float]: ...

@type_check_only
class _HasCDF(_HasPDF, Protocol):
    @property
    def cdf(self, /) -> Callable[..., float]: ...

@type_check_only
class _TDRDist(_HasPDF, Protocol):
    @property
    def dpdf(self, /) -> Callable[..., float]: ...

@type_check_only
class _PINVDist(_HasCDF, Protocol):
    @property
    def logpdf(self, /) -> Callable[..., float]: ...

@type_check_only
class _PPFMethodMixin:
    @overload
    def ppf(self, /, u: onp.ToFloat) -> float: ...
    @overload
    def ppf(self, /, u: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

###

class UNURANError(RuntimeError): ...

class UError(NamedTuple):
    max_error: float
    mean_absolute_error: float

class Method:
    @overload
    def rvs(self, /, size: None = None, random_state: onp.random.ToRNG | None = None) -> float | int: ...
    @overload
    def rvs(self, /, size: int | tuple[int, ...]) -> onp.ArrayND[np.float64 | np.int_]: ...
    def set_random_state(self, /, random_state: onp.random.ToRNG | None = None) -> None: ...

class TransformedDensityRejection(Method):
    def __new__(
        cls,
        dist: _TDRDist,
        *,
        mode: float | None = ...,
        center: float | None = ...,
        domain: tuple[float, float] | None = ...,
        c: float = ...,
        construction_points: onp.ToFloatND = ...,
        use_dars: bool = ...,
        max_squeeze_hat_ratio: float = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...
    @property
    def hat_area(self, /) -> float: ...
    @property
    def squeeze_hat_ratio(self, /) -> float: ...
    @property
    def squeeze_area(self, /) -> float: ...
    @overload
    def ppf_hat(self, /, u: onp.ToFloat) -> float: ...
    @overload
    def ppf_hat(self, /, u: onp.ToScalar | onp.ToArrayND) -> float | onp.ArrayND[np.float64]: ...

class SimpleRatioUniforms(Method):
    def __new__(
        cls,
        dist: _HasPDF,
        *,
        mode: float | None = ...,
        pdf_area: float = ...,
        domain: tuple[float, float] | None = ...,
        cdf_at_mode: float = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...

class NumericalInversePolynomial(_PPFMethodMixin, Method):
    def __new__(
        cls,
        dist: _PINVDist,
        *,
        mode: float | None = ...,
        center: float | None = ...,
        domain: tuple[float, float] | None = ...,
        order: int = ...,
        u_resolution: float = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...
    @property
    def intervals(self, /) -> int: ...
    @overload
    def cdf(self, /, x: onp.ToFloat) -> float: ...
    @overload
    def cdf(self, /, x: onp.ToFloat | onp.ToFloatND) -> float | onp.ArrayND[np.float64]: ...
    def u_error(self, /, sample_size: int = 100_000) -> UError: ...
    def qrvs(
        self, /, size: int | tuple[int, ...] | None = None, d: int | None = None, qmc_engine: stats.qmc.QMCEngine | None = None
    ) -> float | onp.ArrayND[np.float64]: ...

class NumericalInverseHermite(_PPFMethodMixin, Method):
    def __new__(
        cls,
        dist: _HasCDF,
        *,
        domain: tuple[float, float] | None = ...,
        order: int = ...,
        u_resolution: float = ...,
        construction_points: onp.ToFloatND | None = ...,
        max_intervals: int = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...
    @property
    def intervals(self, /) -> int: ...
    @property
    def midpoint_error(self, /) -> float: ...
    def u_error(self, /, sample_size: int = 100_000) -> UError: ...
    def qrvs(
        self, /, size: int | tuple[int, ...] | None = None, d: int | None = None, qmc_engine: stats.qmc.QMCEngine | None = None
    ) -> float | onp.ArrayND[np.float64]: ...

class DiscreteAliasUrn(Method):
    def __new__(
        cls,
        dist: onp.ToFloat | onp.ToFloatND | _HasPMF,
        *,
        domain: tuple[float, float] | None = ...,
        urn_factor: float = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...

class DiscreteGuideTable(_PPFMethodMixin, Method):
    def __new__(
        cls,
        dist: onp.ToFloat | onp.ToFloatND | _HasPMF,
        *,
        domain: tuple[float, float] | None = ...,
        guide_factor: float = ...,
        random_state: onp.random.ToRNG | None = ...,
    ) -> Self: ...
