# TODO(@jorenham): remove `overload-overlap` after dropping support for numpy 2.0
# mypy: disable-error-code="override, explicit-override, overload-overlap"

import abc
import types
from collections.abc import Callable, Iterable, Mapping, Sequence, Set as AbstractSet
from typing import (
    Any,
    ClassVar,
    Final,
    Generic,
    Literal as L,
    Never,
    Protocol,
    Self,
    SupportsIndex,
    TypedDict,
    Unpack,
    final,
    overload,
    override,
    type_check_only,
)
from typing_extensions import ParamSpec, TypeIs, TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._distn_infrastructure import rv_continuous, rv_discrete
from ._probability_distribution import _LMomentMethod, _ProbabilityDistribution
from ._qmc import QMCEngine

__all__ = ["Mixture", "abs", "exp", "log", "make_distribution", "order_statistic", "truncate"]

###

_Tss = ParamSpec("_Tss", default=...)

_FloatT = TypeVar("_FloatT", bound=npc.floating, default=np.float64)
_FloatT_co = TypeVar("_FloatT_co", bound=npc.floating, default=np.float64, covariant=True)
_IntT_co = TypeVar("_IntT_co", bound=npc.integer, default=npc.integer, covariant=True)
_RealT_co = TypeVar("_RealT_co", bound=_Real, default=np.float64 | npc.integer, covariant=True)

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[Any, ...])
_ShapeT1 = TypeVar("_ShapeT1", bound=tuple[int, *tuple[int, ...]], default=tuple[Any, ...])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

_DistT0 = TypeVar("_DistT0", bound=_Dist[_0D])
_DistT_co = TypeVar("_DistT_co", bound=_Dist[tuple[int, ...]], default=UnivariateDistribution, covariant=True)

_AxesT = TypeVar("_AxesT", bound=_Axes, default=Any)

type _ParameterEndpoint = onp.ToFloat | Callable[..., onp.ToFloat] | str
type _ParameterTuple = tuple[_ParameterEndpoint, _ParameterEndpoint]

# NOTE: `TypedDict` with `NotRequired` cannot be used, because `NotRequired` does not
# mean "not required": Other `TypedDict` MUST also include them to be assignable,
# so a `NotRequired` value is REQUIRED. Absolutely ridiculous...
# https://typing.python.org/en/latest/spec/typeddict.html#id4
type _ParameterDict = Mapping[str, tuple[Any, ...]]
type _ParameterSpec = _ParameterDict | _ParameterTuple

@type_check_only
class _DuckDistributionBase(Protocol[_Tss]):
    # NOTE: The dummy `ParamSpec` ensures that `pdf` implementations with required parameter kwargs are assignable.
    @property
    def __make_distribution_version__(self, /) -> str: ...
    @property
    def support(self, /) -> _ParameterSpec: ...
    def pdf(self, x: float, /, *do_not_use_these: _Tss.args, **parameters: _Tss.kwargs) -> float | np.float64: ...

@final
@type_check_only
class _DuckDistributionSingle(_DuckDistributionBase, Protocol):
    @property
    def parameters(self, /) -> Mapping[str, _ParameterSpec]: ...

@final
@type_check_only
class _DuckDistributionMulti(_DuckDistributionBase, Protocol):
    @property
    def parameters(self, /) -> tuple[Mapping[str, _ParameterSpec], ...]: ...
    def process_parameters(self, /) -> Mapping[str, onp.ToFloat]: ...

type _DuckDistributionType = type[_DuckDistributionSingle | _DuckDistributionMulti]

###

type _Real = npc.floating | npc.integer

type _0D = tuple[()]  # noqa: PYI042
type _1D = tuple[int]  # noqa: PYI042
type _2D = tuple[int, int]  # noqa: PYI042
type _3D = tuple[int, int, int]  # noqa: PYI042
type _ND = tuple[int, ...]

type _ToFloatMax1D = onp.ToFloatStrict1D | onp.ToFloat
type _ToFloatMax2D = onp.ToFloatStrict2D | _ToFloatMax1D
type _ToFloatMax3D = onp.ToFloatStrict3D | _ToFloatMax2D

type _ToJustIntMax1D = onp.ToJustIntStrict1D | onp.ToJustInt
type _ToJustIntMax2D = onp.ToJustIntStrict2D | _ToJustIntMax1D
type _ToJustIntMax3D = onp.ToJustIntStrict3D | _ToJustIntMax2D
type _ToJustIntMaxND = onp.ToJustIntND | onp.ToJustInt

type _Null = op.JustObject  # type of `_null`
type _Axes = object  # placeholder for `matplotlib.axes.Axes`

type _DomainRegion = L["domain", "typical"]
type _DomainDrawType = L["in", "out", "on", "nan"]
type _ValidationPolicy = L["skip_all"] | None
type _CachePolicy = L["no_cache"] | None
type _PlotQuantity = L["x", "pdf", "cdf", "ccdf", "icdf", "iccdf", "logpdf", "logcdf", "logccdf", "ilogcdf", "ilogccdf"]
type _SMomentMethod = L["formula", "general", "transform", "normalize", "cache"]
type _CMomentMethod = L["formula", "transform", "quadrature", "cache", "normalize"]
type _RMomentMethod = L["formula", "transform", "quadrature", "cache"]

type _ParamValues = Mapping[str, _ToFloat0ND]
type _ToDomain = tuple[onp.ToFloat | str, onp.ToFloat | str]
type _ToTol = op.JustFloat | _Null
type _DrawProportions = tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat, onp.ToFloat]
type _Elementwise[FloatT: npc.floating] = Callable[[onp.ArrayND[np.float64]], onp.ArrayND[FloatT]]

type _Dist[ShapeT: tuple[int, ...]] = UnivariateDistribution[Any, ShapeT]
type _CDist[ShapeT: tuple[int, ...]] = ContinuousDistribution[np.float64, ShapeT]
type _CDist0[FloatT: npc.floating] = ContinuousDistribution[FloatT, _0D]
type _TransDist[DistT: _Dist[tuple[int, ...]], FloatT: npc.floating, ShapeT: tuple[int, ...]] = TransformedDistribution[
    DistT, FloatT, ShapeT
]
type _LinDist[DistT: _Dist[tuple[int, ...]], FloatT: npc.floating, ShapeT: tuple[int, ...]] = ShiftedScaledDistribution[
    DistT, FloatT, ShapeT
]
type _FoldDist[DistT: _Dist[tuple[int, ...]], FloatT: npc.floating, ShapeT: tuple[int, ...]] = FoldedDistribution[
    DistT, FloatT, ShapeT
]
type _TruncDist[DistT: _Dist[tuple[int, ...]], ShapeT: tuple[int, ...]] = TruncatedDistribution[DistT, ShapeT]

@type_check_only
class _ParamField(Protocol[_FloatT_co, _ShapeT_co]):
    # This actually works (even on mypy)!
    @overload
    def __get__(self: _ParamField[_FloatT, _0D], obj: object, tp: type | None = None, /) -> _FloatT: ...
    @overload
    def __get__(self: _ParamField[_FloatT, _ShapeT1], obj: object, tp: type | None = None, /) -> onp.Array[_ShapeT1, _FloatT]: ...

@type_check_only
class _DistOpts(TypedDict, total=False):
    tol: _ToTol
    validation_policy: _ValidationPolicy
    cache_policy: _CachePolicy

###

_null: Final[_Null] = ...

def _isnull(x: object) -> TypeIs[_Null | None]: ...

class _Domain(abc.ABC, Generic[_XT_co]):  # noqa: UP046
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    # NOTE: This is a `ClassVar[dict[str, float]]` that's overridden as instance attribute in `_SimpleDomain`.
    # https://github.com/scipy/scipy/pull/22139
    symbols: Mapping[str, str] = ...

    @override
    @abc.abstractmethod
    def __str__(self, /) -> str: ...
    @abc.abstractmethod
    def contains(self, /, x: onp.ArrayND[Any]) -> onp.ArrayND[np.bool]: ...
    @abc.abstractmethod
    def draw(self, /, n: int) -> onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def get_numerical_endpoints(self, /, x: _ParamValues) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...

class _Interval(_Domain[_XT_co], Generic[_XT_co]):  # pyrefly: ignore[implicit-abstract-class]  # noqa: UP046
    @override
    @abc.abstractmethod
    def __str__(self, /) -> str: ...

    #
    def __init__(self, /, endpoints: _ToDomain = ..., inclusive: tuple[bool, bool] = (False, False)) -> None: ...
    @override
    def get_numerical_endpoints(  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-param-name-override] # ty: ignore[invalid-method-override]
        self, /, parameter_values: _ParamValues
    ) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
    @override
    def contains(  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
        self, /, item: onp.ArrayND[_Real], parameter_values: _ParamValues | None = None
    ) -> onp.ArrayND[np.bool]: ...
    @override
    def draw(  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
        self,
        /,
        n: int,
        type_: _DomainDrawType,
        min: onp.ArrayND[_Real],
        max: onp.ArrayND[_Real],
        squeezed_base_shape: _ND,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[_XT_co]: ...
    def define_parameters(self, /, *parameters: _Parameter) -> None: ...

class _RealInterval(_Interval[_FloatT_co], Generic[_FloatT_co]):
    @override  # https://github.com/astral-sh/ruff/issues/18372
    def __str__(self, /) -> str: ...  # noqa: PYI029

class _IntegerInterval(_Interval[_IntT_co], Generic[_IntT_co]):
    @override  # https://github.com/astral-sh/ruff/issues/18372
    def __str__(self, /) -> str: ...  # noqa: PYI029

type _ValidateOut0D[RealT: _Real] = tuple[RealT, np.dtype[RealT], onp.Array0D[np.bool]]
type _ValidateOutND[RealT: _Real, _ShapeT1: tuple[int, *tuple[int, ...]]] = tuple[
    onp.ArrayND[RealT, _ShapeT1],
    np.dtype[RealT],
    onp.ArrayND[np.bool, _ShapeT1],
]  # fmt: skip

#
class _Parameter(abc.ABC, Generic[_RealT_co]):
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...
    #
    def __init__(
        self, /, name: str, *, domain: _Domain, symbol: str | None = None, typical: _Domain | _ToDomain | None = None
    ) -> None: ...
    #
    @overload
    @abc.abstractmethod
    def validate(self, /, arr: onp.ToFloat) -> _ValidateOut0D[_RealT_co]: ...
    @overload
    @abc.abstractmethod
    def validate(self, /, arr: onp.ToFloatND) -> _ValidateOutND[_RealT_co, tuple[Any, ...]]: ...
    #
    def draw(
        self,
        /,
        size: _ND | None = None,
        *,
        rng: onp.random.ToRNG | None = None,
        region: _DomainRegion = "domain",
        proportions: _DrawProportions | None = None,
        parameter_values: _ParamValues | None = None,
    ) -> onp.ArrayND[_RealT_co]: ...

class _RealParameter(_Parameter[_FloatT_co], Generic[_FloatT_co]):
    @override
    @overload
    # pyrefly: ignore[bad-override]
    def validate(self, /, arr: onp.ToFloat, parameter_values: _ParamValues) -> _ValidateOut0D[_FloatT_co]: ...
    @overload
    def validate(self, /, arr: onp.ToFloatND, parameter_values: _ParamValues) -> _ValidateOutND[_FloatT_co, tuple[Any, ...]]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

class _Parameterization:
    parameters: Final[Mapping[str, _Parameter]]

    def __init__(self, /, *parameters: _Parameter) -> None: ...
    def __len__(self, /) -> int: ...
    def copy(self, /) -> Self: ...
    def matches(self, /, parameters: AbstractSet[str]) -> bool: ...
    def validation(self, /, parameter_values: Mapping[str, _Parameter]) -> tuple[onp.ArrayND[np.bool], np.dtype[np.float64]]: ...
    def draw(
        self,
        /,
        sizes: _ND | Sequence[_ND] | None = None,
        rng: onp.random.ToRNG | None = None,
        proportions: _DrawProportions | None = None,
        region: _DomainRegion = "domain",
    ) -> dict[str, onp.ArrayND[np.float64]]: ...

###

type _Tuple2[T] = tuple[T, T]

_XT = TypeVar("_XT", bound=npc.number, default=npc.number)
_XT_co = TypeVar("_XT_co", bound=npc.number, default=np.float64, covariant=True)
_ShapeT0_co = TypeVar("_ShapeT0_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

type _BaseDist0[XT: npc.number] = _BaseDistribution[XT, tuple[()]]
type _BaseDist1[XT: npc.number] = _BaseDistribution[XT, tuple[int]]
type _BaseDist2[XT: npc.number] = _BaseDistribution[XT, tuple[int, int]]
type _BaseDist3[XT: npc.number] = _BaseDistribution[XT, tuple[int, int, int]]
type _BaseDist1N[XT: npc.number] = _BaseDistribution[XT, tuple[int, *tuple[int, ...]]]

type _KurtosisConvention = L["non-excess", "excess"]
type _MedianMethod = L["formula", "icdf"]
type _ModeMethod = L["formula", "optimization"]
type _SampleMethod = L["formula", "inverse_transform"]
type _EntropyMethod = L["formula", "logexp", "quadrature"]

type _PXFMethod = L["formula", "logexp"]
type _CDFMethod = L["formula", "logexp", "complement", "quadrature", "subtraction"]
type _CCDFMethod = L["formula", "logexp", "complement", "quadrature", "addition"]
type _ICDFMethod = L["formula", "complement", "inversion"]

type _Float1D = onp.Array1D[np.float64]
type _Float2D = onp.Array2D[np.float64]
type _Float3D = onp.Array3D[np.float64]
type _FloatND[ShapeT: tuple[int, *tuple[int, ...]]] = onp.ArrayND[np.float64, ShapeT]

type _Float1ND = onp.Array[tuple[int, *tuple[Any, ...]], np.float64]
type _Float2ND = onp.Array[tuple[int, int, *tuple[Any, ...]], np.float64]
type _Float3ND = onp.Array[tuple[int, int, int, *tuple[Any, ...]], np.float64]

type _Complex = np.complex128 | np.clongdouble
type _ComplexND[ShapeT: tuple[int, *tuple[int, ...]]] = onp.ArrayND[_Complex, ShapeT]

type _ToFloatND[ShapeT: tuple[int, *tuple[int, ...]]] = onp.CanArrayND[_Real | np.bool, ShapeT]
type _ToFloat0ND = onp.ToFloat | onp.ToFloatND
type _ToFloatMaxND[ShapeT: tuple[int, *tuple[int, ...]]] = _ToFloatND[ShapeT] | _ToFloatMax1D

type _ToQRNG = QMCEngine | onp.random.ToRNG | None

@type_check_only
class _BaseDistribution(_ProbabilityDistribution[_XT_co], Generic[_XT_co, _ShapeT0_co]):
    @override
    @overload
    def support(self: _BaseDist0[_XT], /) -> _Tuple2[_XT]: ...
    @overload
    def support(self: _BaseDistribution[_XT, _ShapeT1], /) -> _Tuple2[onp.Array[_ShapeT1, _XT]]: ...  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]

    #
    @override
    @overload
    def median(self: _BaseDist0[_XT], /, *, method: _MedianMethod | None = None) -> _XT: ...
    @overload
    def median(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _MedianMethod | None = None
    ) -> onp.Array[_ShapeT1, _XT]: ...

    #
    @override
    @overload
    def mode(self: _BaseDist0[_XT], /, *, method: _ModeMethod | None = None) -> _XT: ...
    @overload
    def mode(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _ModeMethod | None = None
    ) -> onp.Array[_ShapeT1, _XT]: ...

    #
    @override
    @overload
    def sample(
        self: _BaseDist0[_XT], /, shape: tuple[()] = (), *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> _XT: ...
    @overload
    def sample(
        self: _BaseDist0[_XT], /, shape: SupportsIndex, *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> onp.Array1D[_XT]: ...
    @overload
    def sample(
        self: _BaseDist0[_XT], /, shape: _ShapeT1, *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...
    @overload
    def sample(
        self: _BaseDistribution[_XT, _ShapeT1],
        /,
        shape: tuple[()] = (),
        *,
        method: _SampleMethod | None = None,
        rng: _ToQRNG = None,
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...
    @overload
    def sample(
        self: _BaseDistribution[_XT, _ShapeT1],
        /,
        shape: SupportsIndex | Iterable[SupportsIndex],
        *,
        method: _SampleMethod | None = None,
        rng: _ToQRNG = None,
    ) -> onp.ArrayND[_XT, _ShapeT1] | onp.ArrayND[_XT]: ...  # first union type is needed on `numpy<2.1`
    @overload
    def sample(
        self, /, shape: SupportsIndex | Iterable[SupportsIndex], *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> _XT_co | onp.ArrayND[_XT_co, _ShapeT1] | onp.ArrayND[_XT_co]: ...  # first union type is needed on `numpy<2.1`

    #
    @override
    @overload
    def mean(self: _BaseDist0[_XT], /, *, method: _RMomentMethod | None = None) -> _XT: ...
    @overload
    def mean(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _RMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @override
    @overload
    def variance(self: _BaseDist0[_XT], /, *, method: _CMomentMethod | None = None) -> _XT: ...
    @overload
    def variance(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _CMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @override
    @overload
    def standard_deviation(self: _BaseDist0[_XT], /, *, method: _CMomentMethod | None = None) -> _XT: ...
    @overload
    def standard_deviation(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _CMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @override
    @overload
    def skewness(self: _BaseDist0[_XT], /, *, method: _SMomentMethod | None = None) -> _XT: ...
    @overload
    def skewness(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _SMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @override
    @overload
    def kurtosis(
        self: _BaseDist0[_XT], /, *, method: _SMomentMethod | None = None, convention: _KurtosisConvention = "non-excess"
    ) -> _XT: ...
    @overload
    def kurtosis(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[_XT, _ShapeT1],
        /,
        *,
        method: _SMomentMethod | None = None,
        convention: _KurtosisConvention = "non-excess",
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @override
    @overload
    def moment(
        self: _BaseDist0[npc.number], /, order: onp.ToInt = 1, kind: L["raw"] = "raw", *, method: _RMomentMethod | None = None
    ) -> np.float64: ...
    @overload
    def moment(
        self: _BaseDist0[npc.number], /, order: onp.ToInt, kind: L["central"], *, method: _CMomentMethod | None = None
    ) -> np.float64: ...
    @overload
    def moment(
        self: _BaseDist0[npc.number], /, order: onp.ToInt = 1, *, kind: L["central"], method: _CMomentMethod | None = None
    ) -> np.float64: ...
    @overload
    def moment(
        self: _BaseDist0[npc.number], /, order: onp.ToInt, kind: L["standardized"], *, method: _SMomentMethod | None = None
    ) -> np.float64: ...
    @overload
    def moment(
        self: _BaseDist0[npc.number], /, order: onp.ToInt = 1, *, kind: L["standardized"], method: _SMomentMethod | None = None
    ) -> np.float64: ...
    @overload
    def moment(
        self: _BaseDistribution[Any, _ShapeT1],
        /,
        order: onp.ToInt = 1,
        kind: L["raw"] = "raw",
        *,
        method: _RMomentMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload
    def moment(
        self: _BaseDistribution[Any, _ShapeT1], /, order: onp.ToInt, kind: L["central"], *, method: _CMomentMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload
    def moment(
        self: _BaseDistribution[Any, _ShapeT1],
        /,
        order: onp.ToInt = 1,
        *,
        kind: L["central"],
        method: _CMomentMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload
    def moment(
        self: _BaseDistribution[Any, _ShapeT1],
        /,
        order: onp.ToInt,
        kind: L["standardized"],
        *,
        method: _SMomentMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload
    def moment(  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        self: _BaseDistribution[Any, _ShapeT1],
        /,
        order: onp.ToInt = 1,
        *,
        kind: L["standardized"],
        method: _SMomentMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...

    #
    @override
    def lmoment(self, /, order: int = 1, *, standardize: bool = True, method: _LMomentMethod | None = None) -> np.float64: ...

    #
    @override
    @overload
    def entropy(self: _BaseDist0[npc.number], /, *, method: _EntropyMethod | None = None) -> np.float64: ...
    @overload
    def entropy(self: _BaseDistribution[Any, _ShapeT1], /, *, method: _EntropyMethod | None = None) -> _FloatND[_ShapeT1]: ...  # ty:ignore[invalid-method-override]

    #
    @override
    @overload
    def logentropy(self: _BaseDist0[npc.number], /, *, method: _EntropyMethod | None = None) -> _Complex: ...
    @overload
    def logentropy(  # ty:ignore[invalid-method-override]
        self: _BaseDistribution[Any, _ShapeT1], /, *, method: _EntropyMethod | None = None
    ) -> _ComplexND[_ShapeT1]: ...

    # NOTE: The signatures of `pdf`, `logpdf`, `pmf`, and `logpmf` are equivalent
    @override
    @overload  # self: T1-d, x: 0-d
    def pdf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def pdf(self: _BaseDist0[npc.number], x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d
    def pdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def pdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def pdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def pdf(
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def pdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def pdf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def pdf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=-d
    def pdf(self: _BaseDist1[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def pdf(self: _BaseDist2[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def pdf(self: _BaseDist2[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def pdf(self: _BaseDist3[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def pdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d
    def logpdf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def logpdf(self: _BaseDist0[npc.number], x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d
    def logpdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def logpdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def logpdf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def logpdf(
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def logpdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def logpdf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def logpdf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d
    def logpdf(self: _BaseDist1[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def logpdf(self: _BaseDist2[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def logpdf(self: _BaseDist2[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def logpdf(self: _BaseDist3[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def logpdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d
    def pmf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def pmf(self: _BaseDist0[npc.number], x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d
    def pmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def pmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def pmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def pmf(
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def pmf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def pmf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def pmf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=-d
    def pmf(self: _BaseDist1[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def pmf(self: _BaseDist2[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def pmf(self: _BaseDist2[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def pmf(self: _BaseDist3[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def pmf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d
    def logpmf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def logpmf(self: _BaseDist0[npc.number], x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d
    def logpmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def logpmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def logpmf(self: _BaseDist0[npc.number], x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def logpmf(
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def logpmf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def logpmf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def logpmf(self: _BaseDist1[npc.number], x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d
    def logpmf(self: _BaseDist1[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def logpmf(self: _BaseDist2[npc.number], x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def logpmf(self: _BaseDist2[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def logpmf(self: _BaseDist3[npc.number], x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def logpmf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    # NOTE: Apart from the `method` type, the signatures of `[log]cdf` and `[log]ccdf` are equivalent
    @override
    @overload  # self: T1-d, x: 0-d, y?: 0-d
    def cdf(
        self: _BaseDistribution[Any, _ShapeT1],
        x: onp.ToFloat,
        y: onp.ToFloat | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d, y?: 0-d
    def cdf(
        self: _BaseDist0[npc.number], x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CDFMethod | None = None
    ) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def cdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def cdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def cdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def cdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict3D,
        y: _ToFloatMax3D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def cdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def cdf(
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1],
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def cdf(
        self: _BaseDist0[npc.number], x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def cdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def cdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloat0ND,
        y: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def cdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def cdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def cdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def cdf(
        self: _BaseDist1[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def cdf(
        self: _BaseDist1[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist2[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def cdf(
        self: _BaseDist2[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def cdf(
        self: _BaseDist2[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def cdf(
        self: _BaseDist2[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def cdf(
        self: _BaseDist3[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def cdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d, y?: 0-d
    def logcdf(
        self: _BaseDistribution[Any, _ShapeT1],
        x: onp.ToFloat,
        y: onp.ToFloat | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d, y?: 0-d
    def logcdf(
        self: _BaseDist0[npc.number], x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CDFMethod | None = None
    ) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def logcdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def logcdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def logcdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict3D,
        y: _ToFloatMax3D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def logcdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def logcdf(
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1],
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def logcdf(
        self: _BaseDist0[npc.number], x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def logcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def logcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloat0ND,
        y: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def logcdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def logcdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def logcdf(
        self: _BaseDist1[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def logcdf(
        self: _BaseDist1[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist2[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist2[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def logcdf(
        self: _BaseDist2[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def logcdf(
        self: _BaseDist2[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def logcdf(
        self: _BaseDist3[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def logcdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d, y?: 0-d
    def ccdf(
        self: _BaseDistribution[Any, _ShapeT1],
        x: onp.ToFloat,
        y: onp.ToFloat | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d, y?: 0-d
    def ccdf(
        self: _BaseDist0[npc.number], x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CCDFMethod | None = None
    ) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def ccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def ccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def ccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def ccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict3D,
        y: _ToFloatMax3D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def ccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def ccdf(
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1],
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def ccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def ccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def ccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloat0ND,
        y: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def ccdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def ccdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def ccdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def ccdf(
        self: _BaseDist1[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def ccdf(
        self: _BaseDist1[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist2[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def ccdf(
        self: _BaseDist2[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def ccdf(
        self: _BaseDist2[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def ccdf(
        self: _BaseDist2[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def ccdf(
        self: _BaseDist3[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def ccdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, x: 0-d, y?: 0-d
    def logccdf(
        self: _BaseDistribution[Any, _ShapeT1],
        x: onp.ToFloat,
        y: onp.ToFloat | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d, y?: 0-d
    def logccdf(
        self: _BaseDist0[npc.number], x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CCDFMethod | None = None
    ) -> np.float64: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def logccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def logccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def logccdf(
        self: _BaseDist0[npc.number],
        x: onp.ToFloatStrict3D,
        y: _ToFloatMax3D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def logccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def logccdf(
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1],
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def logccdf(
        self: _BaseDist0[npc.number], x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def logccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def logccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number],
        x: _ToFloat0ND,
        y: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def logccdf(
        self: _BaseDist1[npc.number],
        x: onp.ToFloatStrict2D,
        y: _ToFloatMax2D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def logccdf(
        self: _BaseDist1[npc.number], x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def logccdf(
        self: _BaseDist1[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def logccdf(
        self: _BaseDist1[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist2[npc.number],
        x: onp.ToFloatStrict1D,
        y: _ToFloatMax1D | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist2[npc.number], x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def logccdf(
        self: _BaseDist2[npc.number], x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def logccdf(
        self: _BaseDist2[npc.number], x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def logccdf(
        self: _BaseDist3[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def logccdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    # NOTE: Apart from the `method` type, the signatures of `i[log]cdf` and `i[log]ccdf` are equivalent to those of `[log]pdf`
    @override
    @overload  # self: T1-d, p: 0-d
    def icdf(
        self: _BaseDistribution[Any, _ShapeT1], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: 0-d
    def icdf(self: _BaseDist0[npc.number], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, p: 1-d
    def icdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, p: 2-d
    def icdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, p: 3-d
    def icdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, p: T1-d
    def icdf(
        self: _BaseDist0[npc.number], p: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def icdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], p: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, p: 1-d
    def icdf(self: _BaseDist1[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, p: 2-d
    def icdf(self: _BaseDist1[npc.number], p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, p: >=-d
    def icdf(self: _BaseDist1[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, p: 1-d
    def icdf(self: _BaseDist2[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, p: >=1-d
    def icdf(self: _BaseDist2[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, p: >=1-d
    def icdf(self: _BaseDist3[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def icdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], p: _ToFloat0ND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, logp: 0-d
    def ilogcdf(
        self: _BaseDistribution[Any, _ShapeT1], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, logp: 0-d
    def ilogcdf(self: _BaseDist0[npc.number], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, logp: 1-d
    def ilogcdf(self: _BaseDist0[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, logp: 2-d
    def ilogcdf(self: _BaseDist0[npc.number], logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, logp: 3-d
    def ilogcdf(self: _BaseDist0[npc.number], logp: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, logp: T1-d
    def ilogcdf(
        self: _BaseDist0[npc.number], logp: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def ilogcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], logp: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, logp: 1-d
    def ilogcdf(self: _BaseDist1[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, logp: 2-d
    def ilogcdf(self: _BaseDist1[npc.number], logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, logp: >=-d
    def ilogcdf(self: _BaseDist1[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, logp: 1-ds
    def ilogcdf(self: _BaseDist2[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, logp: >=1-d
    def ilogcdf(self: _BaseDist2[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, logp: >=1-d
    def ilogcdf(self: _BaseDist3[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def ilogcdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], logp: _ToFloat0ND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, p: 0-d
    def iccdf(
        self: _BaseDistribution[Any, _ShapeT1], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: 0-d
    def iccdf(self: _BaseDist0[npc.number], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, p: 1-d
    def iccdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, p: 2-d
    def iccdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, p: 3-d
    def iccdf(self: _BaseDist0[npc.number], p: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, p: T1-d
    def iccdf(
        self: _BaseDist0[npc.number], p: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def iccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], p: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, p: 1-d
    def iccdf(self: _BaseDist1[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, p: 2-d
    def iccdf(self: _BaseDist1[npc.number], p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, p: >=-d
    def iccdf(self: _BaseDist1[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, p: 1-d
    def iccdf(self: _BaseDist2[npc.number], p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, p: >=1-d
    def iccdf(self: _BaseDist2[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, p: >=1-d
    def iccdf(self: _BaseDist3[npc.number], p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def iccdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], p: _ToFloat0ND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

    #
    @override
    @overload  # self: T1-d, logp: 0-d
    def ilogccdf(
        self: _BaseDistribution[Any, _ShapeT1], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, logp: 0-d
    def ilogccdf(self: _BaseDist0[npc.number], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> np.float64: ...
    @overload  # self: 0-d, logp: 1-d
    def ilogccdf(
        self: _BaseDist0[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, logp: 2-d
    def ilogccdf(
        self: _BaseDist0[npc.number], logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, logp: 3-d
    def ilogccdf(
        self: _BaseDist0[npc.number], logp: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, logp: T1-d
    def ilogccdf(
        self: _BaseDist0[npc.number], logp: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, q: >=1-d
    def ilogccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0[npc.number], logp: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, logp: 1-d
    def ilogccdf(
        self: _BaseDist1[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, logp: 2-d
    def ilogccdf(
        self: _BaseDist1[npc.number], logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, logp: >=-d
    def ilogccdf(self: _BaseDist1[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, logp: 1-d
    def ilogccdf(
        self: _BaseDist2[npc.number], logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, logp: >=1-d
    def ilogccdf(self: _BaseDist2[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, logp: >=1-d
    def ilogccdf(self: _BaseDist3[npc.number], logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def ilogccdf(  # ty:ignore[invalid-method-override]
        self: _BaseDist1N[npc.number], logp: _ToFloat0ND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[tuple[Any, ...]]: ...

#
class UnivariateDistribution(_BaseDistribution[_XT_co, _ShapeT0_co], Generic[_XT_co, _ShapeT0_co]):
    __array_priority__: ClassVar = 1
    _parameterizations: ClassVar[Sequence[_Parameterization]]
    _not_implemented: Final[str]
    _original_parameters: dict[str, _XT_co | onp.ArrayND[_XT_co, _ShapeT0_co]]
    _variable: _Parameter

    def __init__(
        self,
        *,
        tol: float | _Null | None = ...,
        validation_policy: _ValidationPolicy = None,
        cache_policy: _CachePolicy = None,
        **parameters: onp.ToFloat | None,
    ) -> None: ...

    #
    def reset_cache(self) -> None: ...

    #
    @property
    def tol(self, /) -> float | _Null | None: ...
    @tol.setter
    def tol(self, tol: float | _Null | None, /) -> None: ...

    #
    @property
    def cache_policy(self, /) -> _CachePolicy: ...
    @cache_policy.setter
    def cache_policy(self, cache_policy: _CachePolicy, /) -> None: ...

    #
    @property
    def validation_policy(self, /) -> _ValidationPolicy: ...
    @validation_policy.setter
    def validation_policy(self, validation_policy: _ValidationPolicy, /) -> None: ...

    #
    def plot(
        self,
        /,
        x: _PlotQuantity = "x",
        y: _PlotQuantity | None = None,
        *,
        t: tuple[_PlotQuantity, onp.ToJustFloat, onp.ToJustFloat] | None = None,
        ax: _AxesT | None = None,
    ) -> _AxesT: ...

    #
    def __neg__(self, /) -> _LinDist[Self, _FloatT_co, _ShapeT_co]: ...
    def __abs__(self, /) -> _FoldDist[Self, _FloatT_co, _ShapeT_co]: ...

    #
    @overload
    def __add__(self, x: float | npc.integer | np.bool, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __add__(self, x: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __add__(self, x: onp.ToFloat, /) -> _LinDist[Self, npc.floating, _ShapeT_co]: ...
    @overload
    def __add__(self: _DistT0, x: onp.CanArrayND[_FloatT, _ShapeT1], /) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __add__[DistT_1: _Dist[onp.AtMost1D]](self: DistT_1, x: onp.ToFloatStrict1D, /) -> _LinDist[DistT_1, np.float64, _1D]: ...
    @overload
    def __add__[DistT_2: _Dist[onp.AtMost2D]](self: DistT_2, x: onp.ToFloatStrict2D, /) -> _LinDist[DistT_2, np.float64, _2D]: ...
    @overload
    def __add__[DistT_3: _Dist[onp.AtMost3D]](self: DistT_3, x: onp.ToFloatStrict3D, /) -> _LinDist[DistT_3, np.float64, _3D]: ...
    @overload
    def __add__(self, x: onp.ToFloatND, /) -> _LinDist[Self, np.float64, tuple[Any, ...]]: ...
    __radd__ = __add__

    #
    @overload
    def __sub__(self, lshift: float | npc.integer | np.bool, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __sub__(self, lshift: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __sub__(self, lshift: onp.ToFloat, /) -> _LinDist[Self, npc.floating, _ShapeT_co]: ...
    @overload
    def __sub__(
        self: _DistT0, lshift: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __sub__[DistT_1: _Dist[onp.AtMost1D]](
        self: DistT_1, lshift: onp.ToFloatStrict1D, /
    ) -> _LinDist[DistT_1, np.float64, _1D]: ...
    @overload
    def __sub__[DistT_2: _Dist[onp.AtMost2D]](
        self: DistT_2, lshift: onp.ToFloatStrict2D, /
    ) -> _LinDist[DistT_2, np.float64, _2D]: ...
    @overload
    def __sub__[DistT_3: _Dist[onp.AtMost3D]](
        self: DistT_3, lshift: onp.ToFloatStrict3D, /
    ) -> _LinDist[DistT_3, np.float64, _3D]: ...
    @overload
    def __sub__(self, lshift: onp.ToFloatND, /) -> _LinDist[Self, np.float64, tuple[Any, ...]]: ...
    __rsub__ = __sub__

    #
    @overload
    def __mul__(self, scale: float | npc.integer | np.bool, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __mul__(self, scale: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __mul__(self, scale: onp.ToFloat, /) -> _LinDist[Self, npc.floating, _ShapeT_co]: ...
    @overload
    def __mul__(
        self: _DistT0, scale: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __mul__[DistT_1: _Dist[onp.AtMost1D]](
        self: DistT_1, scale: onp.ToFloatStrict1D, /
    ) -> _LinDist[DistT_1, np.float64, _1D]: ...
    @overload
    def __mul__[DistT_2: _Dist[onp.AtMost2D]](
        self: DistT_2, scale: onp.ToFloatStrict2D, /
    ) -> _LinDist[DistT_2, np.float64, _2D]: ...
    @overload
    def __mul__[DistT_3: _Dist[onp.AtMost3D]](
        self: DistT_3, scale: onp.ToFloatStrict3D, /
    ) -> _LinDist[DistT_3, np.float64, _3D]: ...
    @overload
    def __mul__(self, scale: onp.ToFloatND, /) -> _LinDist[Self, np.float64, tuple[Any, ...]]: ...
    __rmul__ = __mul__

    #
    def __pow__(self, exp: onp.ToInt, /) -> MonotonicTransformedDistribution[Self, _ShapeT_co]: ...
    __rpow__ = __pow__

    #
    @overload
    def __truediv__(self, iscale: float | npc.integer | np.bool, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __truediv__(self, iscale: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __truediv__(self, iscale: onp.ToFloat, /) -> _LinDist[Self, npc.floating, _ShapeT_co]: ...
    @overload
    def __truediv__(
        self: _DistT0, iscale: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __truediv__[DistT_1: _Dist[onp.AtMost1D]](
        self: DistT_1, iscale: onp.ToFloatStrict1D, /
    ) -> _LinDist[DistT_1, np.float64, _1D]: ...
    @overload
    def __truediv__[DistT_2: _Dist[onp.AtMost2D]](
        self: DistT_2, iscale: onp.ToFloatStrict2D, /
    ) -> _LinDist[DistT_2, np.float64, _2D]: ...
    @overload
    def __truediv__[DistT_3: _Dist[onp.AtMost3D]](
        self: DistT_3, iscale: onp.ToFloatStrict3D, /
    ) -> _LinDist[DistT_3, np.float64, _3D]: ...
    @overload
    def __truediv__(self, iscale: onp.ToFloatND, /) -> _LinDist[Self, np.float64, tuple[Any, ...]]: ...
    __rtruediv__ = __truediv__

#
class ContinuousDistribution(UnivariateDistribution[_FloatT_co, _ShapeT_co], Generic[_FloatT_co, _ShapeT_co]): ...

#
class DiscreteDistribution(UnivariateDistribution[_RealT_co, _ShapeT_co], Generic[_RealT_co, _ShapeT_co]): ...

# 7 years of asking and >400 upvotes, but still no higher-kinded typing support: https://github.com/python/typing/issues/548
class TransformedDistribution(ContinuousDistribution[_FloatT_co, _ShapeT_co], Generic[_DistT_co, _FloatT_co, _ShapeT_co]):
    _dist: _DistT_co  # readonly

    def __init__(
        self: _TransDist[ContinuousDistribution[_FloatT, _ShapeT], _FloatT, _ShapeT],  # nice trick, eh?
        X: _DistT_co,
        /,
        *args: Never,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

class ShiftedScaledDistribution(
    TransformedDistribution[_DistT_co, _FloatT_co, _ShapeT_co], Generic[_DistT_co, _FloatT_co, _ShapeT_co]
):
    _loc_domain: ClassVar[_RealInterval] = ...
    _loc_param: ClassVar[_RealParameter] = ...
    _scale_domain: ClassVar[_RealInterval] = ...
    _scale_param: ClassVar[_RealParameter] = ...

    loc: _ParamField[_FloatT_co, _ShapeT_co]
    scale: _ParamField[_FloatT_co, _ShapeT_co]

class FoldedDistribution(TransformedDistribution[_DistT_co, _FloatT_co, _ShapeT_co], Generic[_DistT_co, _FloatT_co, _ShapeT_co]):
    @overload
    def __init__[DistT0: _Dist[_0D]](
        self: _FoldDist[DistT0, np.float64, _0D], X: DistT0, /, *args: Never, **kwargs: Unpack[_DistOpts]
    ) -> None: ...
    @overload
    def __init__[DistT1: _Dist[_1D]](
        self: _FoldDist[DistT1, np.float64, _1D], X: DistT1, /, *args: Never, **kwargs: Unpack[_DistOpts]
    ) -> None: ...
    @overload
    def __init__[DistT2: _Dist[_2D]](
        self: _FoldDist[DistT2, np.float64, _2D], X: DistT2, /, *args: Never, **kwargs: Unpack[_DistOpts]
    ) -> None: ...
    @overload
    def __init__[DistT3: _Dist[_3D]](
        self: _FoldDist[DistT3, np.float64, _3D], X: DistT3, /, *args: Never, **kwargs: Unpack[_DistOpts]
    ) -> None: ...
    @overload
    def __init__[DistT: _Dist[tuple[int, ...]]](
        self: _FoldDist[DistT, np.float64, _ND], X: DistT, /, *args: Never, **kwargs: Unpack[_DistOpts]
    ) -> None: ...

class TruncatedDistribution(TransformedDistribution[_DistT_co, np.float64, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]):
    _lb_domain: ClassVar[_RealInterval] = ...
    _lb_param: ClassVar[_RealParameter] = ...
    _ub_domain: ClassVar[_RealInterval] = ...
    _ub_param: ClassVar[_RealParameter] = ...

    lb: _ParamField[np.float64, _ShapeT_co]
    ub: _ParamField[np.float64, _ShapeT_co]

    @overload
    def __init__[DistT0: _Dist[_0D]](
        self: _TruncDist[DistT0, _0D],
        X: DistT0,
        /,
        *args: Never,
        lb: onp.ToFloat = ...,
        ub: onp.ToFloat = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT1: _Dist[_1D]](
        self: _TruncDist[DistT1, _1D],
        X: DistT1,
        /,
        *args: Never,
        lb: _ToFloatMax1D = ...,
        ub: _ToFloatMax1D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT2: _Dist[_2D]](
        self: _TruncDist[DistT2, _2D],
        X: DistT2,
        /,
        *args: Never,
        lb: _ToFloatMax2D = ...,
        ub: _ToFloatMax2D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT3: _Dist[_3D]](
        self: _TruncDist[DistT3, _3D],
        X: DistT3,
        /,
        *args: Never,
        lb: _ToFloatMax3D = ...,
        ub: _ToFloatMax3D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT: _Dist[tuple[int, ...]]](
        self: _TruncDist[DistT, _ND],
        X: DistT,
        /,
        *args: Never,
        lb: _ToFloat0ND = ...,
        ub: _ToFloat0ND = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

# always float64 or longdouble
class OrderStatisticDistribution(TransformedDistribution[_DistT_co, np.float64, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]):
    _r_domain: ClassVar[_IntegerInterval] = ...
    _r_param: ClassVar[_RealParameter] = ...

    _n_domain: ClassVar[_IntegerInterval] = ...
    _n_param: ClassVar[_RealParameter] = ...

    @overload
    def __init__[DistT0: _Dist[_0D]](
        self: OrderStatisticDistribution[DistT0, _0D],
        dist: DistT0,
        /,
        *args: Never,
        r: onp.ToJustInt,
        n: onp.ToJustInt,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT1: _Dist[_1D]](
        self: OrderStatisticDistribution[DistT1, _1D],
        dist: DistT1,
        /,
        *args: Never,
        r: _ToJustIntMax1D,
        n: _ToJustIntMax1D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT2: _Dist[_2D]](
        self: OrderStatisticDistribution[DistT2, _2D],
        dist: DistT2,
        /,
        *args: Never,
        r: _ToJustIntMax2D,
        n: _ToJustIntMax2D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT3: _Dist[_3D]](
        self: OrderStatisticDistribution[DistT3, _3D],
        dist: DistT3,
        /,
        *args: Never,
        r: _ToJustIntMax3D,
        n: _ToJustIntMax3D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__[DistT: _Dist[tuple[int, ...]]](
        self: OrderStatisticDistribution[DistT, _ND],
        X: DistT,
        /,
        *args: Never,
        r: _ToJustIntMaxND,
        n: _ToJustIntMaxND,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

# without HKT there's no reasonable way tot determine the floating scalar type
class MonotonicTransformedDistribution(
    TransformedDistribution[_DistT_co, np.float64, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]
):
    _g: Final[_Elementwise[np.float64]]
    _h: Final[_Elementwise[np.float64]]
    _dh: Final[_Elementwise[np.float64]]
    _logdh: Final[_Elementwise[np.float64]]
    _increasing: Final[bool]
    _repr_pattern: Final[str]
    _str_pattern: Final[str]

    def __init__(
        self: MonotonicTransformedDistribution[_CDist[_ShapeT], _ShapeT],
        X: _DistT_co,
        /,
        *args: Never,
        g: _Elementwise[np.float64],
        h: _Elementwise[np.float64],
        dh: _Elementwise[np.float64],
        logdh: _Elementwise[np.float64] | None = None,
        increasing: bool = True,
        repr_pattern: str | None = None,
        str_pattern: str | None = None,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

class Mixture(_BaseDistribution[_FloatT_co, _0D], Generic[_FloatT_co]):
    _shape: _0D
    _dtype: np.dtype[_FloatT_co]
    _components: Sequence[_CDist0[_FloatT_co]]
    _weights: onp.Array1D[_FloatT_co]
    validation_policy: None

    @property
    def components(self, /) -> list[_CDist0[_FloatT_co]]: ...
    @property
    def weights(self, /) -> onp.Array1D[_FloatT_co]: ...
    #
    def __init__(self, /, components: Sequence[_CDist0[_FloatT_co]], *, weights: onp.ToFloat1D | None = None) -> None: ...
    #
    @override
    def kurtosis(self, /, *, method: _SMomentMethod | None = None) -> np.float64: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

    # always raises NotImplementedError`
    @override
    def lmoment(self, order: int = 1, *, standardize: bool = True, method: _LMomentMethod | None = None) -> Never: ...

###

# still waiting on the intersection type PEP...

@overload
def truncate[DistT0: _Dist[_0D]](X: DistT0, lb: onp.ToFloat = ..., ub: onp.ToFloat = ...) -> _TruncDist[DistT0, _0D]: ...
@overload
def truncate[DistT1: _Dist[_1D]](X: DistT1, lb: _ToFloatMax1D = ..., ub: _ToFloatMax1D = ...) -> _TruncDist[DistT1, _1D]: ...
@overload
def truncate[DistT2: _Dist[_2D]](X: DistT2, lb: _ToFloatMax2D = ..., ub: _ToFloatMax2D = ...) -> _TruncDist[DistT2, _2D]: ...
@overload
def truncate[DistT3: _Dist[_3D]](X: DistT3, lb: _ToFloatMax3D = ..., ub: _ToFloatMax3D = ...) -> _TruncDist[DistT3, _3D]: ...
@overload
def truncate[DistT: _Dist[tuple[int, ...]]](X: DistT, lb: _ToFloat0ND = ..., ub: _ToFloat0ND = ...) -> _TruncDist[DistT, _ND]: ...

#
@overload
def order_statistic[DistT0: _Dist[_0D]](
    X: DistT0, /, *, r: onp.ToJustInt, n: onp.ToJustInt
) -> OrderStatisticDistribution[DistT0, _0D]: ...
@overload
def order_statistic[DistT1: _Dist[_1D]](
    X: DistT1, /, *, r: _ToJustIntMax1D, n: _ToJustIntMax1D
) -> OrderStatisticDistribution[DistT1, _1D]: ...
@overload
def order_statistic[DistT2: _Dist[_2D]](
    X: DistT2, /, *, r: _ToJustIntMax2D, n: _ToJustIntMax2D
) -> OrderStatisticDistribution[DistT2, _2D]: ...
@overload
def order_statistic[DistT3: _Dist[_3D]](
    X: DistT3, /, *, r: _ToJustIntMax3D, n: _ToJustIntMax3D
) -> OrderStatisticDistribution[DistT3, _3D]: ...
@overload
def order_statistic[DistT: _Dist[tuple[int, ...]]](
    X: DistT, /, *, r: _ToJustIntMaxND, n: _ToJustIntMaxND
) -> OrderStatisticDistribution[DistT, _ND]: ...

#
@overload
def abs[DistT0: _Dist[_0D]](X: DistT0, /) -> _FoldDist[DistT0, np.float64, _0D]: ...
@overload
def abs[DistT1: _Dist[_1D]](X: DistT1, /) -> _FoldDist[DistT1, np.float64, _1D]: ...
@overload
def abs[DistT2: _Dist[_2D]](X: DistT2, /) -> _FoldDist[DistT2, np.float64, _2D]: ...
@overload
def abs[DistT3: _Dist[_3D]](X: DistT3, /) -> _FoldDist[DistT3, np.float64, _3D]: ...
@overload
def abs[DistT: _Dist[tuple[int, ...]]](X: DistT, /) -> _FoldDist[DistT, np.float64, _ND]: ...

#
@overload
def exp[DistT0: _Dist[_0D]](X: DistT0, /) -> MonotonicTransformedDistribution[DistT0, _0D]: ...
@overload
def exp[DistT1: _Dist[_1D]](X: DistT1, /) -> MonotonicTransformedDistribution[DistT1, _1D]: ...
@overload
def exp[DistT2: _Dist[_2D]](X: DistT2, /) -> MonotonicTransformedDistribution[DistT2, _2D]: ...
@overload
def exp[DistT3: _Dist[_3D]](X: DistT3, /) -> MonotonicTransformedDistribution[DistT3, _3D]: ...
@overload
def exp[DistT: _Dist[tuple[int, ...]]](X: DistT, /) -> MonotonicTransformedDistribution[DistT, _ND]: ...

#
@overload
def log[DistT0: _Dist[_0D]](X: DistT0, /) -> MonotonicTransformedDistribution[DistT0, _0D]: ...
@overload
def log[DistT1: _Dist[_1D]](X: DistT1, /) -> MonotonicTransformedDistribution[DistT1, _1D]: ...
@overload
def log[DistT2: _Dist[_2D]](X: DistT2, /) -> MonotonicTransformedDistribution[DistT2, _2D]: ...
@overload
def log[DistT3: _Dist[_3D]](X: DistT3, /) -> MonotonicTransformedDistribution[DistT3, _3D]: ...
@overload
def log[DistT: _Dist[tuple[int, ...]]](X: DistT, /) -> MonotonicTransformedDistribution[DistT, _ND]: ...

###
# make_distribution

@type_check_only
class _CustomDistributionMixin:
    def __init__(
        self,
        *,
        tol: float | _Null | None = ...,
        validation_policy: _ValidationPolicy = None,
        cache_policy: _CachePolicy = None,
        **parameters: onp.ToFloat | onp.ToFloatND | None,
    ) -> None: ...

    # NOTE: The signatures of `pdf`, `logpdf`, `pmf`, and `logpmf` are equivalent
    @overload  # x: 0-d
    def pdf(self, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # x: >0-d
    def pdf(self, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> onp.ArrayND[np.float64]: ...

    #
    @overload  # x: 0-d
    def logpdf(self, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # x: >0-d
    def logpdf(self, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> onp.ArrayND[np.float64]: ...

    #
    @overload  # x: 0-d
    def pmf(self, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # x: >0-d
    def pmf(self, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> onp.ArrayND[np.float64]: ...

    #
    @overload  # x: 0-d
    def logpmf(self, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # x: >0-d
    def logpmf(self, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> onp.ArrayND[np.float64]: ...

@type_check_only
class _CustomContinuousDistribution(ContinuousDistribution[np.float64, _ShapeT_co]): ...

@type_check_only
class _CustomDiscreteDistribution(DiscreteDistribution[np.float64, _ShapeT_co]): ...

# Workarounds for the lack of higher-kinded typing (HKT) support in Python.
# See https://github.com/python/typing/issues/548 and https://github.com/jorenham/hkt-survey for details.

@type_check_only
class _CustomContinuousDistributionKind(Protocol):
    @overload
    def __call__(self, /) -> _CustomContinuousDistribution[tuple[()]]: ...  # pyright: ignore[reportOverlappingOverload]
    @overload
    def __call__(self, /, **parameters: onp.ToFloat) -> _CustomContinuousDistribution[tuple[()]]: ...
    @overload
    def __call__(
        self, /, **parameters: onp.ArrayND[_Real, tuple[Never, Never, Never]]
    ) -> _CustomContinuousDistribution[tuple[()]] | _CustomContinuousDistribution[tuple[Any, ...]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatStrict1D) -> _CustomContinuousDistribution[tuple[int]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatStrict2D) -> _CustomContinuousDistribution[tuple[int, int]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatND) -> _CustomContinuousDistribution[tuple[int, *tuple[Any, ...]]]: ...

@type_check_only
class _CustomDiscreteDistributionKind(Protocol):
    @overload
    def __call__(self, /) -> _CustomDiscreteDistribution[tuple[()]]: ...  # pyright: ignore[reportOverlappingOverload]
    @overload
    def __call__(self, /, **parameters: onp.ToFloat) -> _CustomDiscreteDistribution[tuple[()]]: ...
    @overload
    def __call__(
        self, /, **parameters: onp.ArrayND[_Real, tuple[Never, Never, Never]]
    ) -> _CustomDiscreteDistribution[tuple[()]] | _CustomDiscreteDistribution[tuple[Any, ...]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatStrict1D) -> _CustomDiscreteDistribution[tuple[int]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatStrict2D) -> _CustomDiscreteDistribution[tuple[int, int]]: ...
    @overload
    def __call__(self, /, **parameters: onp.ToFloatND) -> _CustomDiscreteDistribution[tuple[int, *tuple[Any, ...]]]: ...

@overload
def make_distribution(dist: _DuckDistributionType) -> type[_CustomContinuousDistribution[tuple[()]]]: ...
@overload
def make_distribution(dist: rv_continuous) -> _CustomContinuousDistributionKind: ...
@overload
def make_distribution(dist: rv_discrete) -> _CustomDiscreteDistributionKind: ...
