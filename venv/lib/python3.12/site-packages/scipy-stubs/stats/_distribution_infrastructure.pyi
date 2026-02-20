# mypy: disable-error-code="override, explicit-override"
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
    TypeAlias,
    TypedDict,
    final,
    overload,
    type_check_only,
)
from typing_extensions import ParamSpec, TypeAliasType, TypeIs, TypeVar, Unpack, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._distn_infrastructure import rv_continuous
from ._probability_distribution import _ProbabilityDistribution
from ._qmc import QMCEngine

__all__ = ["Mixture", "abs", "exp", "log", "make_distribution", "order_statistic", "truncate"]

###

_Tss = ParamSpec("_Tss", default=...)

_FloatT = TypeVar("_FloatT", bound=_Float, default=_Float)
_FloatT_co = TypeVar("_FloatT_co", bound=_Float, default=_Float, covariant=True)
_IntT_co = TypeVar("_IntT_co", bound=_Int, default=_Int, covariant=True)

_RealT = TypeVar("_RealT", bound=_Float | _Int, default=_Float | _Int)
_RealT_co = TypeVar("_RealT_co", bound=_Float | _Int, default=_Float | _Int, covariant=True)

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[Any, ...])
_ShapeT1 = TypeVar("_ShapeT1", bound=tuple[int, *tuple[int, ...]], default=tuple[Any, ...])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

_DistT0 = TypeVar("_DistT0", bound=_Dist[_0D])
_DistT1 = TypeVar("_DistT1", bound=_Dist[_1D])
_DistT_1 = TypeVar("_DistT_1", bound=_Dist[onp.AtMost1D])
_DistT2 = TypeVar("_DistT2", bound=_Dist[_2D])
_DistT_2 = TypeVar("_DistT_2", bound=_Dist[onp.AtMost2D])
_DistT3 = TypeVar("_DistT3", bound=_Dist[_3D])
_DistT_3 = TypeVar("_DistT_3", bound=_Dist[onp.AtMost3D])
_DistT = TypeVar("_DistT", bound=_Dist[tuple[int, ...]])
_DistT_co = TypeVar("_DistT_co", bound=_Dist[tuple[int, ...]], default=UnivariateDistribution, covariant=True)

_AxesT = TypeVar("_AxesT", bound=_Axes, default=Any)

_ParameterEndpoint: TypeAlias = onp.ToFloat | Callable[..., onp.ToFloat] | str
_ParameterTuple: TypeAlias = tuple[_ParameterEndpoint, _ParameterEndpoint]

# NOTE: `TypedDict` with `NotRequired` cannot be used, because `NotRequired` does not
# mean "not required": Other `TypedDict` MUST also include them to be assignable,
# so a `NotRequired` value is REQUIRED. Absolutely ridiculous...
# https://typing.python.org/en/latest/spec/typeddict.html#id4
_ParameterDict: TypeAlias = Mapping[str, tuple[Any, ...]]
_ParameterSpec: TypeAlias = _ParameterDict | _ParameterTuple

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

_DuckDistributionType: TypeAlias = type[_DuckDistributionSingle | _DuckDistributionMulti]

###

_Int: TypeAlias = npc.integer
_Float: TypeAlias = npc.floating
_OutFloat: TypeAlias = np.float64 | np.longdouble

_0D: TypeAlias = tuple[()]  # noqa: PYI042
_1D: TypeAlias = tuple[int]  # noqa: PYI042
_2D: TypeAlias = tuple[int, int]  # noqa: PYI042
_3D: TypeAlias = tuple[int, int, int]  # noqa: PYI042
_ND: TypeAlias = tuple[int, ...]

_ToFloatMax1D: TypeAlias = onp.ToFloatStrict1D | onp.ToFloat
_ToFloatMax2D: TypeAlias = onp.ToFloatStrict2D | _ToFloatMax1D
_ToFloatMax3D: TypeAlias = onp.ToFloatStrict3D | _ToFloatMax2D

_ToJustIntMax1D: TypeAlias = onp.ToJustIntStrict1D | onp.ToJustInt
_ToJustIntMax2D: TypeAlias = onp.ToJustIntStrict2D | _ToJustIntMax1D
_ToJustIntMax3D: TypeAlias = onp.ToJustIntStrict3D | _ToJustIntMax2D
_ToJustIntMaxND: TypeAlias = onp.ToJustIntND | onp.ToJustInt

_Null: TypeAlias = op.JustObject  # type of `_null`
_Axes: TypeAlias = object  # placeholder for `matplotlib.axes.Axes`

_DomainRegion: TypeAlias = L["domain", "typical"]
_DomainDrawType: TypeAlias = L["in", "out", "on", "nan"]
_ValidationPolicy: TypeAlias = L["skip_all"] | None
_CachePolicy: TypeAlias = L["no_cache"] | None
_PlotQuantity: TypeAlias = L["x", "pdf", "cdf", "ccdf", "icdf", "iccdf", "logpdf", "logcdf", "logccdf", "ilogcdf", "ilogccdf"]
_SMomentMethod: TypeAlias = L["formula", "general", "transform", "normalize", "cache"]
_CMomentMethod: TypeAlias = L["formula", "transform", "quadrature", "cache", "normalize"]
_RMomentMethod: TypeAlias = L["formula", "transform", "quadrature", "cache"]

_ParamValues: TypeAlias = Mapping[str, _ToFloat0ND]
_ToDomain: TypeAlias = tuple[onp.ToFloat | str, onp.ToFloat | str]
_ToTol: TypeAlias = op.JustFloat | _Null
_DrawProportions: TypeAlias = tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat, onp.ToFloat]
_Elementwise: TypeAlias = Callable[[onp.ArrayND[np.float64]], onp.ArrayND[_FloatT]]

_Dist: TypeAlias = UnivariateDistribution[Any, _ShapeT]
_CDist: TypeAlias = ContinuousDistribution[_Float, _ShapeT]
_CDist0: TypeAlias = ContinuousDistribution[_FloatT, _0D]
_TransDist: TypeAlias = TransformedDistribution[_DistT, _FloatT, _ShapeT]
_LinDist: TypeAlias = ShiftedScaledDistribution[_DistT, _FloatT, _ShapeT]
_FoldDist: TypeAlias = FoldedDistribution[_DistT, _FloatT, _ShapeT]
_TruncDist: TypeAlias = TruncatedDistribution[_DistT, _ShapeT]

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

class _Domain(abc.ABC, Generic[_XT_co]):
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    # NOTE: This is a `ClassVar[dict[str, float]]` that's overridden as instance attribute in `_SimpleDomain`.
    # https://github.com/scipy/scipy/pull/22139
    symbols: Mapping[str, str] = ...

    @override
    @abc.abstractmethod
    def __str__(self, /) -> str: ...
    @abc.abstractmethod
    def contains(self, /, x: onp.ArrayND[Any]) -> onp.ArrayND[np.bool_]: ...
    @abc.abstractmethod
    def draw(self, /, n: int) -> onp.ArrayND[_XT_co]: ...
    @abc.abstractmethod
    def get_numerical_endpoints(self, /, x: _ParamValues) -> tuple[onp.ArrayND[_OutFloat], onp.ArrayND[_OutFloat]]: ...

class _Interval(_Domain[_XT_co], Generic[_XT_co]):
    @override
    @abc.abstractmethod
    def __str__(self, /) -> str: ...

    #
    def __init__(self, /, endpoints: _ToDomain = ..., inclusive: tuple[bool, bool] = (False, False)) -> None: ...
    @override
    def get_numerical_endpoints(  # pyright: ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-param-name-override] # ty: ignore[invalid-method-override]
        self, /, parameter_values: _ParamValues
    ) -> tuple[onp.ArrayND[_OutFloat], onp.ArrayND[_OutFloat]]: ...
    @override
    def contains(  # pyright: ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-param-name-override] # ty: ignore[invalid-method-override]
        self, /, item: onp.ArrayND[_Int | _Float], parameter_values: _ParamValues | None = None
    ) -> onp.ArrayND[np.bool_]: ...
    @override
    def draw(  # pyright: ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]
        self,
        /,
        n: int,
        type_: _DomainDrawType,
        min: onp.ArrayND[_Float | _Int],
        max: onp.ArrayND[_Float | _Int],
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

_ValidateOut0D = TypeAliasType("_ValidateOut0D", tuple[_RealT, np.dtype[_RealT], onp.Array0D[np.bool_]], type_params=(_RealT,))
_ValidateOutND = TypeAliasType(
    "_ValidateOutND",
    tuple[onp.ArrayND[_RealT, _ShapeT1], np.dtype[_RealT], onp.ArrayND[np.bool_, _ShapeT1]],
    type_params=(_RealT, _ShapeT1),
)

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
    def validate(self, /, arr: onp.ToFloatND) -> _ValidateOutND[_RealT_co]: ...
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
    def validate(self, /, arr: onp.ToFloatND, parameter_values: _ParamValues) -> _ValidateOutND[_FloatT_co]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

class _Parameterization:
    parameters: Final[Mapping[str, _Parameter]]

    def __init__(self, /, *parameters: _Parameter) -> None: ...
    def __len__(self, /) -> int: ...
    def copy(self, /) -> Self: ...
    def matches(self, /, parameters: AbstractSet[str]) -> bool: ...
    def validation(self, /, parameter_values: Mapping[str, _Parameter]) -> tuple[onp.ArrayND[np.bool_], np.dtype[_Float]]: ...
    def draw(
        self,
        /,
        sizes: _ND | Sequence[_ND] | None = None,
        rng: onp.random.ToRNG | None = None,
        proportions: _DrawProportions | None = None,
        region: _DomainRegion = "domain",
    ) -> dict[str, onp.ArrayND[_Float]]: ...

###

_T = TypeVar("_T")
_Tuple2: TypeAlias = tuple[_T, _T]

_XT = TypeVar("_XT", bound=npc.number, default=npc.number)
_XT_co = TypeVar("_XT_co", bound=npc.number, default=np.float64, covariant=True)
_ShapeT0_co = TypeVar("_ShapeT0_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

_BaseDist0: TypeAlias = _BaseDistribution[_XT, tuple[()]]
_BaseDist1: TypeAlias = _BaseDistribution[_XT, tuple[int]]
_BaseDist2: TypeAlias = _BaseDistribution[_XT, tuple[int, int]]
_BaseDist3: TypeAlias = _BaseDistribution[_XT, tuple[int, int, int]]
_BaseDist1N: TypeAlias = _BaseDistribution[_XT, tuple[int, *tuple[int, ...]]]

_KurtosisConvention: TypeAlias = L["non-excess", "excess"]
_MedianMethod: TypeAlias = L["formula", "icdf"]
_ModeMethod: TypeAlias = L["formula", "optimization"]
_SampleMethod: TypeAlias = L["formula", "inverse_transform"]
_EntropyMethod: TypeAlias = L["formula", "logexp", "quadrature"]

_PXFMethod: TypeAlias = L["formula", "logexp"]
_CDFMethod: TypeAlias = L["formula", "logexp", "complement", "quadrature", "subtraction"]
_CCDFMethod: TypeAlias = L["formula", "logexp", "complement", "quadrature", "addition"]
_ICDFMethod: TypeAlias = L["formula", "complement", "inversion"]

_Float1D: TypeAlias = onp.Array1D[_Float]
_Float2D: TypeAlias = onp.Array2D[_Float]
_Float3D: TypeAlias = onp.Array3D[_Float]
_FloatND: TypeAlias = onp.ArrayND[_Float, _ShapeT1]

_Float1ND: TypeAlias = onp.Array[tuple[int, *tuple[Any, ...]], _Float]
# https://github.com/microsoft/pyright/issues/11127
_Float2ND: TypeAlias = onp.Array[tuple[int, int, *tuple[Any, ...]], _Float]  # pyright: ignore[reportInvalidTypeForm]
_Float3ND: TypeAlias = onp.Array[tuple[int, int, int, *tuple[Any, ...]], _Float]  # pyright: ignore[reportInvalidTypeForm]

_Complex: TypeAlias = np.complex128 | np.clongdouble
_ComplexND: TypeAlias = onp.ArrayND[_Complex, _ShapeT1]

_ToFloatND: TypeAlias = onp.CanArrayND[npc.floating | npc.integer | np.bool_, _ShapeT1]
_ToFloat0ND: TypeAlias = onp.ToFloat | onp.ToFloatND
_ToFloatMaxND: TypeAlias = _ToFloatND[_ShapeT1] | _ToFloatMax1D

_ToQRNG: TypeAlias = QMCEngine | onp.random.ToRNG | None

@type_check_only
class _BaseDistribution(_ProbabilityDistribution[_XT_co], Generic[_XT_co, _ShapeT0_co]):
    @overload
    def support(self: _BaseDist0[_XT], /) -> _Tuple2[_XT]: ...
    @overload
    def support(self: _BaseDistribution[_XT, _ShapeT1], /) -> _Tuple2[onp.Array[_ShapeT1, _XT]]: ...  # pyright: ignore[reportIncompatibleMethodOverride]

    #
    @overload
    def median(self: _BaseDist0[_XT], /, *, method: _MedianMethod | None = None) -> _XT: ...
    @overload
    def median(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _MedianMethod | None = None
    ) -> onp.Array[_ShapeT1, _XT]: ...

    #
    @overload
    def mode(self: _BaseDist0[_XT], /, *, method: _ModeMethod | None = None) -> _XT: ...
    @overload
    def mode(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _ModeMethod | None = None
    ) -> onp.Array[_ShapeT1, _XT]: ...

    #
    @overload
    def sample(
        self: _BaseDist0[_XT], /, shape: tuple[()] = (), *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> _XT: ...
    @overload
    def sample(
        self: _BaseDist0[_XT], /, shape: op.CanIndex, *, method: _SampleMethod | None = None, rng: _ToQRNG = None
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
        shape: op.CanIndex | Iterable[op.CanIndex],
        *,
        method: _SampleMethod | None = None,
        rng: _ToQRNG = None,
    ) -> onp.ArrayND[_XT, _ShapeT1] | onp.ArrayND[_XT]: ...  # first union type is needed on `numpy<2.1`
    @overload
    def sample(
        self, /, shape: op.CanIndex | Iterable[op.CanIndex], *, method: _SampleMethod | None = None, rng: _ToQRNG = None
    ) -> _XT_co | onp.ArrayND[_XT_co, _ShapeT1] | onp.ArrayND[_XT_co]: ...  # first union type is needed on `numpy<2.1`

    #
    @overload
    def mean(self: _BaseDist0[_XT], /, *, method: _RMomentMethod | None = None) -> _XT: ...
    @overload
    def mean(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _RMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @overload
    def variance(self: _BaseDist0[_XT], /, *, method: _CMomentMethod | None = None) -> _XT: ...
    @overload
    def variance(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _CMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @overload
    def standard_deviation(self: _BaseDist0[_XT], /, *, method: _CMomentMethod | None = None) -> _XT: ...
    @overload
    def standard_deviation(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _CMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @overload
    def skewness(self: _BaseDist0[_XT], /, *, method: _SMomentMethod | None = None) -> _XT: ...
    @overload
    def skewness(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1], /, *, method: _SMomentMethod | None = None
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @overload
    def kurtosis(
        self: _BaseDist0[_XT], /, *, method: _SMomentMethod | None = None, convention: _KurtosisConvention = "non-excess"
    ) -> _XT: ...
    @overload
    def kurtosis(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[_XT, _ShapeT1],
        /,
        *,
        method: _SMomentMethod | None = None,
        convention: _KurtosisConvention = "non-excess",
    ) -> onp.ArrayND[_XT, _ShapeT1]: ...

    #
    @overload
    def moment(
        self: _BaseDist0, /, order: onp.ToInt = 1, kind: L["raw"] = "raw", *, method: _RMomentMethod | None = None
    ) -> _Float: ...
    @overload
    def moment(self: _BaseDist0, /, order: onp.ToInt, kind: L["central"], *, method: _CMomentMethod | None = None) -> _Float: ...
    @overload
    def moment(
        self: _BaseDist0, /, order: onp.ToInt = 1, *, kind: L["central"], method: _CMomentMethod | None = None
    ) -> _Float: ...
    @overload
    def moment(
        self: _BaseDist0, /, order: onp.ToInt, kind: L["standardized"], *, method: _SMomentMethod | None = None
    ) -> _Float: ...
    @overload
    def moment(
        self: _BaseDist0, /, order: onp.ToInt = 1, *, kind: L["standardized"], method: _SMomentMethod | None = None
    ) -> _Float: ...
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
    def moment(  # pyright: ignore[reportIncompatibleMethodOverride]
        self: _BaseDistribution[Any, _ShapeT1],
        /,
        order: onp.ToInt = 1,
        *,
        kind: L["standardized"],
        method: _SMomentMethod | None = None,
    ) -> _FloatND[_ShapeT1]: ...

    #
    @overload
    def entropy(self: _BaseDist0, /, *, method: _EntropyMethod | None = None) -> _Float: ...
    @overload
    def entropy(self: _BaseDistribution[Any, _ShapeT1], /, *, method: _EntropyMethod | None = None) -> _FloatND[_ShapeT1]: ...

    #
    @overload
    def logentropy(self: _BaseDist0, /, *, method: _EntropyMethod | None = None) -> _Complex: ...
    @overload
    def logentropy(
        self: _BaseDistribution[Any, _ShapeT1], /, *, method: _EntropyMethod | None = None
    ) -> _ComplexND[_ShapeT1]: ...

    #
    # TODO(jorenham): Adjust these, depending on the result of https://github.com/scipy/scipy/issues/22145
    # NOTE: The signatures of `pdf`, `logpdf`, `pmf`, and `logpmf` are equivalent
    @overload  # self: T1-d, x: 0-d
    def pdf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def pdf(self: _BaseDist0, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, x: 1-d
    def pdf(self: _BaseDist0, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def pdf(self: _BaseDist0, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def pdf(self: _BaseDist0, x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def pdf(self: _BaseDist0, x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def pdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def pdf(self: _BaseDist1, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def pdf(self: _BaseDist1, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=-d
    def pdf(self: _BaseDist1, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def pdf(self: _BaseDist2, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def pdf(self: _BaseDist2, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def pdf(self: _BaseDist3, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def pdf(self: _BaseDist1N, x: _ToFloat0ND, /, *, method: _PXFMethod | None = None) -> _FloatND: ...

    #
    @overload  # self: T1-d, x: 0-d
    def logpdf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def logpdf(self: _BaseDist0, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, x: 1-d
    def logpdf(self: _BaseDist0, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def logpdf(self: _BaseDist0, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def logpdf(self: _BaseDist0, x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def logpdf(self: _BaseDist0, x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def logpdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def logpdf(self: _BaseDist1, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def logpdf(self: _BaseDist1, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d
    def logpdf(self: _BaseDist1, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def logpdf(self: _BaseDist2, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def logpdf(self: _BaseDist2, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def logpdf(self: _BaseDist3, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def logpdf(self: _BaseDist1N, x: _ToFloat0ND, /, *, method: _PXFMethod | None = None) -> _FloatND: ...

    #
    @overload  # self: T1-d, x: 0-d
    def pmf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def pmf(self: _BaseDist0, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, x: 1-d
    def pmf(self: _BaseDist0, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def pmf(self: _BaseDist0, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def pmf(self: _BaseDist0, x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def pmf(self: _BaseDist0, x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def pmf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def pmf(self: _BaseDist1, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def pmf(self: _BaseDist1, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=-d
    def pmf(self: _BaseDist1, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def pmf(self: _BaseDist2, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def pmf(self: _BaseDist2, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def pmf(self: _BaseDist3, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def pmf(self: _BaseDist1N, x: _ToFloat0ND, /, *, method: _PXFMethod | None = None) -> _FloatND: ...

    #
    @overload  # self: T1-d, x: 0-d
    def logpmf(
        self: _BaseDistribution[Any, _ShapeT1], x: onp.ToFloat, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: 0-d
    def logpmf(self: _BaseDist0, x: onp.ToFloat, /, *, method: _PXFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, x: 1-d
    def logpmf(self: _BaseDist0, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d
    def logpmf(self: _BaseDist0, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d
    def logpmf(self: _BaseDist0, x: onp.ToFloatStrict3D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d
    def logpmf(self: _BaseDist0, x: _ToFloatND[_ShapeT1], /, *, method: _PXFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d
    def logpmf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _PXFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d
    def logpmf(self: _BaseDist1, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d
    def logpmf(self: _BaseDist1, x: onp.ToFloatStrict2D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d
    def logpmf(self: _BaseDist1, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d
    def logpmf(self: _BaseDist2, x: onp.ToFloatStrict1D, /, *, method: _PXFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d
    def logpmf(self: _BaseDist2, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=1-d
    def logpmf(self: _BaseDist3, x: onp.ToFloatND, /, *, method: _PXFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def logpmf(self: _BaseDist1N, x: _ToFloat0ND, /, *, method: _PXFMethod | None = None) -> _FloatND: ...

    #
    # NOTE: Apart from the `method` type, the signatures of `[log]cdf` and `[log]ccdf` are equivalent
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
    def cdf(self: _BaseDist0, x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CDFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist0, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def cdf(self: _BaseDist0, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def cdf(
        self: _BaseDist0, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def cdf(self: _BaseDist0, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def cdf(
        self: _BaseDist0, x: onp.ToFloatStrict3D, y: _ToFloatMax3D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def cdf(self: _BaseDist0, x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def cdf(
        self: _BaseDist0, x: _ToFloatND[_ShapeT1], y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def cdf(
        self: _BaseDist0, x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def cdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0,
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def cdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloat0ND, y: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist1, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def cdf(self: _BaseDist1, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def cdf(
        self: _BaseDist1, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def cdf(self: _BaseDist1, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def cdf(
        self: _BaseDist1, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def cdf(self: _BaseDist1, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def cdf(
        self: _BaseDist2, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def cdf(self: _BaseDist2, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def cdf(
        self: _BaseDist2, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def cdf(self: _BaseDist2, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def cdf(
        self: _BaseDist3, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def cdf(
        self: _BaseDist1N, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND: ...

    #
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
        self: _BaseDist0, x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist0, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist0, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def logcdf(
        self: _BaseDist0, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def logcdf(
        self: _BaseDist0, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def logcdf(
        self: _BaseDist0, x: onp.ToFloatStrict3D, y: _ToFloatMax3D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def logcdf(
        self: _BaseDist0, x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def logcdf(
        self: _BaseDist0, x: _ToFloatND[_ShapeT1], y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def logcdf(
        self: _BaseDist0, x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def logcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0,
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def logcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloat0ND, y: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _CDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist1, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist1, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def logcdf(
        self: _BaseDist1, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def logcdf(
        self: _BaseDist1, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def logcdf(
        self: _BaseDist1, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def logcdf(self: _BaseDist1, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def logcdf(
        self: _BaseDist2, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def logcdf(
        self: _BaseDist2, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def logcdf(
        self: _BaseDist2, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def logcdf(self: _BaseDist2, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def logcdf(
        self: _BaseDist3, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def logcdf(
        self: _BaseDist1N, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CDFMethod | None = None
    ) -> _FloatND: ...

    #
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
        self: _BaseDist0, x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist0, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def ccdf(self: _BaseDist0, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def ccdf(
        self: _BaseDist0, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def ccdf(self: _BaseDist0, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def ccdf(
        self: _BaseDist0, x: onp.ToFloatStrict3D, y: _ToFloatMax3D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def ccdf(self: _BaseDist0, x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CCDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def ccdf(
        self: _BaseDist0, x: _ToFloatND[_ShapeT1], y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def ccdf(
        self: _BaseDist0, x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def ccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0,
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def ccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloat0ND, y: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist1, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def ccdf(self: _BaseDist1, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def ccdf(
        self: _BaseDist1, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def ccdf(self: _BaseDist1, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def ccdf(
        self: _BaseDist1, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def ccdf(self: _BaseDist1, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def ccdf(
        self: _BaseDist2, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def ccdf(self: _BaseDist2, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def ccdf(
        self: _BaseDist2, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def ccdf(self: _BaseDist2, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def ccdf(
        self: _BaseDist3, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def ccdf(
        self: _BaseDist1N, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND: ...

    #
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
        self: _BaseDist0, x: onp.ToFloat, y: onp.ToFloat | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float: ...
    @overload  # self: 0-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist0, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist0, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float1D: ...
    @overload  # self: 0-d, x: 2-d, y?: <=2-d
    def logccdf(
        self: _BaseDist0, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: <=2-d, y: 2-d
    def logccdf(
        self: _BaseDist0, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 0-d, x: 3-d, y?: <=3-d
    def logccdf(
        self: _BaseDist0, x: onp.ToFloatStrict3D, y: _ToFloatMax3D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: <=3-d, y: 3-d
    def logccdf(
        self: _BaseDist0, x: _ToFloatMax3D, y: onp.ToFloatStrict3D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 0-d, x: T1-d, y?: T1-d | <=1-d
    def logccdf(
        self: _BaseDist0, x: _ToFloatND[_ShapeT1], y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: T1-d | <=1-d, y: T1-d
    def logccdf(
        self: _BaseDist0, x: _ToFloatMaxND[_ShapeT1], y: _ToFloatND[_ShapeT1], /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, x: >=1-d, y?: >=0-d
    def logccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0,
        x: _ToFloatND[_ShapeT1] | onp.ToFloatND,
        y: _ToFloat0ND | None = None,
        /,
        *,
        method: _CCDFMethod | None = None,
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 0-d, x: >=0-d, y: >=1-d
    def logccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, x: _ToFloat0ND, y: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist1, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist1, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float2D: ...
    @overload  # self: 1-d, x: 2-d, y?: <=2-d
    def logccdf(
        self: _BaseDist1, x: onp.ToFloatStrict2D, y: _ToFloatMax2D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: <=2-d, y: 2-d
    def logccdf(
        self: _BaseDist1, x: _ToFloatMax2D, y: onp.ToFloatStrict2D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 1-d, x: >=1-d, y?: >=0-d
    def logccdf(
        self: _BaseDist1, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float2ND: ...
    @overload  # self: 1-d, x: >=0-d, y: >=1-d
    def logccdf(self: _BaseDist1, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, x: 1-d, y?: <=1-d
    def logccdf(
        self: _BaseDist2, x: onp.ToFloatStrict1D, y: _ToFloatMax1D | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: <=1-d, y: 1-d
    def logccdf(
        self: _BaseDist2, x: _ToFloatMax1D, y: onp.ToFloatStrict1D, /, *, method: _CCDFMethod | None = None
    ) -> _Float3D: ...
    @overload  # self: 2-d, x: >=1-d, y?: >=0-d
    def logccdf(
        self: _BaseDist2, x: onp.ToFloatND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: 2-d, x: >=0-d, y: >=1-d
    def logccdf(self: _BaseDist2, x: _ToFloat0ND, y: onp.ToFloatND, /, *, method: _CCDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, x: >=0-d, y?: >=0-d
    def logccdf(
        self: _BaseDist3, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _Float3ND: ...
    @overload  # self: >=1-d, x: >=0-d, y?: >=0-d
    def logccdf(
        self: _BaseDist1N, x: _ToFloat0ND, y: _ToFloat0ND | None = None, /, *, method: _CCDFMethod | None = None
    ) -> _FloatND: ...

    # NOTE: Apart from the `method` type, the signatures of `i[log]cdf` and `i[log]ccdf` are equivalent to those of `[log]pdf`
    @overload  # self: T1-d, p: 0-d
    def icdf(
        self: _BaseDistribution[Any, _ShapeT1], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: 0-d
    def icdf(self: _BaseDist0, p: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, p: 1-d
    def icdf(self: _BaseDist0, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, p: 2-d
    def icdf(self: _BaseDist0, p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, p: 3-d
    def icdf(self: _BaseDist0, p: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, p: T1-d
    def icdf(self: _BaseDist0, p: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def icdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, p: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, p: 1-d
    def icdf(self: _BaseDist1, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, p: 2-d
    def icdf(self: _BaseDist1, p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, p: >=-d
    def icdf(self: _BaseDist1, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, p: 1-d
    def icdf(self: _BaseDist2, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, p: >=1-d
    def icdf(self: _BaseDist2, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, p: >=1-d
    def icdf(self: _BaseDist3, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def icdf(self: _BaseDist1N, p: _ToFloat0ND, /, *, method: _ICDFMethod | None = None) -> _FloatND: ...
    #
    @overload  # self: T1-d, logp: 0-d
    def ilogcdf(
        self: _BaseDistribution[Any, _ShapeT1], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, logp: 0-d
    def ilogcdf(self: _BaseDist0, logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, logp: 1-d
    def ilogcdf(self: _BaseDist0, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, logp: 2-d
    def ilogcdf(self: _BaseDist0, logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, logp: 3-d
    def ilogcdf(self: _BaseDist0, logp: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, logp: T1-d
    def ilogcdf(self: _BaseDist0, logp: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def ilogcdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, logp: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, logp: 1-d
    def ilogcdf(self: _BaseDist1, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, logp: 2-d
    def ilogcdf(self: _BaseDist1, logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, logp: >=-d
    def ilogcdf(self: _BaseDist1, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, logp: 1-ds
    def ilogcdf(self: _BaseDist2, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, logp: >=1-d
    def ilogcdf(self: _BaseDist2, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, logp: >=1-d
    def ilogcdf(self: _BaseDist3, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def ilogcdf(self: _BaseDist1N, logp: _ToFloat0ND, /, *, method: _ICDFMethod | None = None) -> _FloatND: ...

    #
    @overload  # self: T1-d, p: 0-d
    def iccdf(
        self: _BaseDistribution[Any, _ShapeT1], p: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: 0-d
    def iccdf(self: _BaseDist0, p: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, p: 1-d
    def iccdf(self: _BaseDist0, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, p: 2-d
    def iccdf(self: _BaseDist0, p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, p: 3-d
    def iccdf(self: _BaseDist0, p: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, p: T1-d
    def iccdf(self: _BaseDist0, p: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, p: >=1-d
    def iccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, p: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, p: 1-d
    def iccdf(self: _BaseDist1, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, p: 2-d
    def iccdf(self: _BaseDist1, p: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, p: >=-d
    def iccdf(self: _BaseDist1, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, p: 1-d
    def iccdf(self: _BaseDist2, p: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, p: >=1-d
    def iccdf(self: _BaseDist2, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, p: >=1-d
    def iccdf(self: _BaseDist3, p: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def iccdf(self: _BaseDist1N, p: _ToFloat0ND, /, *, method: _ICDFMethod | None = None) -> _FloatND: ...
    #
    @overload  # self: T1-d, logp: 0-d
    def ilogccdf(
        self: _BaseDistribution[Any, _ShapeT1], logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, logp: 0-d
    def ilogccdf(self: _BaseDist0, logp: onp.ToFloat, /, *, method: _ICDFMethod | None = None) -> _Float: ...
    @overload  # self: 0-d, logp: 1-d
    def ilogccdf(self: _BaseDist0, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float1D: ...
    @overload  # self: 0-d, logp: 2-d
    def ilogccdf(self: _BaseDist0, logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 0-d, logp: 3-d
    def ilogccdf(self: _BaseDist0, logp: onp.ToFloatStrict3D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 0-d, logp: T1-d
    def ilogccdf(self: _BaseDist0, logp: _ToFloatND[_ShapeT1], /, *, method: _ICDFMethod | None = None) -> _FloatND[_ShapeT1]: ...
    @overload  # self: 0-d, q: >=1-d
    def ilogccdf(  # first union type is needed on `numpy<2.1`
        self: _BaseDist0, logp: _ToFloatND[_ShapeT1] | onp.ToFloatND, /, *, method: _ICDFMethod | None = None
    ) -> _FloatND[_ShapeT1] | _Float1ND: ...
    @overload  # self: 1-d, logp: 1-d
    def ilogccdf(self: _BaseDist1, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float2D: ...
    @overload  # self: 1-d, logp: 2-d
    def ilogccdf(self: _BaseDist1, logp: onp.ToFloatStrict2D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 1-d, logp: >=-d
    def ilogccdf(self: _BaseDist1, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float2ND: ...
    @overload  # self: 2-d, logp: 1-d
    def ilogccdf(self: _BaseDist2, logp: onp.ToFloatStrict1D, /, *, method: _ICDFMethod | None = None) -> _Float3D: ...
    @overload  # self: 2-d, logp: >=1-d
    def ilogccdf(self: _BaseDist2, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: 3-d, logp: >=1-d
    def ilogccdf(self: _BaseDist3, logp: onp.ToFloatND, /, *, method: _ICDFMethod | None = None) -> _Float3ND: ...
    @overload  # self: >=1-d
    def ilogccdf(self: _BaseDist1N, logp: _ToFloat0ND, /, *, method: _ICDFMethod | None = None) -> _FloatND: ...

#
class UnivariateDistribution(_BaseDistribution[_XT_co], Generic[_XT_co, _ShapeT0_co]):
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
    def __neg__(self, /) -> _LinDist[Self, _FloatT_co, _ShapeT_co]: ...
    def __abs__(self, /) -> _FoldDist[Self, _FloatT_co, _ShapeT_co]: ...

    #
    @overload
    def __add__(self, x: float | _Int | np.bool_, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __add__(self, x: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __add__(self, x: onp.ToFloat, /) -> _LinDist[Self, _Float, _ShapeT_co]: ...
    @overload
    def __add__(self: _DistT0, x: onp.CanArrayND[_FloatT, _ShapeT1], /) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __add__(self: _DistT_1, x: onp.ToFloatStrict1D, /) -> _LinDist[_DistT_1, _Float, _1D]: ...
    @overload
    def __add__(self: _DistT_2, x: onp.ToFloatStrict2D, /) -> _LinDist[_DistT_2, _Float, _2D]: ...
    @overload
    def __add__(self: _DistT_3, x: onp.ToFloatStrict3D, /) -> _LinDist[_DistT_3, _Float, _3D]: ...
    @overload
    def __add__(self, x: onp.ToFloatND, /) -> _LinDist[Self]: ...
    __radd__ = __add__

    #
    @overload
    def __sub__(self, lshift: float | _Int | np.bool_, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __sub__(self, lshift: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __sub__(self, lshift: onp.ToFloat, /) -> _LinDist[Self, _Float, _ShapeT_co]: ...
    @overload
    def __sub__(
        self: _DistT0, lshift: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __sub__(self: _DistT_1, lshift: onp.ToFloatStrict1D, /) -> _LinDist[_DistT_1, _Float, _1D]: ...
    @overload
    def __sub__(self: _DistT_2, lshift: onp.ToFloatStrict2D, /) -> _LinDist[_DistT_2, _Float, _2D]: ...
    @overload
    def __sub__(self: _DistT_3, lshift: onp.ToFloatStrict3D, /) -> _LinDist[_DistT_3, _Float, _3D]: ...
    @overload
    def __sub__(self, lshift: onp.ToFloatND, /) -> _LinDist[Self]: ...
    __rsub__ = __sub__

    #
    @overload
    def __mul__(self, scale: float | _Int | np.bool_, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __mul__(self, scale: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __mul__(self, scale: onp.ToFloat, /) -> _LinDist[Self, _Float, _ShapeT_co]: ...
    @overload
    def __mul__(
        self: _DistT0, scale: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __mul__(self: _DistT_1, scale: onp.ToFloatStrict1D, /) -> _LinDist[_DistT_1, _Float, _1D]: ...
    @overload
    def __mul__(self: _DistT_2, scale: onp.ToFloatStrict2D, /) -> _LinDist[_DistT_2, _Float, _2D]: ...
    @overload
    def __mul__(self: _DistT_3, scale: onp.ToFloatStrict3D, /) -> _LinDist[_DistT_3, _Float, _3D]: ...
    @overload
    def __mul__(self, scale: onp.ToFloatND, /) -> _LinDist[Self]: ...
    __rmul__ = __mul__

    #
    def __pow__(self, exp: onp.ToInt, /) -> MonotonicTransformedDistribution[Self, _ShapeT_co]: ...
    __rpow__ = __pow__

    #
    @overload
    def __truediv__(self, iscale: float | _Int | np.bool_, /) -> _LinDist[Self, np.float64 | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __truediv__(self, iscale: _FloatT, /) -> _LinDist[Self, _FloatT | _FloatT_co, _ShapeT_co]: ...
    @overload
    def __truediv__(self, iscale: onp.ToFloat, /) -> _LinDist[Self, _Float, _ShapeT_co]: ...
    @overload
    def __truediv__(
        self: _DistT0, iscale: onp.CanArrayND[_FloatT, _ShapeT1], /
    ) -> _LinDist[_DistT0, _FloatT | _FloatT_co, _ShapeT1]: ...
    @overload
    def __truediv__(self: _DistT_1, iscale: onp.ToFloatStrict1D, /) -> _LinDist[_DistT_1, _Float, _1D]: ...
    @overload
    def __truediv__(self: _DistT_2, iscale: onp.ToFloatStrict2D, /) -> _LinDist[_DistT_2, _Float, _2D]: ...
    @overload
    def __truediv__(self: _DistT_3, iscale: onp.ToFloatStrict3D, /) -> _LinDist[_DistT_3, _Float, _3D]: ...
    @overload
    def __truediv__(self, iscale: onp.ToFloatND, /) -> _LinDist[Self]: ...
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

class ShiftedScaledDistribution(_TransDist[_DistT_co, _FloatT_co, _ShapeT_co], Generic[_DistT_co, _FloatT_co, _ShapeT_co]):
    _loc_domain: ClassVar[_RealInterval] = ...
    _loc_param: ClassVar[_RealParameter] = ...
    _scale_domain: ClassVar[_RealInterval] = ...
    _scale_param: ClassVar[_RealParameter] = ...

    loc: _ParamField[_FloatT_co, _ShapeT_co]
    scale: _ParamField[_FloatT_co, _ShapeT_co]

class FoldedDistribution(_TransDist[_DistT_co, _FloatT_co, _ShapeT_co], Generic[_DistT_co, _FloatT_co, _ShapeT_co]):
    @overload
    def __init__(self: _FoldDist[_DistT0, _Float, _0D], X: _DistT0, /, *args: Never, **kwargs: Unpack[_DistOpts]) -> None: ...
    @overload
    def __init__(self: _FoldDist[_DistT1, _Float, _1D], X: _DistT1, /, *args: Never, **kwargs: Unpack[_DistOpts]) -> None: ...
    @overload
    def __init__(self: _FoldDist[_DistT2, _Float, _2D], X: _DistT2, /, *args: Never, **kwargs: Unpack[_DistOpts]) -> None: ...
    @overload
    def __init__(self: _FoldDist[_DistT3, _Float, _3D], X: _DistT3, /, *args: Never, **kwargs: Unpack[_DistOpts]) -> None: ...
    @overload
    def __init__(self: _FoldDist[_DistT, _Float, _ND], X: _DistT, /, *args: Never, **kwargs: Unpack[_DistOpts]) -> None: ...

class TruncatedDistribution(_TransDist[_DistT_co, _Float, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]):
    _lb_domain: ClassVar[_RealInterval] = ...
    _lb_param: ClassVar[_RealParameter] = ...
    _ub_domain: ClassVar[_RealInterval] = ...
    _ub_param: ClassVar[_RealParameter] = ...

    lb: _ParamField[_Float, _ShapeT_co]
    ub: _ParamField[_Float, _ShapeT_co]

    @overload
    def __init__(
        self: _TruncDist[_DistT0, _0D],
        X: _DistT0,
        /,
        *args: Never,
        lb: onp.ToFloat = ...,
        ub: onp.ToFloat = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: _TruncDist[_DistT1, _1D],
        X: _DistT1,
        /,
        *args: Never,
        lb: _ToFloatMax1D = ...,
        ub: _ToFloatMax1D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: _TruncDist[_DistT2, _2D],
        X: _DistT2,
        /,
        *args: Never,
        lb: _ToFloatMax2D = ...,
        ub: _ToFloatMax2D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: _TruncDist[_DistT3, _3D],
        X: _DistT3,
        /,
        *args: Never,
        lb: _ToFloatMax3D = ...,
        ub: _ToFloatMax3D = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: _TruncDist[_DistT, _ND],
        X: _DistT,
        /,
        *args: Never,
        lb: _ToFloat0ND = ...,
        ub: _ToFloat0ND = ...,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

# always float64 or longdouble
class OrderStatisticDistribution(_TransDist[_DistT_co, _OutFloat, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]):
    # these should actually be integral; but the `_IntegerDomain` isn't finished yet
    _r_domain: ClassVar[_RealInterval] = ...
    _r_param: ClassVar[_RealParameter] = ...
    _n_domain: ClassVar[_RealInterval] = ...
    _n_param: ClassVar[_RealParameter] = ...

    @overload
    def __init__(
        self: OrderStatisticDistribution[_DistT0, _0D],
        dist: _DistT0,
        /,
        *args: Never,
        r: onp.ToJustInt,
        n: onp.ToJustInt,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: OrderStatisticDistribution[_DistT1, _1D],
        dist: _DistT1,
        /,
        *args: Never,
        r: _ToJustIntMax1D,
        n: _ToJustIntMax1D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: OrderStatisticDistribution[_DistT2, _2D],
        dist: _DistT2,
        /,
        *args: Never,
        r: _ToJustIntMax2D,
        n: _ToJustIntMax2D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: OrderStatisticDistribution[_DistT3, _3D],
        dist: _DistT3,
        /,
        *args: Never,
        r: _ToJustIntMax3D,
        n: _ToJustIntMax3D,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...
    @overload
    def __init__(
        self: OrderStatisticDistribution[_DistT, _ND],
        X: _DistT,
        /,
        *args: Never,
        r: _ToJustIntMaxND,
        n: _ToJustIntMaxND,
        **kwargs: Unpack[_DistOpts],
    ) -> None: ...

# without HKT there's no reasonable way tot determine the floating scalar type
class MonotonicTransformedDistribution(_TransDist[_DistT_co, _Float, _ShapeT_co], Generic[_DistT_co, _ShapeT_co]):
    _g: Final[_Elementwise]
    _h: Final[_Elementwise]
    _dh: Final[_Elementwise]
    _logdh: Final[_Elementwise]
    _increasing: Final[bool]
    _repr_pattern: Final[str]
    _str_pattern: Final[str]

    def __init__(
        self: MonotonicTransformedDistribution[_CDist[_ShapeT], _ShapeT],
        X: _DistT_co,
        /,
        *args: Never,
        g: _Elementwise,
        h: _Elementwise,
        dh: _Elementwise,
        logdh: _Elementwise | None = None,
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
    def kurtosis(self, /, *, method: _SMomentMethod | None = None) -> _OutFloat: ...  # pyright: ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]

###

# still waiting on the intersection type PEP...

@overload
def truncate(X: _DistT0, lb: onp.ToFloat = ..., ub: onp.ToFloat = ...) -> _TruncDist[_DistT0, _0D]: ...
@overload
def truncate(X: _DistT1, lb: _ToFloatMax1D = ..., ub: _ToFloatMax1D = ...) -> _TruncDist[_DistT1, _1D]: ...
@overload
def truncate(X: _DistT2, lb: _ToFloatMax2D = ..., ub: _ToFloatMax2D = ...) -> _TruncDist[_DistT2, _2D]: ...
@overload
def truncate(X: _DistT3, lb: _ToFloatMax3D = ..., ub: _ToFloatMax3D = ...) -> _TruncDist[_DistT3, _3D]: ...
@overload
def truncate(X: _DistT, lb: _ToFloat0ND = ..., ub: _ToFloat0ND = ...) -> _TruncDist[_DistT, _ND]: ...

#
@overload
def order_statistic(X: _DistT0, /, *, r: onp.ToJustInt, n: onp.ToJustInt) -> OrderStatisticDistribution[_DistT0, _0D]: ...
@overload
def order_statistic(X: _DistT1, /, *, r: _ToJustIntMax1D, n: _ToJustIntMax1D) -> OrderStatisticDistribution[_DistT1, _1D]: ...
@overload
def order_statistic(X: _DistT2, /, *, r: _ToJustIntMax2D, n: _ToJustIntMax2D) -> OrderStatisticDistribution[_DistT2, _2D]: ...
@overload
def order_statistic(X: _DistT3, /, *, r: _ToJustIntMax3D, n: _ToJustIntMax3D) -> OrderStatisticDistribution[_DistT3, _3D]: ...
@overload
def order_statistic(X: _DistT, /, *, r: _ToJustIntMaxND, n: _ToJustIntMaxND) -> OrderStatisticDistribution[_DistT, _ND]: ...

#
@overload
def abs(X: _DistT0, /) -> _FoldDist[_DistT0, _Float, _0D]: ...
@overload
def abs(X: _DistT1, /) -> _FoldDist[_DistT1, _Float, _1D]: ...
@overload
def abs(X: _DistT2, /) -> _FoldDist[_DistT2, _Float, _2D]: ...
@overload
def abs(X: _DistT3, /) -> _FoldDist[_DistT3, _Float, _3D]: ...
@overload
def abs(X: _DistT, /) -> _FoldDist[_DistT, _Float, _ND]: ...

#
@overload
def exp(X: _DistT0, /) -> MonotonicTransformedDistribution[_DistT0, _0D]: ...
@overload
def exp(X: _DistT1, /) -> MonotonicTransformedDistribution[_DistT1, _1D]: ...
@overload
def exp(X: _DistT2, /) -> MonotonicTransformedDistribution[_DistT2, _2D]: ...
@overload
def exp(X: _DistT3, /) -> MonotonicTransformedDistribution[_DistT3, _3D]: ...
@overload
def exp(X: _DistT, /) -> MonotonicTransformedDistribution[_DistT, _ND]: ...

#
@overload
def log(X: _DistT0, /) -> MonotonicTransformedDistribution[_DistT0, _0D]: ...
@overload
def log(X: _DistT1, /) -> MonotonicTransformedDistribution[_DistT1, _1D]: ...
@overload
def log(X: _DistT2, /) -> MonotonicTransformedDistribution[_DistT2, _2D]: ...
@overload
def log(X: _DistT3, /) -> MonotonicTransformedDistribution[_DistT3, _3D]: ...
@overload
def log(X: _DistT, /) -> MonotonicTransformedDistribution[_DistT, _ND]: ...

# NOTE: These currently don't support >0-d parameters, and it looks like they always return float64, regardless of dtype
@type_check_only
class CustomContinuousDistribution(ContinuousDistribution[np.float64, _0D]):
    _dtype: np.dtype[_Float]  # ignored

# TODO(jorenham): support for `rv_discrete`
def make_distribution(dist: rv_continuous | _DuckDistributionType) -> type[CustomContinuousDistribution]: ...
