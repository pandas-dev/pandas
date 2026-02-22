# NOTE: this is needed because of the >50 LSP violations...
# mypy: disable-error-code="override"
# pyright: reportIncompatibleMethodOverride = false
import types
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, Final, Generic, Literal as L, Self, SupportsIndex, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, Unpack, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape
from scipy.integrate._typing import QuadOpts as _QuadOpts
from scipy.stats._censored_data import CensoredData

_T = TypeVar("_T")
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[int, ...])
_ArgT = TypeVar("_ArgT", bound=_ToFloatOrND, default=_ToFloatOrND)
_FloatNDT = TypeVar("_FloatNDT", bound=_FloatOrND, default=_FloatOrND)

_FloatNDT_co = TypeVar("_FloatNDT_co", bound=_FloatOrND, default=_FloatOrND, covariant=True)
_RVT = TypeVar("_RVT", bound=rv_generic, default=rv_generic)
_RVT_co = TypeVar("_RVT_co", bound=rv_generic, default=rv_generic, covariant=True)
_CRVT_co = TypeVar("_CRVT_co", bound=rv_continuous, default=rv_continuous, covariant=True)
_DRVT_co = TypeVar("_DRVT_co", bound=rv_discrete, default=rv_discrete, covariant=True)
_XKT_co = TypeVar("_XKT_co", bound=_CoFloat, covariant=True, default=_CoFloat)
_PKT_co = TypeVar("_PKT_co", bound=_Floating, covariant=True, default=_Floating)

_Tuple2: TypeAlias = tuple[_T, _T]
_Tuple3: TypeAlias = tuple[_T, _T, _T]
_Tuple4: TypeAlias = tuple[_T, _T, _T, _T]

_Floating: TypeAlias = np.float64 | np.float32 | np.float16  # longdouble often results in trouble
_CoFloat: TypeAlias = _Floating | npc.integer

_Bool: TypeAlias = bool | np.bool_
_Int: TypeAlias = int | np.int32 | np.int64
_Float: TypeAlias = float | np.float64

_Float0D: TypeAlias = onp.Array0D[np.float64]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_BoolND: TypeAlias = onp.ArrayND[np.bool_]
_IntND: TypeAlias = onp.ArrayND[np.int_]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_CoFloatND: TypeAlias = onp.ArrayND[_CoFloat]

_BoolOrND: TypeAlias = _Bool | _BoolND
_IntOrND: TypeAlias = _Int | _IntND
_FloatOrND: TypeAlias = _Float | _FloatND

# pyright bug workaround on `numpy<2.1` (note the weird shape-type)
_Float1ND: TypeAlias = onp.ArrayND[np.float64, tuple[int] | tuple[Any, ...]]
_FloatOr1ND: TypeAlias = _Float | _Float1ND

_ToFloatOrND: TypeAlias = onp.ToFloat | onp.ToFloatND

_Expectant: TypeAlias = Callable[[float], onp.ToFloat]

# there are at most 4 + 2 args
_RVArgs: TypeAlias = (
    tuple[()]
    | tuple[_ArgT]
    | tuple[_ArgT, _ArgT]
    | tuple[_ArgT, _ArgT, _ArgT]
    | tuple[_ArgT, _ArgT, _ArgT, _ArgT]
    | tuple[_ArgT, _ArgT, _ArgT, _ArgT, _ArgT]
    | tuple[_ArgT, _ArgT, _ArgT, _ArgT, _ArgT, _ArgT]
)
_RVKwds: TypeAlias = dict[str, _ToFloatOrND]

_Moment1: TypeAlias = L["m", "v", "s", "k"]
_Moment2: TypeAlias = L[
    "mv", "ms", "mk",
    "vm", "vs", "vk",
    "sm", "sv", "sk",
    "km", "kv", "ks",
]  # fmt: skip
_Moment3: TypeAlias = L[
    "mvs", "mvk", "msv", "msk", "mkv", "mks",
    "vms", "vmk", "vsm", "vsk", "vkm", "vks",
    "smv", "smk", "svm", "svk", "skm", "skv",
    "kmv", "kms", "kvm", "kvs", "ksm", "ksv",
]  # fmt: skip
_Moment4: TypeAlias = L[
    "mvsk", "mvks", "msvk", "mskv", "mkvs", "mksv",
    "vmsk", "vmks", "vsmk", "vskm", "vkms", "vksm",
    "smvk", "smkv", "svmk", "svkm", "skmv", "skvm",
    "kmvs", "kmsv", "kvms", "kvsm", "ksmv", "ksvm",
]  # fmt: skip

_MomentType: TypeAlias = L[0, 1]
_FitMethod: TypeAlias = L["MLE", "MM"]

###

docheaders: Final[dict[str, str]] = ...
docdict: Final[dict[str, str]] = ...
docdict_discrete: Final[dict[str, str]] = ...
parse_arg_template: Final[str] = ...

class rv_frozen(Generic[_RVT_co, _FloatNDT_co]):
    dist: _RVT_co
    args: _RVArgs[_FloatNDT_co]
    kwds: _RVKwds

    @property
    def random_state(self, /) -> onp.random.RNG: ...
    @random_state.setter
    def random_state(self, seed: onp.random.ToRNG | None, /) -> None: ...

    #
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    #
    @overload
    def __init__(self: rv_frozen[_RVT, _Float], /, dist: _RVT) -> None: ...
    @overload
    def __init__(self, /, dist: _RVT_co, *args: _FloatNDT_co, **kwds: _FloatNDT_co) -> None: ...
    @overload
    def __init__(self, /, dist: _RVT_co, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> None: ...

    #
    @overload
    def cdf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def cdf(self, /, x: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def logcdf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def logcdf(self, /, x: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def sf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def sf(self, /, x: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def logsf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def logsf(self, /, x: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def ppf(self, /, q: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def ppf(self, /, q: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def isf(self, /, q: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def isf(self, /, q: onp.ToFloatND) -> _FloatND: ...

    #
    @overload
    def rvs(self, /, size: tuple[()] | None = None, random_state: onp.random.ToRNG | None = None) -> _FloatNDT_co: ...
    @overload
    def rvs(
        self, /, size: op.CanIndex | tuple[op.CanIndex, *tuple[op.CanIndex, ...]], random_state: onp.random.ToRNG | None = None
    ) -> _FloatND: ...
    @overload
    def rvs(self, /, size: AnyShape | None = None, random_state: onp.random.ToRNG | None = None) -> _FloatOrND: ...

    #
    def median(self, /) -> _FloatNDT_co: ...
    def mean(self, /) -> _FloatNDT_co: ...
    def var(self, /) -> _FloatNDT_co: ...
    def std(self, /) -> _FloatNDT_co: ...
    def entropy(self, /) -> _FloatNDT_co: ...

    #
    @overload
    def stats(self, /, moments: _Moment1) -> _FloatNDT_co: ...
    @overload
    def stats(self, /, moments: _Moment2 = "mv") -> _Tuple2[_FloatNDT_co]: ...
    @overload
    def stats(self, /, moments: _Moment3) -> _Tuple3[_FloatNDT_co]: ...
    @overload
    def stats(self, /, moments: _Moment4) -> _Tuple4[_FloatNDT_co]: ...

    # NOTE: Even though `order` defaults to `None`, this will raise  a `TypeError`.
    def moment(self, /, order: onp.ToInt | None = None) -> _FloatNDT_co: ...

    # NOTE: Will raise a `TypeError` with n-D parameters.
    def expect(
        self: rv_frozen[_RVT, _Float],
        /,
        func: _Expectant | None = None,
        lb: onp.ToFloat | None = None,
        ub: onp.ToFloat | None = None,
        conditional: _Bool = False,
        **kwds: Unpack[_QuadOpts],
    ) -> _Float: ...

    #
    def support(self, /) -> _Tuple2[_FloatNDT_co]: ...
    def interval(self, /, confidence: onp.ToFloat | None = None) -> _Tuple2[_FloatNDT_co]: ...

# undocumented
class rv_continuous_frozen(rv_frozen[_CRVT_co, _FloatNDT_co], Generic[_CRVT_co, _FloatNDT_co]):
    @overload
    def pdf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def pdf(self, /, x: onp.ToFloatND) -> _FloatND: ...
    #
    @overload
    def logpdf(self, /, x: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def logpdf(self, /, x: onp.ToFloatND) -> _FloatND: ...

# undocumented
class rv_discrete_frozen(rv_frozen[_DRVT_co, _FloatNDT_co], Generic[_DRVT_co, _FloatNDT_co]):
    @overload
    def pmf(self, /, k: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def pmf(self, /, k: onp.ToFloatND) -> _FloatND: ...
    #
    @overload
    def logpmf(self, /, k: onp.ToFloat) -> _FloatNDT_co: ...
    @overload
    def logpmf(self, /, k: onp.ToFloatND) -> _FloatND: ...

# NOTE: Because of the limitations of `ParamSpec`, there is no proper way to annotate specific "positional or keyword arguments".
# Considering the Liskov Substitution Principle, the only remaining option is to annotate `*args, and `**kwargs` as `Any`.
class rv_generic:
    @property
    def random_state(self, /) -> onp.random.RNG: ...
    @random_state.setter
    def random_state(self, seed: onp.random.ToRNG | None, /) -> None: ...

    #
    def __init__(self, /, seed: onp.random.ToRNG | None = None) -> None: ...
    def _attach_methods(self, /) -> None: ...
    def _attach_argparser_methods(self, /) -> None: ...
    def _construct_argparser(
        self, /, meths_to_inspect: Iterable[Callable[..., Any]], locscale_in: str, locscale_out: str
    ) -> None: ...
    def _construct_doc(self, /, docdict: dict[str, str], shapes_vals: tuple[float, ...] | None = None) -> None: ...
    def _construct_default_doc(
        self,
        /,
        longname: str | None = None,
        docdict: dict[str, str] | None = None,
        discrete: L["continuous", "discrete"] = "continuous",
    ) -> None: ...

    #
    @overload
    def __call__(self, /) -> rv_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> rv_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> rv_frozen[Self]: ...

    #
    @overload
    def freeze(self, /) -> rv_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> rv_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> rv_frozen[Self]: ...
    #
    def _stats(self, /, *args: onp.ToFloat, **kwds: object) -> _Tuple4[_Float | None] | _Tuple4[_FloatND | None]: ...
    def _munp(self, /, n: onp.ToInt | onp.ToIntND, *args: onp.ToFloat) -> _FloatND: ...

    #
    def _argcheck_rvs(
        self, /, *args: onp.ToFloat, size: onp.ToInt | onp.ToIntND | None = None
    ) -> tuple[list[_CoFloatND], _CoFloatND, _CoFloatND, tuple[int, ...] | tuple[np.int_, ...]]: ...
    def _argcheck(self, /, *args: onp.ToFloat) -> _BoolOrND: ...

    #
    def _get_support(self, /, *args: onp.ToFloat, **kwargs: onp.ToFloat) -> _Tuple2[_FloatOrND]: ...
    def _support_mask(self, /, x: _CoFloatND, *args: onp.ToFloat) -> _BoolND: ...
    def _open_support_mask(self, /, x: _CoFloatND, *args: onp.ToFloat) -> _BoolOrND: ...

    #
    def _rvs(
        self, /, *args: onp.ToFloat, size: AnyShape | None = None, random_state: onp.random.ToRNG | None = None
    ) -> _FloatOrND: ...

    #
    def _logcdf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _sf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _logsf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _ppf(self, /, q: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _isf(self, /, q: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    #
    @overload
    def rvs(
        self, /, *args: onp.ToFloat, random_state: onp.random.ToRNG | None, discrete: onp.ToTrue, **kwds: onp.ToFloat
    ) -> _Int: ...
    @overload
    def rvs(
        self, /, *args: onp.ToFloat, random_state: onp.random.ToRNG | None, discrete: onp.ToTrue, **kwds: _ToFloatOrND
    ) -> _IntOrND: ...
    @overload
    def rvs(
        self,
        /,
        *args: onp.ToFloat,
        random_state: onp.random.ToRNG | None,
        discrete: onp.ToFalse | None = ...,
        **kwds: onp.ToFloat,
    ) -> _Float: ...
    @overload
    def rvs(
        self,
        /,
        *args: _ToFloatOrND,
        random_state: onp.random.ToRNG | None,
        discrete: onp.ToFalse | None = ...,
        **kwds: _ToFloatOrND,
    ) -> _FloatOrND: ...

    #
    @overload
    def stats(self, /, *args: onp.ToFloat, moment: _Moment1, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def stats(self, /, *args: onp.ToFloat, moment: _Moment2 = ..., **kwds: onp.ToFloat) -> _Tuple2[_Float]: ...
    @overload
    def stats(self, /, *args: onp.ToFloat, moment: _Moment3, **kwds: onp.ToFloat) -> _Tuple3[_Float]: ...
    @overload
    def stats(self, /, *args: onp.ToFloat, moment: _Moment4, **kwds: onp.ToFloat) -> _Tuple4[_Float]: ...
    @overload
    def stats(self, /, *args: _ToFloatOrND, moment: _Moment1, **kwds: _ToFloatOrND) -> _FloatOrND: ...
    @overload
    def stats(self, /, *args: _ToFloatOrND, moment: _Moment2 = ..., **kwds: _ToFloatOrND) -> _Tuple2[_FloatOrND]: ...
    @overload
    def stats(self, /, *args: _ToFloatOrND, moment: _Moment3, **kwds: _ToFloatOrND) -> _Tuple3[_FloatOrND]: ...
    @overload
    def stats(self, /, *args: _ToFloatOrND, moment: _Moment4, **kwds: _ToFloatOrND) -> _Tuple4[_FloatOrND]: ...

    #
    @overload
    def entropy(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def entropy(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def moment(self, /, order: onp.ToInt, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def moment(self, /, order: onp.ToInt, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def median(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def median(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def mean(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def mean(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def var(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def var(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def std(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def std(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _FloatOrND: ...

    #
    @overload
    def interval(self, /, confidence: onp.ToFloat, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Tuple2[_Float]: ...
    @overload
    def interval(self, /, confidence: _ToFloatOrND, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _Tuple2[_FloatOrND]: ...

    #
    @overload
    def support(self, /, *args: onp.ToFloat, **kwds: onp.ToFloat) -> _Tuple2[_Float]: ...
    @overload
    def support(self, /, *args: _ToFloatOrND, **kwds: _ToFloatOrND) -> _Tuple2[_FloatOrND]: ...

    #
    @overload
    def nnlf(self, /, theta: onp.ToFloat1D, x: onp.ToFloatStrict1D) -> _Float | _Float0D: ...
    @overload
    def nnlf(self, /, theta: onp.ToFloat1D, x: onp.ToFloatStrict2D) -> _Float1D: ...
    @overload
    def nnlf(self, /, theta: onp.ToFloat1D, x: onp.ToFloatStrict3D) -> _Float2D: ...
    @overload
    def nnlf(self, /, theta: onp.ToFloat1D, x: onp.ToFloatND) -> _FloatOrND: ...

    #
    def _nnlf(self, /, x: _CoFloatND, *args: onp.ToFloat) -> _Float | _Float0D: ...
    def _penalized_nnlf(self, /, theta: Sequence[Any], x: _CoFloatND) -> _Float | _Float0D: ...
    def _penalized_nlpsf(self, /, theta: Sequence[Any], x: _CoFloatND) -> _Float | _Float0D: ...

class _ShapeInfo:
    name: Final[str]
    integrality: Final[bool]
    domain: Final[Sequence[_Float]]  # in practice always a list of size two

    def __init__(
        self, /, name: str, integrality: bool = False, domain: Sequence[_Float] = ..., inclusive: Sequence[bool] = (True, True)
    ) -> None: ...

@type_check_only
class _rv_mixin:
    name: Final[str]
    a: Final[_Float]
    b: Final[_Float]
    badvalue: Final[_Float]
    shapes: Final[str]

    def generic_moment(self, /, n: onp.ToInt | onp.ToIntND, *args: onp.ToFloat) -> _FloatND: ...

    #
    def _shape_info(self, /) -> list[_ShapeInfo]: ...
    def _param_info(self, /) -> list[_ShapeInfo]: ...
    def _attach_methods(self, /) -> None: ...
    def _logpxf(self, /, x: _CoFloatND, *args: onp.ToFloat) -> _FloatND: ...
    def _cdf_single(self, /, x: onp.ToFloat, *args: onp.ToFloat) -> _Float: ...
    def _cdfvec(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _cdf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _ppfvec(self, /, q: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    @overload
    def _unpack_loc_scale(self, /, theta: Sequence[onp.ToFloat]) -> tuple[onp.ToFloat, onp.ToFloat, tuple[onp.ToFloat, ...]]: ...
    @overload
    def _unpack_loc_scale(
        self, /, theta: Sequence[onp.ToFloatND]
    ) -> tuple[onp.ToFloatND, onp.ToFloatND, tuple[onp.ToFloatND, ...]]: ...

class rv_continuous(_rv_mixin, rv_generic):
    moment_type: Final[_MomentType]
    xtol: Final[_Float]

    def __init__(
        self,
        /,
        momtype: _MomentType = 1,
        a: _Float | None = None,
        b: _Float | None = None,
        xtol: _Float = 1e-14,
        badvalue: _Float | None = None,
        name: str | None = None,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...

    #
    @overload
    def __call__(self, /) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def __call__(
        self, /, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def __call__(
        self, /, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> rv_continuous_frozen[Self]: ...

    #
    @overload
    def freeze(self, /) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def freeze(
        self, /, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def freeze(
        self, /, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> rv_continuous_frozen[Self]: ...

    #
    def _pdf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...
    def _logpdf(self, /, x: _FloatNDT, *args: onp.ToFloat) -> _FloatNDT: ...

    #
    @overload
    def pdf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def pdf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def pdf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def logpdf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def logpdf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logpdf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def cdf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def cdf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def logcdf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logcdf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def sf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def sf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def sf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def logsf(
        self, /, x: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def logsf(
        self,
        /,
        x: onp.CanArrayND[_CoFloat, _ShapeT],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logsf(
        self, /, x: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def ppf(
        self, /, q: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def ppf(
        self,
        /,
        q: onp.Array[_ShapeT, _CoFloat],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def ppf(
        self, /, q: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    @overload
    def isf(
        self, /, q: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, **kwds: onp.ToFloat
    ) -> _Float: ...
    @overload
    def isf(
        self,
        /,
        q: onp.Array[_ShapeT, _CoFloat],
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        **kwds: onp.ToFloat,
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def isf(
        self, /, q: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1, **kwds: _ToFloatOrND
    ) -> _FloatOr1ND: ...

    #
    def _nnlf_and_penalty(self, /, x: _FloatND, args: Sequence[onp.ToFloat]) -> _Float: ...

    #
    def _reduce_func(
        self, /, args: tuple[onp.ToFloat, ...], kwds: Mapping[str, onp.ToFloat], data: _ToFloatOrND | None = None
    ) -> tuple[
        list[_Float],
        Callable[[list[onp.ToFloat], _CoFloatND], _Float],
        Callable[[list[onp.ToFloat], _CoFloatND], list[_Float]],
        list[_Float],
    ]: ...

    #
    def _moment_error(self, /, theta: list[onp.ToFloat], x: _CoFloatND, data_moments: onp.ToFloat1D) -> _Float: ...

    #
    def _fitstart(
        self, /, data: _FloatND, args: tuple[onp.ToFloat, ...] | None = None
    ) -> tuple[*tuple[_Float, ...], _Float, _Float]: ...

    #
    def _fit_loc_scale_support(self, /, data: _ToFloatOrND, *args: onp.ToFloat) -> _Tuple2[np.int32 | np.int64 | _Float]: ...
    def fit_loc_scale(self, /, data: _ToFloatOrND, *args: onp.ToFloat) -> _Tuple2[_Float]: ...

    #
    def fit(
        self,
        /,
        data: _ToFloatOrND | CensoredData[np.float64],
        *args: onp.ToFloat,
        optimizer: Callable[[_FloatND, tuple[float, ...], tuple[float, ...], bool], tuple[onp.ToFloat, ...]] | None = ...,
        method: _FitMethod = "MLE",
        **kwds: onp.ToFloat,
    ) -> tuple[_Float, ...]: ...

    #
    def expect(
        self,
        /,
        func: _Expectant | None = None,
        args: tuple[onp.ToFloat, ...] = (),
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        lb: onp.ToFloat | None = None,
        ub: onp.ToFloat | None = None,
        conditional: op.CanBool = False,
        **kwds: Unpack[_QuadOpts],
    ) -> _Float: ...

    #
    @override
    @overload  # size: int | (int, )
    def rvs(
        self,
        /,
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        size: SupportsIndex | tuple[SupportsIndex],
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat,
    ) -> _Float1D: ...
    @overload  # size: (int, int)
    def rvs(
        self,
        /,
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        size: tuple[SupportsIndex, SupportsIndex],
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat,
    ) -> _Float2D: ...
    @overload  # size: ()  (default)
    def rvs(
        self,
        /,
        *args: onp.ToFloat,
        loc: onp.ToFloat = 0,
        scale: onp.ToFloat = 1,
        size: tuple[()] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat,
    ) -> np.float64: ...
    @overload  # size: int | (int, [int, ...])
    def rvs(
        self,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, *tuple[SupportsIndex, ...]],
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # arg0: Nd
    def rvs(
        self,
        arg0: onp.ToFloatND,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # arg1: Nd
    def rvs(
        self,
        arg0: onp.ToFloat | onp.ToFloatND,
        arg1: onp.ToFloatND,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # arg2: Nd
    def rvs(
        self,
        arg0: onp.ToFloat | onp.ToFloatND,
        arg1: onp.ToFloat | onp.ToFloatND,
        arg2: onp.ToFloatND,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # loc: Nd
    def rvs(
        self,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloatND,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # scale: Nd
    def rvs(
        self,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloatND,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloat | onp.ToFloatND,
    ) -> _FloatND: ...
    @overload  # kwargs: Nd  (technically kwargs can still be omitted, but it seems to work either way)
    def rvs(
        self,
        /,
        *args: onp.ToFloat | onp.ToFloatND,
        loc: onp.ToFloat | onp.ToFloatND = 0,
        scale: onp.ToFloat | onp.ToFloatND = 1,
        size: SupportsIndex | tuple[SupportsIndex, ...] = (),
        random_state: onp.random.ToRNG | None = None,
        **kwds: onp.ToFloatND,
    ) -> _FloatND: ...

class rv_discrete(_rv_mixin, rv_generic):
    inc: Final[int]
    moment_tol: Final[float]

    @overload
    def __new__(
        cls,
        a: onp.ToFloat = 0,
        b: onp.ToFloat = ...,
        name: str | None = None,
        badvalue: _Float | None = None,
        moment_tol: _Float = 1e-8,
        values: None = None,
        inc: int | np.int_ = 1,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> Self: ...
    # NOTE: The return types of the following overloads is ignored by mypy
    @overload
    def __new__(
        cls,
        a: onp.ToFloat,
        b: onp.ToFloat,
        name: str | None,
        badvalue: _Float | None,
        moment_tol: _Float,
        values: _Tuple2[onp.ToFloatND],
        inc: int | np.int_ = 1,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> rv_sample: ...
    @overload
    def __new__(
        cls,
        a: onp.ToFloat = 0,
        b: onp.ToFloat = ...,
        name: str | None = None,
        badvalue: _Float | None = None,
        moment_tol: _Float = 1e-8,
        *,
        values: _Tuple2[onp.ToFloatND],
        inc: int | np.int_ = 1,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> rv_sample: ...

    #
    def __init__(
        self,
        /,
        a: onp.ToFloat = 0,
        b: onp.ToFloat = ...,
        name: str | None = None,
        badvalue: _Float | None = None,
        moment_tol: _Float = 1e-8,
        # mypy workaround: `values` can only be None
        values: _Tuple2[onp.ToFloatND] | None = None,
        inc: int | np.int_ = 1,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...

    #
    @overload
    def __call__(self, /) -> rv_discrete_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> rv_discrete_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> rv_discrete_frozen[Self]: ...

    #
    @overload
    def freeze(self, /) -> rv_discrete_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> rv_discrete_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> rv_discrete_frozen[Self]: ...

    #
    @overload
    def pmf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def pmf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def pmf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def logpmf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def logpmf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logpmf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def cdf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def cdf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def cdf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def logcdf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def logcdf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logcdf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def sf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def sf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def sf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def logsf(self, /, k: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def logsf(
        self, /, k: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logsf(self, /, k: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def ppf(self, /, q: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def ppf(
        self, /, q: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def ppf(self, /, q: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    @overload
    def isf(self, /, q: onp.ToFloat, *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat) -> _Float: ...
    @overload
    def isf(
        self, /, q: onp.Array[_ShapeT, _CoFloat], *args: onp.ToFloat, loc: onp.ToFloat = 0, **kwds: onp.ToFloat
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def isf(self, /, q: _ToFloatOrND, *args: _ToFloatOrND, loc: _ToFloatOrND = 0, **kwds: _ToFloatOrND) -> _FloatOr1ND: ...

    #
    def expect(
        self,
        /,
        func: Callable[[onp.Array1D[np.int_]], onp.ToFloatND] | None = None,
        args: tuple[onp.ToFloat, ...] = (),
        loc: onp.ToFloat = 0,
        lb: onp.ToInt | None = None,
        ub: onp.ToInt | None = None,
        conditional: op.CanBool = False,
        maxcount: onp.ToInt = 1000,
        tolerance: onp.ToFloat = 1e-10,
        chunksize: onp.ToInt = 32,
    ) -> _Float: ...

    #
    @override
    def rvs(
        self,
        /,
        *args: _ToFloatOrND,
        loc: _ToFloatOrND = 0,
        size: AnyShape = 1,
        random_state: onp.random.ToRNG | None = None,
        **kwds: _ToFloatOrND,
    ) -> _IntOrND: ...

# returned by `rv_discrete.__new__` if `values` is specified
class rv_sample(rv_discrete, Generic[_XKT_co, _PKT_co]):
    xk: onp.Array1D[_XKT_co]
    pk: onp.Array1D[_PKT_co]
    qvals: onp.Array1D[_PKT_co]

    def __init__(
        self,
        /,
        a: onp.ToFloat = 0,
        b: onp.ToFloat = ...,
        name: str | None = None,
        badvalue: _Float | None = None,
        moment_tol: _Float = 1e-8,
        # never None in practice, but required by stubtest
        values: _Tuple2[onp.ToFloatND] | None = None,
        inc: int | np.int_ = 1,
        longname: str | None = None,
        shapes: str | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...

    #
    def _entropy(self, /) -> _Float: ...
    vecentropy: Final = _entropy

    #
    @override
    def generic_moment(self, /, n: onp.ToInt | onp.ToIntND) -> _FloatND: ...

# private helper subtypes
@type_check_only
class _rv_continuous_0(rv_continuous):
    # overrides of rv_generic
    @overload  # loc: 0-d, scale: 0-d, moments: 1 (positional)
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloat, moment: _Moment1) -> _Float: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 1 (keyword)
    def stats(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, *, moment: _Moment1) -> _Float: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 2 (default)
    def stats(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, moment: _Moment2 = "mv") -> _Tuple2[_Float]: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 3 (positional)
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloat, moment: _Moment3) -> _Tuple3[_Float]: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 3 (keyword)
    def stats(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, *, moment: _Moment3) -> _Tuple3[_Float]: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 4 (positional)
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloat, moment: _Moment4) -> _Tuple4[_Float]: ...
    @overload  # loc: 0-d, scale: 0-d, moments: 4 (keyword)
    def stats(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, *, moment: _Moment4) -> _Tuple4[_Float]: ...

    #
    @overload  # loc: 0-d, scale: n-d (positional), moments: 1
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloatND, moment: _Moment1) -> _FloatND: ...
    @overload  # loc: 0-d, scale: n-d (positional), moments: 2 (default)
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloatND, moment: _Moment2 = "mv") -> _Tuple2[_FloatND]: ...
    @overload  # loc: 0-d, scale: n-d (positional), moments: 3
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloatND, moment: _Moment3) -> _Tuple3[_FloatND]: ...
    @overload  # loc: 0-d, scale: n-d (positional), moments: 4
    def stats(self, /, loc: onp.ToFloat, scale: onp.ToFloatND, moment: _Moment4) -> _Tuple4[_FloatND]: ...

    #
    @overload  # loc: 0-d, scale: n-d (keyword), moments: 1
    def stats(self, /, loc: onp.ToFloat = 0, *, scale: onp.ToFloatND, moment: _Moment1) -> _FloatND: ...
    @overload  # loc: 0-d, scale: n-d (keyword), moments: 2 (default)
    def stats(self, /, loc: onp.ToFloat = 0, *, scale: onp.ToFloatND, moment: _Moment2 = "mv") -> _Tuple2[_FloatND]: ...
    @overload  # loc: 0-d, scale: n-d (keyword), moments: 3
    def stats(self, /, loc: onp.ToFloat = 0, *, scale: onp.ToFloatND, moment: _Moment3) -> _Tuple3[_FloatND]: ...
    @overload  # loc: 0-d, scale: n-d (keyword), moments: 4
    def stats(self, /, loc: onp.ToFloat = 0, *, scale: onp.ToFloatND, moment: _Moment4) -> _Tuple4[_FloatND]: ...

    #
    @overload  # loc: n-d, scale: ?-d, moments: 1 (positional)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND, moment: _Moment1) -> _FloatND: ...
    @overload  # loc: n-d, scale: ?-d, moments: 1 (keyword)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND = 1, *, moment: _Moment1) -> _FloatND: ...
    @overload  # loc: n-d, scale: ?-d, moments: 2 (default)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND = 1, moment: _Moment2 = "mv") -> _Tuple2[_FloatND]: ...
    @overload  # loc: n-d, scale: ?-d, moments: 3 (positional)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND, moment: _Moment3) -> _Tuple3[_FloatND]: ...
    @overload  # loc: n-d, scale: ?-d, moments: 3 (keyword)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND = 1, *, moment: _Moment3) -> _Tuple3[_FloatND]: ...
    @overload  # loc: n-d, scale: ?-d, moments: 4 (positional)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND, moment: _Moment4) -> _Tuple4[_FloatND]: ...
    @overload  # loc: n-d, scale: ?-d, moments: 4 (keyword)
    def stats(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND = 1, *, moment: _Moment4) -> _Tuple4[_FloatND]: ...

    #
    @override
    @overload
    def entropy(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def entropy(self, /, loc: onp.ToFloat, scale: onp.ToFloatND) -> _FloatND: ...
    @overload
    def entropy(self, /, loc: onp.ToFloat = 0, *, scale: onp.ToFloatND) -> _FloatND: ...
    @overload
    def entropy(self, /, loc: onp.ToFloatND, scale: _ToFloatOrND = 1) -> _FloatND: ...

    #
    @override
    @overload
    def moment(self, /, order: onp.ToInt, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def moment(self, /, order: onp.ToInt, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOrND: ...

    #
    @override
    @overload
    def median(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def median(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOrND: ...

    #
    @override
    @overload
    def mean(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def mean(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOrND: ...
    #
    @override
    @overload
    def var(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def var(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOrND: ...

    #
    @override
    @overload
    def std(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def std(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOrND: ...

    #
    @override
    @overload
    def interval(self, /, confidence: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Tuple2[_Float]: ...
    @overload
    def interval(
        self, /, confidence: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1
    ) -> _Tuple2[_Float] | _Tuple2[_FloatND]: ...

    #
    @override
    @overload
    def support(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Tuple2[_Float]: ...
    @overload
    def support(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _Tuple2[_Float] | _Tuple2[_FloatND]: ...

    # overrides of rv_continuous
    @override
    @overload
    def __call__(self, /) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def __call__(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> rv_continuous_frozen[Self]: ...

    #
    @override
    @overload
    def freeze(self, /) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> rv_continuous_frozen[Self, _Float]: ...
    @overload
    def freeze(self, /, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> rv_continuous_frozen[Self]: ...

    #
    @overload
    def pdf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def pdf(self, /, x: onp.ToFloat, loc: _ToFloatOrND, scale: onp.ToFloatND) -> _FloatND: ...
    @overload
    def pdf(self, /, x: onp.ToFloat, loc: _ToFloatOrND = 0, *, scale: onp.ToFloatND) -> _FloatND: ...
    @overload
    def pdf(self, /, x: onp.ToFloat, loc: onp.ToFloatND, scale: _ToFloatOrND) -> _FloatND: ...
    @overload
    def pdf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def pdf(self, /, x: onp.ToFloatND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _Float1ND: ...
    @overload
    def pdf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def logpdf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def logpdf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logpdf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def cdf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def cdf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def cdf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def logcdf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def logcdf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logcdf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def sf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def sf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def sf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def logsf(self, /, x: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def logsf(
        self, /, x: onp.CanArrayND[_CoFloat, _ShapeT], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def logsf(self, /, x: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def ppf(self, /, q: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def ppf(
        self, /, q: onp.Array[_ShapeT, _CoFloat], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def ppf(self, /, q: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    @overload
    def isf(self, /, q: onp.ToFloat, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1) -> _Float: ...
    @overload
    def isf(
        self, /, q: onp.Array[_ShapeT, _CoFloat], loc: onp.ToFloat = 0, scale: onp.ToFloat = 1
    ) -> onp.Array[_ShapeT, np.float64]: ...
    @overload
    def isf(self, /, q: _ToFloatOrND, loc: _ToFloatOrND = 0, scale: _ToFloatOrND = 1) -> _FloatOr1ND: ...

    #
    @override
    def rvs(
        self, /, loc: onp.ToFloat = 0, scale: onp.ToFloat = 1, size: AnyShape = 1, random_state: onp.random.ToRNG | None = None
    ) -> _FloatOrND: ...

    #
    @overload
    def _unpack_loc_scale(self, /, theta: Sequence[onp.ToFloat]) -> tuple[onp.ToFloat, onp.ToFloat, tuple[()]]: ...
    @overload
    def _unpack_loc_scale(self, /, theta: Sequence[onp.ToFloatND]) -> tuple[onp.ToFloatND, onp.ToFloatND, tuple[()]]: ...

# undocumented
def argsreduce(cond: _BoolND, *args: _ToFloatOrND) -> list[_CoFloatND]: ...

# undocumented
def get_distribution_names(namespace_pairs: Iterable[tuple[str, type]], rv_base_class: type) -> _Tuple2[list[str]]: ...
