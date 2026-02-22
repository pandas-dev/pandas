import types
from typing import Any, ClassVar, Final, Generic, Literal, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._covariance import _PSD, Covariance

__all__ = [
    "dirichlet",
    "dirichlet_multinomial",
    "invwishart",
    "matrix_normal",
    "matrix_t",
    "multinomial",
    "multivariate_hypergeom",
    "multivariate_normal",
    "multivariate_t",
    "normal_inverse_gamma",
    "ortho_group",
    "random_correlation",
    "random_table",
    "special_ortho_group",
    "uniform_direction",
    "unitary_group",
    "vonmises_fisher",
    "wishart",
]

_ScalarT = TypeVar("_ScalarT", bound=np.generic, default=np.float64)
_ScalarT_co = TypeVar("_ScalarT_co", bound=np.generic, default=np.float64, covariant=True)

# TODO(@jorenham): rename as {}T_co
_RVG_co = TypeVar("_RVG_co", bound=multi_rv_generic, default=multi_rv_generic, covariant=True)
_RVF_co = TypeVar("_RVF_co", bound=multi_rv_frozen, covariant=True)

_Scalar_f: TypeAlias = npc.floating
_Scalar_uif: TypeAlias = npc.integer | _Scalar_f
_ToFloatMax2D: TypeAlias = onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D
_ToJustFloat: TypeAlias = float | _Scalar_f

_Array1ND: TypeAlias = onp.Array[tuple[int, *tuple[Any, ...]], _ScalarT]
# https://github.com/microsoft/pyright/issues/11127
_Array2ND: TypeAlias = onp.Array[tuple[int, int, *tuple[Any, ...]], _ScalarT]  # pyright: ignore[reportInvalidTypeForm]
_Array3ND: TypeAlias = onp.Array[tuple[int, int, int, *tuple[Any, ...]], _ScalarT]  # pyright: ignore[reportInvalidTypeForm]

_ScalarOrArray_f8: TypeAlias = np.float64 | _Array1ND
_AnyCov: TypeAlias = Covariance | onp.ToFloat2D | onp.ToFloat

@type_check_only
class rng_mixin:
    @property
    def random_state(self, /) -> onp.random.RNG: ...
    @random_state.setter
    def random_state(self, /, seed: onp.random.ToRNG | None) -> None: ...

###

class multi_rv_generic(rng_mixin):
    def __init__(self, /, seed: onp.random.ToRNG | None = None) -> None: ...
    def _get_random_state(self, /, random_state: onp.random.ToRNG | None) -> onp.random.RNG: ...

class multi_rv_frozen(rng_mixin, Generic[_RVG_co]):
    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    _dist: _RVG_co

class multivariate_normal_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
    ) -> multivariate_normal_frozen: ...
    def logpdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat1D | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat1D | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> _ScalarOrArray_f8: ...
    def logcdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> _ScalarOrArray_f8: ...
    def cdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> _ScalarOrArray_f8: ...

    # returns a scalar or array depending on the *value* of size and *mean*
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat | None = None,
        cov: onp.ToFloat = 1,
        size: Literal[1] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat1D | None,
        cov: Covariance | onp.ToFloat2D,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat1D | None = None,
        *,
        cov: Covariance | onp.ToFloat2D,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def rvs(
        self, /, mean: onp.ToFloat1D | None, cov: _AnyCov, size: tuple[int, ...], random_state: onp.random.ToRNG | None = None
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        *,
        size: tuple[int, ...],
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _ScalarOrArray_f8: ...

    #
    def entropy(self, /, mean: onp.ToFloat1D | None = None, cov: _AnyCov = 1) -> np.float64: ...
    def fit(
        self, /, x: onp.ToFloatND, fix_mean: onp.ToFloat1D | None = None, fix_cov: onp.ToFloat2D | None = None
    ) -> tuple[onp.Array1D[np.float64], onp.Array2D[np.float64]]: ...
    def marginal(
        self,
        dimensions: int | onp.ToInt1D,
        mean: onp.ToFloat1D | None = None,
        cov: onp.ToFloat | onp.ToFloat2D = 1,
        allow_singular: bool = False,
    ) -> multivariate_normal_frozen: ...

# TODO(@jorenham): Generic shape-type for mean and cov, so that we can determine whether the methods return scalars or arrays.
# https://github.com/scipy/scipy-stubs/issues/406
class multivariate_normal_frozen(multi_rv_frozen[multivariate_normal_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    dim: Final[int]
    allow_singular: Final[bool]
    maxpts: Final[int]
    abseps: Final[float]
    releps: Final[float]
    cov_object: Final[Covariance]
    mean: onp.Array1D[np.float64]

    @property
    def cov(self, /) -> onp.Array2D[np.float64]: ...

    #
    def __init__(
        self,
        /,
        mean: onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
        maxpts: onp.ToJustInt | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
    ) -> None: ...

    #
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def logcdf(
        self, /, x: onp.ToFloatND, *, lower_limit: onp.ToFloat1D | None = None, rng: onp.random.ToRNG | None = None
    ) -> _ScalarOrArray_f8: ...
    def cdf(
        self, /, x: onp.ToFloatND, *, lower_limit: onp.ToFloat1D | None = None, rng: onp.random.ToRNG | None = None
    ) -> _ScalarOrArray_f8: ...

    # If mean and cov are 0-D, and size = 1, this returns a scalar. But without knowing the shape of mean and cov at type-check
    # time, we cannot determine that. See https://github.com/scipy/scipy-stubs/issues/406 for more info.
    @overload
    def rvs(self, /, size: int = 1, random_state: onp.random.ToRNG | None = None) -> _ScalarOrArray_f8: ...
    @overload
    def rvs(self, /, size: tuple[int, ...], random_state: onp.random.ToRNG | None = None) -> onp.ArrayND[np.float64]: ...

    #
    def entropy(self, /) -> np.float64: ...
    def marginal(self, dimensions: int | onp.ToInt1D) -> multivariate_normal_frozen: ...

class matrix_normal_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
        seed: onp.random.ToRNG | None = None,
    ) -> matrix_normal_frozen: ...

    #
    def logpdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
    ) -> _ScalarOrArray_f8: ...

    # If `size > 1` the output is 3-D, otherwise 2-D.
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
        size: Literal[1] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None,
        rowcov: onp.ToFloat2D | onp.ToFloat,
        colcov: onp.ToFloat2D | onp.ToFloat,
        size: int,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
        *,
        size: int,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...

    #
    def entropy(self, /, rowcov: _AnyCov = 1, colcov: _AnyCov = 1) -> np.float64: ...

class matrix_normal_frozen(multi_rv_frozen[matrix_normal_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    rowpsd: Final[_PSD]
    colpsd: Final[_PSD]

    def __init__(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat2D | onp.ToFloat = 1,
        colcov: onp.ToFloat2D | onp.ToFloat = 1,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...
    def logpdf(self, /, X: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, X: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    @overload
    def rvs(self, /, size: Literal[1] = 1, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, size: int, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...

    #
    def entropy(self, /) -> np.float64: ...

class matrix_t_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        row_spread: onp.ToFloat2D | float = 1,
        col_spread: onp.ToFloat2D | float = 1,
        df: float | None = None,
        seed: onp.random.ToRNG | None = None,
    ) -> matrix_t_frozen: ...

    #
    def logpdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        row_spread: onp.ToFloat2D | float = 1,
        col_spread: onp.ToFloat2D | float = 1,
        df: float = 1,
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        row_spread: onp.ToFloat2D | float = 1,
        col_spread: onp.ToFloat2D | float = 1,
        df: float = 1,
    ) -> _ScalarOrArray_f8: ...

    # If `size > 1` the output is 3-D, otherwise 2-D.
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        row_spread: onp.ToFloat2D | float = 1,
        col_spread: onp.ToFloat2D | float = 1,
        df: float = 1,
        size: Literal[1] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        row_spread: onp.ToFloat2D | float = 1,
        col_spread: onp.ToFloat2D | float = 1,
        df: float = 1,
        *,
        size: int,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...

class matrix_t_frozen(multi_rv_frozen[matrix_t_gen]):
    rowpsd: Final[_PSD]
    colpsd: Final[_PSD]

    mean: Final[onp.Array2D[np.float64]]
    df: Final[float]

    def __init__(
        self,
        /,
        mean: onp.ToFloat2D,
        row_spread: onp.ToFloat2D | float,
        col_spread: onp.ToFloat2D | float,
        df: float,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...

    #
    def logpdf(self, /, X: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, X: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    @overload
    def rvs(self, /, size: Literal[1] = 1, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, size: int, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...

class dirichlet_gen(multi_rv_generic):
    def __call__(self, /, alpha: onp.ToFloat1D, seed: onp.random.ToRNG | None = None) -> dirichlet_frozen: ...
    def logpdf(self, /, x: onp.ToFloatND, alpha: onp.ToFloat1D) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND, alpha: onp.ToFloat1D) -> _ScalarOrArray_f8: ...
    def mean(self, /, alpha: onp.ToFloat1D) -> onp.Array1D[np.float64]: ...
    def var(self, /, alpha: onp.ToFloat1D) -> onp.Array1D[np.float64]: ...
    def cov(self, /, alpha: onp.ToFloat1D) -> onp.Array2D[np.float64]: ...
    def entropy(self, /, alpha: onp.ToFloat1D) -> np.float64: ...
    @overload
    def rvs(
        self, /, alpha: onp.ToFloat1D, size: tuple[()], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self, /, alpha: onp.ToFloat1D, size: op.CanIndex | tuple[op.CanIndex] = 1, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self, /, alpha: onp.ToFloat1D, size: tuple[op.CanIndex, op.CanIndex], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        alpha: onp.ToFloat1D,
        size: tuple[op.CanIndex, op.CanIndex, op.CanIndex, *tuple[op.CanIndex, ...]],
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND: ...

class dirichlet_frozen(multi_rv_frozen[dirichlet_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    alpha: Final[onp.Array1D[_Scalar_uif]]

    def __init__(self, /, alpha: onp.ToFloat1D, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> onp.Array1D[np.float64]: ...
    def var(self, /) -> onp.Array1D[np.float64]: ...
    def cov(self, /) -> onp.Array2D[np.float64]: ...
    def entropy(self, /) -> np.float64: ...
    @overload
    def rvs(self, /, size: tuple[()], random_state: onp.random.ToRNG | None = None) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self, /, size: op.CanIndex | tuple[op.CanIndex] = 1, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self, /, size: tuple[op.CanIndex, op.CanIndex], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        size: tuple[op.CanIndex, op.CanIndex, op.CanIndex, *tuple[op.CanIndex, ...]],
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND: ...

class wishart_gen(multi_rv_generic):
    def __call__(
        self, /, df: onp.ToFloat | None = None, scale: _ToFloatMax2D | None = None, seed: onp.random.ToRNG | None = None
    ) -> wishart_frozen: ...

    #
    def logpdf(self, /, x: onp.ToFloatND, df: onp.ToFloat, scale: _ToFloatMax2D) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND, df: onp.ToFloat, scale: _ToFloatMax2D) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | onp.Array2D[np.float64]: ...
    def mode(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | None: ...
    def var(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | onp.Array2D[np.float64]: ...
    def entropy(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64: ...

    #
    @overload
    def rvs(
        self, /, df: onp.ToFloat, scale: onp.ToFloat, size: Literal[1] = 1, random_state: onp.random.ToRNG | None = None
    ) -> np.float64: ...
    @overload
    def rvs(
        self, /, df: onp.ToFloat, scale: _ToFloatMax2D, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None
    ) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        df: onp.ToFloat,
        scale: _ToFloatMax2D,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64] | np.float64: ...

class wishart_frozen(multi_rv_frozen[wishart_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    dim: Final[int]
    df: Final[onp.ToFloat]
    scale: Final[onp.Array2D[np.float64]]
    C: Final[onp.Array2D[np.float64]]
    log_det_scale: Final[float]

    def __init__(self, /, df: onp.ToFloat, scale: _ToFloatMax2D, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
    def mode(self, /) -> np.float64 | None: ...
    def var(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
    def entropy(self, /) -> np.float64: ...

    #
    @overload
    def rvs(self, /, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(self, /, size: int | tuple[int, ...] = 1, random_state: onp.random.ToRNG | None = None) -> _ScalarOrArray_f8: ...

class invwishart_gen(wishart_gen):
    @override
    def __call__(  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
        self, /, df: onp.ToFloat | None = None, scale: _ToFloatMax2D | None = None, seed: onp.random.ToRNG | None = None
    ) -> invwishart_frozen: ...  # ty: ignore[invalid-method-override]
    @override
    def mean(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | None: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @override
    def mode(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | onp.Array2D[np.float64]: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @override
    def var(self, /, df: onp.ToFloat, scale: _ToFloatMax2D) -> np.float64 | None: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]

class invwishart_frozen(multi_rv_frozen[invwishart_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    def __init__(self, /, df: onp.ToFloat, scale: _ToFloatMax2D, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
    def mode(self, /) -> np.float64 | None: ...
    def var(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
    def entropy(self, /) -> np.float64: ...

    #
    @overload
    def rvs(self, /, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(self, /, size: int | tuple[int, ...] = 1, random_state: onp.random.ToRNG | None = None) -> _ScalarOrArray_f8: ...

# NOTE: `n` and `p` are broadcast-able (although this breaks `.rvs()` at runtime...)
class multinomial_gen(multi_rv_generic):
    def __call__(
        self, /, n: onp.ToJustIntND, p: onp.ToJustFloatND, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen: ...

    #
    def logpmf(self, /, x: onp.ToFloatND, n: onp.ToJustIntND, p: onp.ToJustFloatND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND, n: onp.ToJustIntND, p: onp.ToJustFloatND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /, n: onp.ToJustIntND, p: onp.ToJustFloatND) -> _Array1ND: ...
    def cov(self, /, n: onp.ToJustIntND, p: onp.ToJustFloatND) -> _Array2ND: ...
    def entropy(self, /, n: onp.ToJustIntND, p: onp.ToJustFloatND) -> _ScalarOrArray_f8: ...

    #
    @overload
    def rvs(
        self, /, n: int | onp.ToJustIntND, p: onp.ToJustFloatND, size: tuple[()], random_state: onp.random.ToRNG | None = None
    ) -> _Array1ND: ...
    @overload
    def rvs(
        self,
        /,
        n: onp.ToJustIntND,
        p: onp.ToJustFloatND,
        size: int | onp.AtLeast1D | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND: ...

#
class multinomial_frozen(multi_rv_frozen[multinomial_gen]):
    def __init__(self, /, n: onp.ToJustIntND, p: onp.ToJustFloatND, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    def logpmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /) -> _Array1ND: ...
    def cov(self, /) -> _Array2ND: ...
    def entropy(self, /) -> _ScalarOrArray_f8: ...

    #
    @overload
    def rvs(self, /, size: tuple[()], random_state: onp.random.ToRNG | None = None) -> _Array1ND: ...
    @overload
    def rvs(self, /, size: onp.AtLeast1D | int = 1, random_state: onp.random.ToRNG | None = None) -> _Array2ND: ...

@type_check_only
class _group_rv_gen_mixin(Generic[_RVF_co, _ScalarT_co]):
    # NOTE: Contrary to what the `dim` default suggests, it is required.
    def __call__(self, /, dim: int | None = None, seed: onp.random.ToRNG | None = None) -> _RVF_co: ...

    #
    @overload
    def rvs(
        self, /, dim: int, size: Literal[0, 1] = 1, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[_ScalarT_co]: ...
    @overload
    def rvs(self, /, dim: int, size: int, random_state: onp.random.ToRNG | None = None) -> _Array2ND[_ScalarT_co]: ...

@type_check_only
class _group_rv_frozen_mixin(Generic[_ScalarT_co]):
    __class_getitem__: ClassVar[None] = None

    dim: Final[int]

    # NOTE: Contrary to what the `dim` default suggests, it is required.
    def __init__(self, /, dim: int | None = None, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    @overload
    def rvs(self, /, size: Literal[0, 1] = 1, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[_ScalarT_co]: ...
    @overload
    def rvs(self, /, size: int, random_state: onp.random.ToRNG | None = None) -> _Array2ND[_ScalarT_co]: ...

class special_ortho_group_gen(_group_rv_gen_mixin[special_ortho_group_frozen], multi_rv_generic): ...
class special_ortho_group_frozen(_group_rv_frozen_mixin, multi_rv_frozen[special_ortho_group_gen]): ...  # type: ignore[misc]
class ortho_group_gen(_group_rv_gen_mixin[ortho_group_frozen], multi_rv_generic): ...
class ortho_group_frozen(_group_rv_frozen_mixin, multi_rv_frozen[ortho_group_gen]): ...  # type: ignore[misc]
class unitary_group_gen(_group_rv_gen_mixin[unitary_group_frozen, np.complex128], multi_rv_generic): ...
class unitary_group_frozen(_group_rv_frozen_mixin[np.complex128], multi_rv_frozen[unitary_group_gen]): ...  # type: ignore[misc]

class uniform_direction_gen(multi_rv_generic):
    def __call__(self, /, dim: int | None = None, seed: onp.random.ToRNG | None = None) -> uniform_direction_frozen: ...

    #
    @overload
    def rvs(self, /, dim: int, size: None = None, random_state: onp.random.ToRNG | None = None) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self, /, dim: int, size: int | tuple[int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, dim: int, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array3ND[np.float64]: ...

class uniform_direction_frozen(multi_rv_frozen[uniform_direction_gen]):
    dim: Final[int]

    def __init__(self, /, dim: int | None = None, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    @overload
    def rvs(self, /, size: None = None, random_state: onp.random.ToRNG | None = None) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(self, /, size: int | tuple[int], random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...

class random_correlation_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        eigs: onp.ToFloat1D,
        seed: onp.random.ToRNG | None = None,
        tol: _ToJustFloat = 1e-13,
        diag_tol: _ToJustFloat = 1e-7,
    ) -> random_correlation_frozen: ...

    #
    def rvs(
        self,
        /,
        eigs: onp.ToFloat1D,
        random_state: onp.random.ToRNG | None = None,
        tol: _ToJustFloat = 1e-13,
        diag_tol: _ToJustFloat = 1e-7,
    ) -> onp.Array2D[np.float64]: ...

class random_correlation_frozen(multi_rv_frozen[random_correlation_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    eigs: Final[onp.Array1D[np.float64]]
    tol: Final[float]
    diag_tol: Final[float]

    def __init__(
        self,
        /,
        eigs: onp.ToFloat1D,
        seed: onp.random.ToRNG | None = None,
        tol: _ToJustFloat = 1e-13,
        diag_tol: _ToJustFloat = 1e-7,
    ) -> None: ...

    #
    def rvs(self, /, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.float64]: ...

class multivariate_t_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: onp.ToJustInt = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
    ) -> multivariate_t_frozen: ...

    #
    def logpdf(
        self, /, x: onp.ToFloatND, loc: onp.ToFloat1D | None = None, shape: onp.ToFloat | onp.ToFloat2D = 1, df: int = 1
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self,
        /,
        x: onp.ToFloatND,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        allow_singular: bool = False,
    ) -> _ScalarOrArray_f8: ...
    def cdf(
        self,
        /,
        x: onp.ToFloatND,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        allow_singular: bool = False,
        *,
        maxpts: int | None = None,
        lower_limit: onp.ToFloat1D | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _ScalarOrArray_f8: ...

    #
    def entropy(self, /, loc: onp.ToFloat1D | None = None, shape: onp.ToFloat | onp.ToFloat2D = 1, df: int = 1) -> np.float64: ...

    #
    @overload
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        *,
        size: tuple[()],
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None,
        shape: onp.ToFloat | onp.ToFloat2D,
        df: int,
        size: tuple[()],
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        size: int | tuple[int] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        *,
        size: onp.AtLeast2D,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND: ...
    @overload
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None,
        shape: onp.ToFloat | onp.ToFloat2D,
        df: int,
        size: onp.AtLeast2D,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND: ...

    #
    def marginal(
        self,
        dimensions: int | onp.ToInt1D,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        allow_singular: bool = False,
    ) -> multivariate_t_frozen: ...

class multivariate_t_frozen(multi_rv_frozen[multivariate_t_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    dim: Final[int]
    df: Final[int]
    loc: Final[onp.Array1D[np.float64]]
    shape: Final[onp.Array2D[np.float64]]
    shape_info: Final[_PSD]

    def __init__(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...

    #
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def cdf(
        self,
        /,
        x: onp.ToFloatND,
        *,
        maxpts: int | None = None,
        lower_limit: onp.ToFloat1D | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _ScalarOrArray_f8: ...

    #
    def entropy(self, /) -> np.float64: ...

    #
    @overload
    def rvs(self, /, size: tuple[()], random_state: onp.random.ToRNG | None = None) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(self, /, size: int | tuple[int] = 1, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array3ND: ...

    #
    def marginal(self, dimensions: int | onp.ToInt1D) -> multivariate_t_frozen: ...

# NOTE: `m` and `n` are broadcastable (but doing so will break `.rvs()` at runtime...)
class multivariate_hypergeom_gen(multi_rv_generic):
    def __call__(
        self, /, m: onp.ToJustIntND, n: onp.ToJustInt | onp.ToJustIntND, seed: onp.random.ToRNG | None = None
    ) -> multivariate_hypergeom_frozen: ...
    def logpmf(self, /, x: onp.ToFloatND, m: onp.ToJustIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array1ND: ...
    def var(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array1ND: ...
    def cov(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array2ND: ...
    @overload
    def rvs(
        self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND, size: tuple[()], random_state: onp.random.ToRNG | None = None
    ) -> _Array1ND: ...
    @overload
    def rvs(
        self,
        /,
        m: onp.ToIntND,
        n: onp.ToJustInt | onp.ToJustIntND,
        size: int | onp.AtLeast1D | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND: ...

class multivariate_hypergeom_frozen(multi_rv_frozen[multivariate_hypergeom_gen]):
    def __init__(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> _Array1ND: ...
    def var(self, /) -> _Array1ND: ...
    def cov(self, /) -> _Array2ND: ...
    @overload
    def rvs(self, /, size: tuple[()], random_state: onp.random.ToRNG | None = None) -> _Array1ND: ...
    @overload
    def rvs(self, /, size: int | tuple[int] = 1, random_state: onp.random.ToRNG | None = None) -> _Array2ND: ...

_RandomTableRVSMethod: TypeAlias = Literal["boyett", "patefield"]

class random_table_gen(multi_rv_generic):
    def __call__(
        self, /, row: onp.ToJustIntND, col: onp.ToJustIntND, *, seed: onp.random.ToRNG | None = None
    ) -> random_table_frozen: ...

    #
    def logpmf(self, /, x: onp.ToFloatND, row: onp.ToJustIntND, col: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND, row: onp.ToJustIntND, col: onp.ToJustIntND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /, row: onp.ToJustIntND, col: onp.ToJustIntND) -> onp.Array2D[np.float64]: ...

    #
    @overload
    def rvs(
        self,
        /,
        row: onp.ToJustIntND,
        col: onp.ToJustIntND,
        *,
        size: None = None,
        method: _RandomTableRVSMethod | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        row: onp.ToJustIntND,
        col: onp.ToJustIntND,
        *,
        size: int | tuple[int],
        method: _RandomTableRVSMethod | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        row: onp.ToJustIntND,
        col: onp.ToJustIntND,
        *,
        size: onp.AtLeast1D,
        method: _RandomTableRVSMethod | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND[np.float64]: ...

class random_table_frozen(multi_rv_frozen[random_table_gen]):
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    def __init__(self, /, row: onp.ToJustIntND, col: onp.ToJustIntND, *, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    def logpmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /) -> onp.Array2D[np.float64]: ...

    #
    @overload
    def rvs(
        self, /, size: None = None, method: _RandomTableRVSMethod | None = None, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self, /, size: int | tuple[int], method: _RandomTableRVSMethod | None = None, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self, /, size: onp.AtLeast1D, method: _RandomTableRVSMethod | None = None, random_state: onp.random.ToRNG | None = None
    ) -> _Array3ND[np.float64]: ...

class dirichlet_multinomial_gen(multi_rv_generic):
    def __call__(
        self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND, seed: onp.random.ToRNG | None = None
    ) -> dirichlet_multinomial_frozen: ...
    def logpmf(self, /, x: onp.ToIntND, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToIntND, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array1ND: ...
    def var(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array1ND: ...
    def cov(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array2ND: ...

class dirichlet_multinomial_frozen(multi_rv_frozen[dirichlet_multinomial_gen]):
    alpha: _Array1ND
    n: _Array1ND[np.int_]  # broadcasted against alpha

    def __init__(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpmf(self, /, x: onp.ToIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> _Array1ND: ...
    def var(self, /) -> _Array1ND: ...
    def cov(self, /) -> _Array2ND: ...

class vonmises_fisher_gen(multi_rv_generic):
    def __call__(
        self, /, mu: onp.ToFloat1D | None = None, kappa: int = 1, seed: onp.random.ToRNG | None = None
    ) -> vonmises_fisher_frozen: ...
    def logpdf(self, /, x: onp.ToFloatND, mu: onp.ToFloat1D | None = None, kappa: int = 1) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND, mu: onp.ToFloat1D | None = None, kappa: int = 1) -> _ScalarOrArray_f8: ...
    def entropy(self, /, mu: onp.ToFloat1D | None = None, kappa: int = 1) -> np.float64: ...
    def rvs(
        self,
        /,
        mu: onp.ToFloat1D | None = None,
        kappa: int = 1,
        size: op.CanIndex | tuple[op.CanIndex, *tuple[op.CanIndex, ...]] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND: ...
    def fit(self, /, x: onp.ToFloatND) -> tuple[onp.Array1D[np.float64], float]: ...

class vonmises_fisher_frozen(multi_rv_frozen[vonmises_fisher_gen]):
    def __init__(self, /, mu: onp.ToFloat1D | None = None, kappa: int = 1, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def entropy(self, /) -> np.float64: ...
    def rvs(
        self,
        /,
        size: op.CanIndex | tuple[op.CanIndex, *tuple[op.CanIndex, ...]] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND: ...

class normal_inverse_gamma_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
        seed: onp.random.ToRNG | None = None,
    ) -> normal_inverse_gamma_frozen: ...
    def logpdf(
        self,
        /,
        x: onp.ToFloat | onp.ToFloatND,
        s2: onp.ToFloat | onp.ToFloatND,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self,
        /,
        x: onp.ToFloat | onp.ToFloatND,
        s2: onp.ToFloat | onp.ToFloatND,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
    ) -> _ScalarOrArray_f8: ...
    def mean(
        self,
        /,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
    ) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...
    def var(
        self,
        /,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
    ) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...
    def rvs(
        self,
        /,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
        size: op.CanIndex | tuple[op.CanIndex, ...] | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...

#
class normal_inverse_gamma_frozen(multi_rv_frozen[normal_inverse_gamma_gen]):
    def __init__(
        self,
        /,
        mu: onp.ToFloat | onp.ToFloatND = 0,
        lmbda: onp.ToFloat | onp.ToFloatND = 1,
        a: onp.ToFloat | onp.ToFloatND = 1,
        b: onp.ToFloat | onp.ToFloatND = 1,
        seed: onp.random.ToRNG | None = None,
    ) -> None: ...
    def logpdf(self, /, x: onp.ToFloat | onp.ToFloatND, s2: onp.ToFloat | onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloat | onp.ToFloatND, s2: onp.ToFloat | onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...
    def var(self, /) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...
    def rvs(
        self, /, size: op.CanIndex | tuple[op.CanIndex, ...] | None = None, random_state: onp.random.ToRNG | None = None
    ) -> tuple[_ScalarOrArray_f8, _ScalarOrArray_f8]: ...

multivariate_normal: Final[multivariate_normal_gen] = ...
matrix_normal: Final[matrix_normal_gen] = ...
matrix_t: Final[matrix_t_gen] = ...
dirichlet: Final[dirichlet_gen] = ...
wishart: Final[wishart_gen] = ...
invwishart: Final[invwishart_gen] = ...
multinomial: Final[multinomial_gen] = ...
special_ortho_group: Final[special_ortho_group_gen] = ...
ortho_group: Final[ortho_group_gen] = ...
random_correlation: Final[random_correlation_gen] = ...
unitary_group: Final[unitary_group_gen] = ...
multivariate_t: Final[multivariate_t_gen] = ...
multivariate_hypergeom: Final[multivariate_hypergeom_gen] = ...
random_table: Final[random_table_gen] = ...
uniform_direction: Final[uniform_direction_gen] = ...
dirichlet_multinomial: Final[dirichlet_multinomial_gen] = ...
vonmises_fisher: Final[vonmises_fisher_gen] = ...
normal_inverse_gamma: Final[normal_inverse_gamma_gen] = ...
