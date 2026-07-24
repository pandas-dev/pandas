# there are `overload-overlap` false positivies on `numpy<2.1` that can't selectively be ignored or worked around
# mypy: disable-error-code=overload-overlap

import types
from typing import Any, ClassVar, Final, Generic, Literal, Never, SupportsIndex, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
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

###

type _ToFloatMax2D = onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D
type _ToJustFloat = float | npc.floating

type _Array1ND[ScalarT: np.generic] = onp.Array[tuple[int, *tuple[Any, ...]], ScalarT]
type _Array2ND[ScalarT: np.generic] = onp.Array[tuple[int, int, *tuple[Any, ...]], ScalarT]
type _Array3ND[ScalarT: np.generic] = onp.Array[tuple[int, int, int, *tuple[Any, ...]], ScalarT]

type _ScalarOrArray_f8 = np.float64 | _Array1ND[np.float64]
type _AnyCov = Covariance | onp.ToFloat2D | onp.ToFloat

type _ToIntStrict1D = onp.ToArrayStrict1D[int, npc.integer]
type _ToIntStrict2D = onp.ToArrayStrict2D[int, npc.integer]
type _ToIntND = onp.ToArrayND[int, npc.integer]

type _ToFloatStrict1D = onp.ToArrayStrict1D[float, npc.floating]
type _ToFloatStrict2D = onp.ToArrayStrict2D[float, npc.floating]
type _ToFloatStrict3D = onp.ToArrayStrict3D[float, npc.floating]
type _ToFloatND = onp.ToArrayND[float, npc.floating]

# workaround for https://github.com/microsoft/pyright/issues/10232
type _JustAnyShape = tuple[Never, Never, Never, Never]
# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
type _WorkaroundForPyright = tuple[int] | tuple[Any, ...]

_ScalarT_co = TypeVar("_ScalarT_co", bound=np.generic, default=np.float64, covariant=True)

# TODO(@jorenham): rename as {}T_co
_RVG_co = TypeVar("_RVG_co", bound=multi_rv_generic, default=multi_rv_generic, covariant=True)
_RVF_co = TypeVar("_RVF_co", bound=multi_rv_frozen, covariant=True)

_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

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
    @overload
    def __call__(
        self,
        /,
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
    ) -> multivariate_normal_frozen[tuple[()]]: ...
    @overload
    def __call__(
        self, /, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False, seed: onp.random.ToRNG | None = None
    ) -> multivariate_normal_frozen[tuple[int]]: ...

    #
    @overload
    def logpdf(
        self, /, x: onp.ToFloat, mean: onp.ToFloat | onp.ToFloat1D | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64: ...
    @overload
    def logpdf[ShapeT: tuple[int, ...]](
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def logpdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def logpdf(
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def logpdf(
        self, /, x: onp.ToFloatStrict1D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64: ...
    @overload
    def logpdf(
        self, /, x: onp.ToFloatStrict2D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def logpdf(
        self, /, x: onp.ToFloatStrict3D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def logpdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def pdf(
        self, /, x: onp.ToFloat, mean: onp.ToFloat | onp.ToFloat1D | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64: ...
    @overload
    def pdf[ShapeT: tuple[int, ...]](
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def pdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def pdf(
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def pdf(
        self, /, x: onp.ToFloatStrict1D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64: ...
    @overload
    def pdf(
        self, /, x: onp.ToFloatStrict2D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def pdf(
        self, /, x: onp.ToFloatStrict3D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def pdf(
        self, /, x: onp.ToFloatND, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloat,
        mean: onp.ToFloat | onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def logcdf[ShapeT: tuple[int, ...]](
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloatStrict1D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloatStrict2D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloatStrict3D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | onp.ToFloat3D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def logcdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloat,
        mean: onp.ToFloat | onp.ToFloat1D | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def cdf[ShapeT: tuple[int, ...]](
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloatStrict1D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloatStrict2D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloatStrict3D,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | onp.ToFloat3D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def cdf(
        self,
        /,
        x: onp.ToFloatND,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
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
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    def entropy(self, /, mean: onp.ToFloat | onp.ToFloat1D | None = None, cov: _AnyCov = 1) -> np.float64: ...

    #
    def fit(
        self,
        /,
        x: onp.ToFloatND,
        fix_mean: onp.ToFloat | onp.ToFloat1D | None = None,
        fix_cov: onp.ToFloat1D | onp.ToFloat2D | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.Array2D[np.float64]]: ...

    #
    @overload
    def marginal(
        self, dimensions: int, mean: onp.ToFloat | onp.ToFloat1D | None = None, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> multivariate_normal_frozen[tuple[()]]: ...
    @overload
    def marginal(
        self, dimensions: onp.ToInt1D, mean: onp.ToFloat1D, cov: _AnyCov = 1, allow_singular: bool = False
    ) -> multivariate_normal_frozen[tuple[int]]: ...

#
class multivariate_normal_frozen(multi_rv_frozen[multivariate_normal_gen], Generic[_ShapeT_co]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    dim: Final[int]
    allow_singular: Final[bool]
    maxpts: Final[int]
    abseps: Final[float]
    releps: Final[float]
    cov_object: Final[Covariance[np.float64]]
    mean: Final[onp.Array1D[np.float64]]

    @property
    def cov(self, /) -> onp.Array2D[np.float64]: ...

    #
    @overload
    def __init__(
        self: multivariate_normal_frozen[tuple[()]],
        /,
        mean: onp.ToFloat | None = None,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
    ) -> None: ...
    @overload
    def __init__(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        mean: onp.ToFloat1D,
        cov: _AnyCov = 1,
        allow_singular: bool = False,
        seed: onp.random.ToRNG | None = None,
        maxpts: int | None = None,
        abseps: float = 1e-5,
        releps: float = 1e-5,
    ) -> None: ...

    #
    @overload
    def logpdf(self, /, x: onp.ToFloat) -> np.float64: ...
    @overload
    def logpdf[ShapeT: tuple[int, ...]](
        self: multivariate_normal_frozen[tuple[()]], /, x: onp.ArrayND[npc.floating | npc.integer, ShapeT]
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def logpdf(
        self: multivariate_normal_frozen[tuple[()]], /, x: onp.ToFloatND
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def logpdf(
        self: multivariate_normal_frozen[tuple[int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def logpdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict1D) -> np.float64: ...
    @overload
    def logpdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload
    def logpdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload
    def logpdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatND) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def pdf(self, /, x: onp.ToFloat) -> np.float64: ...
    @overload
    def pdf[ShapeT: tuple[int, ...]](
        self: multivariate_normal_frozen[tuple[()]], /, x: onp.ArrayND[npc.floating | npc.integer, ShapeT]
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def pdf(
        self: multivariate_normal_frozen[tuple[()]], /, x: onp.ToFloatND
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def pdf(
        self: multivariate_normal_frozen[tuple[int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def pdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict1D) -> np.float64: ...
    @overload
    def pdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload
    def pdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload
    def pdf(self: multivariate_normal_frozen[tuple[int]], /, x: onp.ToFloatND) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def logcdf(
        self, /, x: onp.ToFloat, *, lower_limit: onp.ToFloat | None = None, rng: onp.random.ToRNG | None = None
    ) -> np.float64: ...
    @overload
    def logcdf[ShapeT: tuple[int, ...]](
        self: multivariate_normal_frozen[tuple[()]],
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[()]],
        /,
        x: onp.ToFloatND,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict1D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict2D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict3D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | onp.ToFloat3D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def logcdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatND,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def cdf(
        self, /, x: onp.ToFloat, *, lower_limit: onp.ToFloat | None = None, rng: onp.random.ToRNG | None = None
    ) -> np.float64: ...
    @overload
    def cdf[ShapeT: tuple[int, ...]](
        self: multivariate_normal_frozen[tuple[()]],
        /,
        x: onp.ArrayND[npc.floating | npc.integer, ShapeT],
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[()]],
        /,
        x: onp.ToFloatND,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict1D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict2D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatStrict3D,
        *,
        lower_limit: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D | onp.ToFloat3D | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def cdf(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        x: onp.ToFloatND,
        *,
        lower_limit: onp.ToFloat | onp.ToFloatND | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def rvs(
        self: multivariate_normal_frozen[tuple[()]],
        /,
        size: Literal[1] | tuple[()] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload
    def rvs(
        self: multivariate_normal_frozen[tuple[()]], /, size: int | tuple[int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self: multivariate_normal_frozen[tuple[int]],
        /,
        size: Literal[1] | tuple[()] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def rvs(
        self: multivariate_normal_frozen[tuple[int]], /, size: int | tuple[int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(self, /, size: onp.AtLeast0D, random_state: onp.random.ToRNG | None = None) -> onp.ArrayND[np.float64]: ...

    #
    def entropy(self, /) -> np.float64: ...

    #
    @overload
    def marginal(self, dimensions: int) -> multivariate_normal_frozen[tuple[()]]: ...
    @overload
    def marginal(
        self: multivariate_normal_frozen[tuple[int]], dimensions: onp.ToInt1D
    ) -> multivariate_normal_frozen[tuple[int]]: ...

class matrix_normal_gen(multi_rv_generic):
    def __call__(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        seed: onp.random.ToRNG | None = None,
    ) -> matrix_normal_frozen: ...

    #
    def logpdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
    ) -> _ScalarOrArray_f8: ...
    def pdf(
        self,
        /,
        X: onp.ToFloatND,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
    ) -> _ScalarOrArray_f8: ...

    # If `size > 1` the output is 3-D, otherwise 2-D.
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        size: Literal[1] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D,
        size: int,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        *,
        size: int,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...

    #
    def entropy(
        self, /, rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1, colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1
    ) -> np.float64: ...

class matrix_normal_frozen(multi_rv_frozen[matrix_normal_gen]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    rowpsd: Final[_PSD]
    colpsd: Final[_PSD]

    def __init__(
        self,
        /,
        mean: onp.ToFloat2D | None = None,
        rowcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
        colcov: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D = 1,
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
        self,
        /,
        alpha: onp.ToFloat1D,
        size: SupportsIndex | tuple[SupportsIndex] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self, /, alpha: onp.ToFloat1D, size: tuple[SupportsIndex, SupportsIndex], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        alpha: onp.ToFloat1D,
        size: tuple[SupportsIndex, SupportsIndex, SupportsIndex, *tuple[SupportsIndex, ...]],
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND[np.float64]: ...

class dirichlet_frozen(multi_rv_frozen[dirichlet_gen]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    alpha: Final[onp.Array1D[npc.integer | npc.floating]]

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
        self, /, size: SupportsIndex | tuple[SupportsIndex] = 1, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def rvs(
        self, /, size: tuple[SupportsIndex, SupportsIndex], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def rvs(
        self,
        /,
        size: tuple[SupportsIndex, SupportsIndex, SupportsIndex, *tuple[SupportsIndex, ...]],
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array3ND[np.float64]: ...

class wishart_gen(multi_rv_generic):
    def __call__(
        self, /, df: onp.ToFloat | None = None, scale: _ToFloatMax2D | None = None, seed: onp.random.ToRNG | None = None
    ) -> wishart_frozen: ...

    #
    def logpdf(self, /, x: onp.ToFloatND, df: onp.ToFloat, scale: _ToFloatMax2D) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND, df: onp.ToFloat, scale: _ToFloatMax2D) -> _ScalarOrArray_f8: ...

    #
    @overload
    def mode(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64 | None: ...
    @overload
    def mode(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64] | None: ...

    #
    @overload
    def mean(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64: ...
    @overload
    def mean(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...

    #
    @overload
    def var(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64: ...
    @overload
    def var(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...

    #
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
    # pyrefly: ignore [bad-override]
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
    def mode(self, /) -> np.float64 | onp.Array2D[np.float64] | None: ...
    def mean(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
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

    #
    @override
    @overload
    def mode(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64: ...
    @overload
    def mode(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...

    #
    @override  # type:ignore[override]
    @overload
    def mean(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64 | None: ...  # pyrefly:ignore[bad-override]
    @overload
    def mean(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64] | None: ...  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]

    #
    @override  # type:ignore[override]
    @overload
    def var(self, /, df: onp.ToFloat, scale: onp.ToFloat) -> np.float64 | None: ...  # pyrefly:ignore[bad-override]
    @overload
    def var(self, /, df: onp.ToFloat, scale: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64] | None: ...  # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]

class invwishart_frozen(multi_rv_frozen[invwishart_gen]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    def __init__(self, /, df: onp.ToFloat, scale: _ToFloatMax2D, seed: onp.random.ToRNG | None = None) -> None: ...

    #
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...

    #
    def mode(self, /) -> np.float64 | onp.Array2D[np.float64]: ...
    def mean(self, /) -> np.float64 | onp.Array2D[np.float64] | None: ...
    def var(self, /) -> np.float64 | onp.Array2D[np.float64] | None: ...
    def entropy(self, /) -> np.float64: ...

    #
    @overload
    def rvs(self, /, size: onp.AtLeast2D, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.float64]: ...
    @overload
    def rvs(self, /, size: int | tuple[int, ...] = 1, random_state: onp.random.ToRNG | None = None) -> _ScalarOrArray_f8: ...

class multinomial_gen(multi_rv_generic):
    @overload  # n: 0d, p: Nd -> Nd
    def __call__[ShapeT: tuple[int, ...]](
        self, /, n: int, p: onp.ArrayND[npc.floating, ShapeT], seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[ShapeT]: ...
    @overload  # n: 0d, p: 1d -> 1d
    def __call__(
        self, /, n: int, p: _ToFloatStrict1D, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[int]]: ...
    @overload  # n: 0d, p: 2d -> 2d
    def __call__(
        self, /, n: int, p: _ToFloatStrict2D, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[int, int]]: ...
    @overload  # n: 0d, p: 3d -> 3d
    def __call__(
        self, /, n: int, p: _ToFloatStrict3D, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[int, int, int]]: ...
    @overload  # n: Nd, p: ?d -> ?d  (pyright workaround)
    def __call__(
        self, /, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape], seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[Any, ...]]: ...
    @overload  # n: 1d, p: {1,2}d -> 2d
    def __call__(
        self, /, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[int, int]]: ...
    @overload  # n: 1d, p: 3d -> 3d
    def __call__(
        self, /, n: _ToIntStrict1D, p: _ToFloatStrict3D, seed: onp.random.ToRNG | None = None
    ) -> multinomial_frozen[tuple[int, int, int]]: ...
    @overload  # n: 2d, p: {1,2,3}d -> 3d
    def __call__(
        self,
        /,
        n: _ToIntStrict2D,
        p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D,
        seed: onp.random.ToRNG | None = None,
    ) -> multinomial_frozen[tuple[int, int, int]]: ...
    @overload  # ?d, ?d -> ?d  (fallback)
    def __call__(self, /, n: int | _ToIntND, p: _ToFloatND, seed: onp.random.ToRNG | None = None) -> multinomial_frozen: ...

    # unlike `pmf`, `logpmf` returns 0d arrays instead of scalars
    @overload  # x: ?d  (pyright workaround)
    def logpmf(
        self, /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape], n: int | _ToIntND, p: _ToFloatND
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 1d, n: 0d, p: 1d -> 0d
    def logpmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict1D) -> onp.Array0D[np.float64]: ...
    @overload  # x: 1d, n: 0d, p: 2d -> 1d
    def logpmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # x: 1d, n: 0d, p: 3d -> 2d
    def logpmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 1d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def logpmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 1d, n: 1d, p: {1,2}d -> 1d
    def logpmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # x: 1d, n: 1d, p: 3d -> 2d
    def logpmf(self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 1d, n: 2d, p: {1,2,3}d -> 2d
    def logpmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: 0d, p: {1,2}d -> 1d
    def logpmf(self, /, x: onp.ToFloatStrict2D, n: int, p: _ToFloatStrict1D | _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # x: 2d, n: 0d, p: 3d -> 2d
    def logpmf(self, /, x: onp.ToFloatStrict2D, n: int, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def logpmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 2d, n: 1d, p: {1,2}d -> 1d
    def logpmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # x: 2d, n: 1d, p: 3d -> 2d
    def logpmf(self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: 2d, p: {1,2,3}d -> 2d
    def logpmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: 0d, p: {1,2,3}d -> 2d
    def logpmf(
        self, /, x: onp.ToFloatStrict3D, n: int, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def logpmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 3d, n: 1d, p: {1,2,3}d -> 2d
    def logpmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: 2d, p: {1,2,3}d -> 2d
    def logpmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: ?d, ?d, ?d -> ?d  (fallback)
    def logpmf(self, /, x: onp.ToFloatND, n: int | _ToIntND, p: _ToFloatND) -> onp.ArrayND[np.float64]: ...

    # TODO(@jorenham): more shape-type overloads
    @overload  # x: ?d  (pyright workaround)
    def pmf(
        self, /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape], n: int | _ToIntND, p: _ToFloatND
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # x: 1d, n: 0d, p: 1d -> 0d
    def pmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict1D) -> np.float64: ...
    @overload  # x: 1d, n: 0d, p: 2d -> 1d
    def pmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # x: 1d, n: 0d, p: 3d -> 2d
    def pmf(self, /, x: onp.ToFloatStrict1D, n: int, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 1d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def pmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 1d, n: 1d, p: {1,2}d -> 1d
    def pmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # x: 1d, n: 1d, p: 3d -> 2d
    def pmf(self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 1d, n: 2d, p: {1,2,3}d -> 2d
    def pmf(
        self, /, x: onp.ToFloatStrict1D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: 0d, p: {1,2}d -> 1d
    def pmf(self, /, x: onp.ToFloatStrict2D, n: int, p: _ToFloatStrict1D | _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # x: 2d, n: 0d, p: 3d -> 2d
    def pmf(self, /, x: onp.ToFloatStrict2D, n: int, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def pmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 2d, n: 1d, p: {1,2}d -> 1d
    def pmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # x: 2d, n: 1d, p: 3d -> 2d
    def pmf(self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # x: 2d, n: 2d, p: {1,2,3}d -> 2d
    def pmf(
        self, /, x: onp.ToFloatStrict2D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: 0d, p: {1,2,3}d -> 2d
    def pmf(
        self, /, x: onp.ToFloatStrict3D, n: int, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: Nd, p: ?d -> ?d  (pyright workaround)
    def pmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # x: 3d, n: 1d, p: {1,2,3}d -> 2d
    def pmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: 3d, n: 2d, p: {1,2,3}d -> 2d
    def pmf(
        self, /, x: onp.ToFloatStrict3D, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # x: ?d, ?d, ?d -> ?d  (fallback)
    def pmf(self, /, x: onp.ToFloatND, n: int | _ToIntND, p: _ToFloatND) -> onp.ArrayND[np.float64] | Any: ...

    # keep in sync with `__call__`
    @overload  # n: 0d, p: Nd -> Nd
    def mean[ShapeT: tuple[int, ...]](
        self, /, n: int, p: onp.ArrayND[npc.floating, ShapeT]
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload  # n: 0d, p: 1d -> 1d
    def mean(self, /, n: int, p: _ToFloatStrict1D) -> onp.Array1D[np.float64]: ...
    @overload  # n: 0d, p: 2d -> 2d
    def mean(self, /, n: int, p: _ToFloatStrict2D) -> onp.Array2D[np.float64]: ...
    @overload  # n: 0d, p: 3d -> 3d
    def mean(self, /, n: int, p: _ToFloatStrict3D) -> onp.Array3D[np.float64]: ...
    @overload  # n: Nd, p: ?d -> ?d  (pyright workaround)
    def mean(self, /, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]) -> onp.ArrayND[np.float64]: ...
    @overload  # n: 1d, p: {1,2}d -> 2d
    def mean(self, /, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D) -> onp.Array2D[np.float64]: ...
    @overload  # n: 1d, p: 3d -> 3d
    def mean(self, /, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array3D[np.float64]: ...
    @overload  # n: 2d, p: {1,2,3}d -> 3d
    def mean(
        self, /, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array3D[np.float64]: ...
    @overload  # ?d, ?d -> ?d  (fallback)
    def mean(self, /, n: int | _ToIntND, p: _ToFloatND) -> onp.ArrayND[np.float64, _WorkaroundForPyright]: ...

    #
    @overload  # n: 0d, p: ?d -> ?d  (pyright workaround)
    def cov(self, /, n: int, p: onp.ArrayND[npc.floating, _JustAnyShape]) -> onp.ArrayND[np.float64]: ...
    @overload  # n: 0d, p: 1d -> 2d
    def cov(self, /, n: int, p: _ToFloatStrict1D) -> onp.Array2D[np.float64]: ...
    @overload  # n: 0d, p: 2d -> 3d
    def cov(self, /, n: int, p: _ToFloatStrict2D) -> onp.Array3D[np.float64]: ...
    @overload  # n: 0d, p: 3d -> 4d
    def cov(self, /, n: int, p: _ToFloatStrict3D) -> onp.ArrayND[np.float64, tuple[int, int, int, int]]: ...
    @overload  # n: ?d, p: Nd -> ?d  (pyright workaround)
    def cov(self, /, n: onp.ArrayND[npc.integer, _JustAnyShape], p: _ToFloatND) -> onp.ArrayND[np.float64]: ...
    @overload  # n: Nd, p: ?d -> ?d  (pyright workaround)
    def cov(self, /, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]) -> onp.ArrayND[np.float64]: ...
    @overload  # n: 1d, p: {1,2}d -> 3d
    def cov(self, /, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D) -> onp.Array3D[np.float64]: ...
    @overload  # n: 1d, p: 3d -> 3d
    def cov(self, /, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.ArrayND[np.float64, tuple[int, int, int, int]]: ...
    @overload  # n: 2d, p: {1,2,3}d -> 4d
    def cov(
        self, /, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.ArrayND[np.float64, tuple[int, int, int, int]]: ...
    @overload  # ?d, ?d -> ?d  (fallback)
    def cov(self, /, n: int | _ToIntND, p: _ToFloatND) -> onp.ArrayND[np.float64]: ...

    # keep in sync with `cov` (but return `ndim-2`)
    # `entropy` returns 0d arrays instead of bare scalars.
    @overload  # n: 0d, p: ?d -> ?d  (pyright workaround)
    def entropy(self, /, n: int, p: onp.ArrayND[npc.floating, _JustAnyShape]) -> onp.ArrayND[np.float64]: ...
    @overload  # n: 0d, p: 1d -> 0d
    def entropy(self, /, n: int, p: _ToFloatStrict1D) -> onp.Array0D[np.float64]: ...
    @overload  # n: 0d, p: 2d -> 1d
    def entropy(self, /, n: int, p: _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # n: 0d, p: 3d -> 2d
    def entropy(self, /, n: int, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # n: ?d, p: Nd -> ?d  (pyright workaround)
    def entropy(self, /, n: onp.ArrayND[npc.integer, _JustAnyShape], p: _ToFloatND) -> onp.ArrayND[np.float64]: ...
    @overload  # n: Nd, p: ?d -> ?d  (pyright workaround)
    def entropy(self, /, n: _ToIntND, p: onp.ArrayND[npc.floating, _JustAnyShape]) -> onp.ArrayND[np.float64]: ...
    @overload  # n: 1d, p: {1,2}d -> 1d
    def entropy(self, /, n: _ToIntStrict1D, p: _ToFloatStrict1D | _ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # n: 1d, p: 3d -> 2d
    def entropy(self, /, n: _ToIntStrict1D, p: _ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # n: 2d, p: {1,2,3}d -> 2d
    def entropy(
        self, /, n: _ToIntStrict2D, p: _ToFloatStrict1D | _ToFloatStrict2D | _ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # ?d, ?d -> ?d  (fallback)
    def entropy(self, /, n: int | _ToIntND, p: _ToFloatND) -> onp.ArrayND[np.float64]: ...

    # this will currently (1.17.1) raise a `TypeError` for `ndim(n)>0` or `ndim(p)>1`
    @overload
    def rvs(
        self, /, n: int, p: _ToFloatStrict1D, size: tuple[()] | None = None, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array1D[np.int_]: ...
    @overload
    def rvs(
        self, /, n: int, p: _ToFloatStrict1D, size: int | tuple[int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.int_]: ...
    @overload
    def rvs(
        self, /, n: int, p: _ToFloatStrict1D, size: tuple[int, int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.int_]: ...
    @overload
    def rvs(
        self, /, n: int, p: _ToFloatStrict1D, size: tuple[int, ...], random_state: onp.random.ToRNG | None = None
    ) -> onp.ArrayND[np.int_]: ...

# `_ShapeT_co` corresponds to the shape of the mean
class multinomial_frozen(multi_rv_frozen[multinomial_gen], Generic[_ShapeT_co]):
    def __init__(self, /, n: int | _ToIntND, p: _ToFloatND, seed: onp.random.ToRNG | None = None) -> None: ...

    # unlike `pmf`, `logpmf` returns 0d arrays instead of scalars
    @overload  # self: ?d, x: 1d -> >=0d  (pyright workaround)
    def logpmf(self: multinomial_frozen[_JustAnyShape], /, x: onp.ToFloatStrict1D) -> onp.ArrayND[np.float64]: ...
    @overload  # self: 1d, x: ?d -> >=0d  (pyright workaround)
    def logpmf(
        self: multinomial_frozen[tuple[int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # self: 1d, x: 1d -> 0d
    def logpmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict1D) -> onp.Array0D[np.float64]: ...
    @overload  # self: 1d, x: 2d -> 1d
    def logpmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # self: 1d, x: 3d -> 2d
    def logpmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # self: 2d, x: ?d -> >=1d  (pyright workaround)
    def logpmf(
        self: multinomial_frozen[tuple[int, int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # self: 2d, x: {1,2}d  -> 1d
    def logpmf(
        self: multinomial_frozen[tuple[int, int]], /, x: onp.ToFloatStrict1D | onp.ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # self: 2d, x: 3d  -> 2d
    def logpmf(self: multinomial_frozen[tuple[int, int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # self: 3d, x: ?d -> >=2d  (pyright workaround)
    def logpmf(
        self: multinomial_frozen[tuple[int, int, int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # self: 3d, x: {1,2,3}d  -> 2d
    def logpmf(
        self: multinomial_frozen[tuple[int, int, int]], /, x: onp.ToFloatStrict1D | onp.ToFloatStrict2D | onp.ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # fallback
    def logpmf(self, /, x: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

    #
    @overload  # self: ?d, x: 1d -> >=0d  (pyright workaround)
    def pmf(self: multinomial_frozen[_JustAnyShape], /, x: onp.ToFloatStrict1D) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # self: 1d, x: ?d -> >=0d  (pyright workaround)
    def pmf(
        self: multinomial_frozen[tuple[int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload  # self: 1d, x: 1d -> 0d
    def pmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict1D) -> np.float64: ...
    @overload  # self: 1d, x: 2d -> 1d
    def pmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict2D) -> onp.Array1D[np.float64]: ...
    @overload  # self: 1d, x: 3d -> 2d
    def pmf(self: multinomial_frozen[tuple[int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # self: 2d, x: ?d -> >=1d  (pyright workaround)
    def pmf(
        self: multinomial_frozen[tuple[int, int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> onp.ArrayND[np.float64, onp.AtLeast1D[Any]]: ...
    @overload  # self: 2d, x: {1,2}d  -> 1d
    def pmf(
        self: multinomial_frozen[tuple[int, int]], /, x: onp.ToFloatStrict1D | onp.ToFloatStrict2D
    ) -> onp.Array1D[np.float64]: ...
    @overload  # self: 2d, x: 3d  -> 2d
    def pmf(self: multinomial_frozen[tuple[int, int]], /, x: onp.ToFloatStrict3D) -> onp.Array2D[np.float64]: ...
    @overload  # self: 3d, x: ?d -> >=2d  (pyright workaround)
    def pmf(
        self: multinomial_frozen[tuple[int, int, int]], /, x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]
    ) -> onp.ArrayND[np.float64, onp.AtLeast2D[Any]]: ...
    @overload  # self: 3d, x: {1,2,3}d  -> 2d
    def pmf(
        self: multinomial_frozen[tuple[int, int, int]], /, x: onp.ToFloatStrict1D | onp.ToFloatStrict2D | onp.ToFloatStrict3D
    ) -> onp.Array2D[np.float64]: ...
    @overload  # fallback
    def pmf(self, /, x: onp.ToFloatND) -> onp.ArrayND[np.float64] | Any: ...

    #
    def mean(self, /) -> onp.ArrayND[np.float64, _ShapeT_co]: ...

    #
    @overload  # ?d -> >=2d  (pyright workaround)
    def cov(self: multinomial_frozen[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload  # 1d -> 2d
    def cov(self: multinomial_frozen[tuple[int]], /) -> onp.Array2D[np.float64]: ...
    @overload  # 2d -> 3d
    def cov(self: multinomial_frozen[tuple[int, int]], /) -> onp.Array3D[np.float64]: ...
    @overload  # 3d -> 4d  (`optype.numpy` has no `Array4D` at the moment)
    def cov(self: multinomial_frozen[tuple[int, int, int]], /) -> onp.ArrayND[np.float64, tuple[int, int, int, int]]: ...
    @overload  # fallback
    def cov(self, /) -> onp.ArrayND[np.float64]: ...

    # `entropy` returns 0d arrays instead of bare scalars.
    @overload  # ?d -> >=0d  (pyright workaround)
    def entropy(self: multinomial_frozen[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload  # 1d -> 0d
    def entropy(self: multinomial_frozen[tuple[int]], /) -> onp.Array0D[np.float64]: ...
    @overload  # 2d -> 1d
    def entropy(self: multinomial_frozen[tuple[int, int]], /) -> onp.Array1D[np.float64]: ...
    @overload  # 3d -> 2d
    def entropy(self: multinomial_frozen[tuple[int, int, int]], /) -> onp.Array2D[np.float64]: ...
    @overload  # fallback
    def entropy(self, /) -> onp.ArrayND[np.float64]: ...

    # this will currently (1.17.1) raise a `TypeError` for >1d `_ShapeT_co`
    @overload
    def rvs(
        self: multinomial_frozen[tuple[int]], /, size: tuple[()], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array1D[np.int_]: ...
    @overload
    def rvs(
        self: multinomial_frozen[tuple[int]], /, size: int | tuple[int] = 1, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array2D[np.int_]: ...
    @overload
    def rvs(
        self: multinomial_frozen[tuple[int]], /, size: tuple[int, int], random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.int_]: ...
    @overload
    def rvs(
        self: multinomial_frozen[tuple[int]], /, size: tuple[int, ...], random_state: onp.random.ToRNG | None = None
    ) -> onp.ArrayND[np.int_]: ...

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

# pyrefly: ignore [inconsistent-inheritance]
class special_ortho_group_frozen(_group_rv_frozen_mixin, multi_rv_frozen[special_ortho_group_gen]): ...  # type: ignore[misc]
class ortho_group_gen(_group_rv_gen_mixin[ortho_group_frozen], multi_rv_generic): ...

# pyrefly: ignore [inconsistent-inheritance]
class ortho_group_frozen(_group_rv_frozen_mixin, multi_rv_frozen[ortho_group_gen]): ...  # type: ignore[misc]
class unitary_group_gen(_group_rv_gen_mixin[unitary_group_frozen, np.complex128], multi_rv_generic): ...

# pyrefly: ignore [inconsistent-inheritance]
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
    # pyrefly: ignore [bad-override]
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
    @overload  # shape=<0d>, loc=<0d>, size=1 (default)
    def rvs(
        self,
        /,
        loc: onp.ToFloat | None = None,
        shape: onp.ToFloat = 1,
        df: int = 1,
        size: Literal[1] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> np.float64: ...
    @overload  # shape=<2d> (positional)
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None,
        shape: onp.ToFloat2D,
        df: int = 1,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # shape=<2d> (keyword)
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        *,
        shape: onp.ToFloat2D,
        df: int = 1,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # size=(int, ...)
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        *,
        size: tuple[int, ...],
        random_state: onp.random.ToRNG | None = None,
    ) -> onp.ArrayND[np.float64]: ...
    @overload  # fallback
    def rvs(
        self,
        /,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        size: int | tuple[int, ...] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    def marginal(
        self,
        dimensions: int | onp.ToInt1D,
        loc: onp.ToFloat1D | None = None,
        shape: onp.ToFloat | onp.ToFloat2D = 1,
        df: int = 1,
        allow_singular: bool = False,
    ) -> multivariate_t_frozen: ...

# TODO: make generic: https://github.com/scipy/scipy-stubs/issues/406
class multivariate_t_frozen(multi_rv_frozen[multivariate_t_gen]):
    # pyrefly: ignore [bad-override]
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
    def entropy(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def rvs(self, /, size: int = 1, random_state: onp.random.ToRNG | None = None) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def rvs(self, /, size: tuple[int, ...], random_state: onp.random.ToRNG | None = None) -> onp.ArrayND[np.float64]: ...

    #
    def marginal(self, dimensions: int | onp.ToInt1D) -> multivariate_t_frozen: ...

# NOTE: `m` and `n` are broadcastable (but doing so will break `.rvs()` at runtime...)
class multivariate_hypergeom_gen(multi_rv_generic):
    def __call__(
        self, /, m: onp.ToJustIntND, n: onp.ToJustInt | onp.ToJustIntND, seed: onp.random.ToRNG | None = None
    ) -> multivariate_hypergeom_frozen: ...
    def logpmf(self, /, x: onp.ToFloatND, m: onp.ToJustIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array1ND[np.float64]: ...
    def var(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array1ND[np.float64]: ...
    def cov(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND) -> _Array2ND[np.float64]: ...
    def rvs(
        self,
        /,
        m: onp.ToIntND,
        n: onp.ToJustInt | onp.ToJustIntND,
        size: int | tuple[int, ...] | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array1ND[np.int_]: ...

class multivariate_hypergeom_frozen(multi_rv_frozen[multivariate_hypergeom_gen]):
    def __init__(self, /, m: onp.ToIntND, n: onp.ToJustInt | onp.ToJustIntND, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> _Array1ND[np.float64]: ...
    def var(self, /) -> _Array1ND[np.float64]: ...
    def cov(self, /) -> _Array2ND[np.float64]: ...

    #
    @overload
    def rvs(self, /, size: tuple[()], random_state: onp.random.ToRNG | None = None) -> _Array1ND[np.int_]: ...
    @overload
    def rvs(self, /, size: int | tuple[int] = 1, random_state: onp.random.ToRNG | None = None) -> _Array2ND[np.int_]: ...

type _RandomTableRVSMethod = Literal["boyett", "patefield"]

class random_table_gen(multi_rv_generic):
    def __call__(
        self, /, row: onp.ToJustIntND, col: onp.ToJustIntND, *, seed: onp.random.ToRNG | None = None
    ) -> random_table_frozen: ...

    #
    def logpmf(self, /, x: onp.ToFloatND, row: onp.ToJustIntND, col: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToFloatND, row: onp.ToJustIntND, col: onp.ToJustIntND) -> _ScalarOrArray_f8: ...

    #
    def mean(self, /, row: onp.ToJustIntND, col: onp.ToJustIntND) -> onp.Array2D[np.float64]: ...

    # NOTE: the docstring example incorrectly suggest a float64 return dtype
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
    ) -> onp.Array2D[np.int_]: ...
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
    ) -> onp.Array3D[np.int_]: ...
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
    ) -> _Array3ND[np.int_]: ...

class random_table_frozen(multi_rv_frozen[random_table_gen]):
    # pyrefly: ignore [bad-override]
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
    ) -> onp.Array2D[np.int_]: ...
    @overload
    def rvs(
        self, /, size: int | tuple[int], method: _RandomTableRVSMethod | None = None, random_state: onp.random.ToRNG | None = None
    ) -> onp.Array3D[np.int_]: ...
    @overload
    def rvs(
        self, /, size: onp.AtLeast1D, method: _RandomTableRVSMethod | None = None, random_state: onp.random.ToRNG | None = None
    ) -> _Array3ND[np.int_]: ...

class dirichlet_multinomial_gen(multi_rv_generic):
    def __call__(
        self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND, seed: onp.random.ToRNG | None = None
    ) -> dirichlet_multinomial_frozen: ...
    def logpmf(self, /, x: onp.ToIntND, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToIntND, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array1ND[np.float64]: ...
    def var(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array1ND[np.float64]: ...
    def cov(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND) -> _Array2ND[np.float64]: ...

class dirichlet_multinomial_frozen(multi_rv_frozen[dirichlet_multinomial_gen]):
    alpha: _Array1ND[np.float64]
    n: _Array1ND[np.int_]  # broadcasted against alpha

    def __init__(self, /, alpha: onp.ToFloatND, n: onp.ToJustIntND, seed: onp.random.ToRNG | None = None) -> None: ...
    def logpmf(self, /, x: onp.ToIntND) -> _ScalarOrArray_f8: ...
    def pmf(self, /, x: onp.ToIntND) -> _ScalarOrArray_f8: ...
    def mean(self, /) -> _Array1ND[np.float64]: ...
    def var(self, /) -> _Array1ND[np.float64]: ...
    def cov(self, /) -> _Array2ND[np.float64]: ...

class vonmises_fisher_gen(multi_rv_generic):
    def __call__(
        self, /, mu: onp.ToFloat1D | None = None, kappa: onp.ToFloat = 1, seed: onp.random.ToRNG | None = None
    ) -> vonmises_fisher_frozen: ...
    def logpdf(self, /, x: onp.ToFloatND, mu: onp.ToFloat1D | None = None, kappa: onp.ToFloat = 1) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND, mu: onp.ToFloat1D | None = None, kappa: onp.ToFloat = 1) -> _ScalarOrArray_f8: ...
    def entropy(self, /, mu: onp.ToFloat1D | None = None, kappa: onp.ToFloat = 1) -> np.float64: ...
    def rvs(
        self,
        /,
        mu: onp.ToFloat1D | None = None,
        kappa: onp.ToFloat = 1,
        size: SupportsIndex | tuple[SupportsIndex, *tuple[SupportsIndex, ...]] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...
    def fit(self, /, x: onp.ToFloatND) -> tuple[onp.Array1D[np.float64], float]: ...

class vonmises_fisher_frozen(multi_rv_frozen[vonmises_fisher_gen]):
    def __init__(
        self, /, mu: onp.ToFloat1D | None = None, kappa: onp.ToFloat = 1, seed: onp.random.ToRNG | None = None
    ) -> None: ...
    def logpdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def pdf(self, /, x: onp.ToFloatND) -> _ScalarOrArray_f8: ...
    def entropy(self, /) -> np.float64: ...
    def rvs(
        self,
        /,
        size: SupportsIndex | tuple[SupportsIndex, *tuple[SupportsIndex, ...]] = 1,
        random_state: onp.random.ToRNG | None = None,
    ) -> _Array2ND[np.float64]: ...

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
        size: SupportsIndex | tuple[SupportsIndex, ...] | None = None,
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
        self, /, size: SupportsIndex | tuple[SupportsIndex, ...] | None = None, random_state: onp.random.ToRNG | None = None
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
