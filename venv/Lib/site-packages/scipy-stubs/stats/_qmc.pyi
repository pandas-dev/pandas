import abc
import numbers
from collections.abc import Callable, Mapping
from typing import (
    ClassVar,
    Concatenate,
    Final,
    Generic,
    Literal,
    Protocol,
    Self,
    SupportsIndex,
    overload,
    override,
    type_check_only,
)
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

from scipy.spatial.distance import _Metric

__all__ = [
    "Halton",
    "LatinHypercube",
    "MultinomialQMC",
    "MultivariateNormalQMC",
    "PoissonDisk",
    "QMCEngine",
    "Sobol",
    "discrepancy",
    "geometric_discrepancy",
    "scale",
    "update_discrepancy",
]

###

type _Real = npc.floating | npc.integer

type _MethodQMC = Literal["random-cd", "lloyd"]
type _MethodDisc = Literal["CD", "WD", "MD", "L2-star"]
type _MethodDist = Literal["mindist", "mst"]
type _HyperSphere = Literal["volume", "surface"]

_InexactT_co = TypeVar("_InexactT_co", bound=npc.inexact, default=np.float64, covariant=True)
_EngineT_co = TypeVar("_EngineT_co", bound=QMCEngine[npc.inexact], default=Sobol, covariant=True)

@type_check_only
class _Optimizer(Protocol):
    def __call__[InexactT: npc.inexact](self, sample: onp.ArrayND[InexactT], /) -> onp.Array2D[InexactT]: ...

@type_check_only
class _HypersphereMethod(Protocol):
    def __call__(
        self, /, center: onp.Array2D[np.float64], radius: onp.ToFloat, candidates: onp.ToJustInt = 1
    ) -> onp.Array2D[np.float64]: ...

@type_check_only
class _QMCDistribution(Generic[_EngineT_co]):
    engine: _EngineT_co

    def __init__(self, /, *, engine: _EngineT_co | None = None, rng: onp.random.ToRNG | None = None) -> None: ...
    def random[InexactT: npc.inexact](
        self: _QMCDistribution[QMCEngine[InexactT]], /, n: onp.ToJustInt = 1
    ) -> onp.Array2D[InexactT]: ...

###

class QMCEngine(abc.ABC, Generic[_InexactT_co]):
    d: Final[int | np.int32 | np.int64]
    optimization_method: Callable[Concatenate[onp.ArrayND[_InexactT_co], ...], onp.Array2D[_InexactT_co]]
    rng_seed: Final[np.random.Generator]

    rng: np.random.Generator
    num_generated: int

    @overload
    @abc.abstractmethod
    def __init__(
        self, /, d: onp.ToJustInt, *, optimization: _MethodQMC | None = None, rng: onp.random.ToRNG | None = None
    ) -> None: ...
    @overload  # will be deprecated in the future
    @abc.abstractmethod
    def __init__(self, /, d: onp.ToJustInt, *, optimization: _MethodQMC | None = None, seed: onp.random.ToRNG | None) -> None: ...

    #
    def _initialize(
        self, /, d: onp.ToJustInt, *, optimization: _MethodQMC | None = None, rng: onp.random.ToRNG | None = None
    ) -> None: ...

    #
    @abc.abstractmethod
    def _random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[_InexactT_co]: ...
    def random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[_InexactT_co]: ...
    def integers(
        self,
        /,
        l_bounds: onp.ToInt | onp.ToIntND,
        *,
        u_bounds: onp.ToInt | onp.ToIntND | None = None,
        n: onp.ToJustInt = 1,
        endpoint: bool = False,
        workers: onp.ToJustInt = 1,
    ) -> onp.Array2D[np.int64]: ...

    #
    def reset(self, /) -> Self: ...
    def fast_forward(self, /, n: onp.ToJustInt) -> Self: ...

class Halton(QMCEngine[np.float64]):
    base: Final[list[int]]
    scramble: Final[bool]
    _permutations: Final[list[onp.Array2D[np.int_]]]

    @overload
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        scramble: bool = True,
        optimization: _MethodQMC | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self, /, d: onp.ToJustInt, *, scramble: bool = True, optimization: _MethodQMC | None = None, seed: onp.random.ToRNG | None
    ) -> None: ...

    #
    def _initialize_permutations(self, /) -> None: ...
    @override
    def _random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[np.float64]: ...

class LatinHypercube(QMCEngine[np.float64]):
    scramble: Final[bool]
    lhs_method: Final[Callable[[onp.ToInt], onp.Array2D[np.float64]]]

    @overload
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        scramble: bool = True,
        strength: int = 1,
        optimization: _MethodQMC | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        scramble: bool = True,
        strength: int = 1,
        optimization: _MethodQMC | None = None,
        seed: onp.random.ToRNG | None,
    ) -> None: ...

    #
    @override
    def _random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[np.float64]: ...
    def _random_lhs(self, /, n: onp.ToJustInt = 1) -> onp.Array2D[np.float64]: ...
    def _random_oa_lhs(self, /, n: onp.ToJustInt = 4) -> onp.Array2D[np.float64]: ...

class Sobol(QMCEngine[np.float64]):
    MAXDIM: ClassVar[int] = 21_201

    dtype_i: Final[type[np.uint32 | np.uint64]]
    scramble: Final[bool]

    bits: int | npc.integer
    maxn: int | npc.integer

    @overload
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        scramble: bool = True,
        bits: onp.ToJustInt | None = None,
        optimization: _MethodQMC | None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        scramble: bool = True,
        bits: onp.ToJustInt | None = None,
        optimization: _MethodQMC | None = None,
        seed: onp.random.ToRNG | None,
    ) -> None: ...

    #
    def _scramble(self, /) -> None: ...
    @override
    def _random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[np.float64]: ...
    def random_base2(self, /, m: onp.ToJustInt) -> onp.Array2D[np.float64]: ...

class PoissonDisk(QMCEngine[np.float64]):
    hypersphere_method: Final[_HypersphereMethod]
    radius_factor: Final[float]
    radius: Final[onp.ToFloat]
    radius_squared: Final[onp.ToFloat]
    ncandidates: Final[onp.ToInt]

    l_bounds: onp.Array1D[np.float64]
    u_bounds: onp.Array1D[np.float64]

    cell_size: Final[np.float64]
    grid_size: Final[onp.Array1D[np.int_]]

    sample_pool: list[onp.Array1D]
    sample_grid: onp.Array2D[np.float32]

    @overload
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        radius: onp.ToFloat = 0.05,
        hypersphere: _HyperSphere = "volume",
        ncandidates: onp.ToJustInt = 30,
        optimization: _MethodQMC | None = None,
        rng: onp.random.ToRNG | None = None,
        l_bounds: onp.ToFloat1D | None = None,
        u_bounds: onp.ToFloat1D | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self,
        /,
        d: onp.ToJustInt,
        *,
        radius: onp.ToFloat = 0.05,
        hypersphere: _HyperSphere = "volume",
        ncandidates: onp.ToJustInt = 30,
        optimization: _MethodQMC | None = None,
        seed: onp.random.ToRNG | None,
        l_bounds: onp.ToFloat1D | None = None,
        u_bounds: onp.ToFloat1D | None = None,
    ) -> None: ...

    #
    def _initialize_grid_pool(self, /) -> None: ...
    @override
    def _random(self, /, n: onp.ToJustInt = 1, *, workers: onp.ToJustInt = 1) -> onp.Array2D[np.float64]: ...

    #
    def _hypersphere_volume_sample(
        self, /, center: onp.Array2D[np.float64], radius: onp.ToFloat, candidates: onp.ToJustInt = 1
    ) -> onp.Array2D[np.float64]: ...
    def _hypersphere_surface_sample(
        self, /, center: onp.Array2D[np.float64], radius: onp.ToFloat, candidates: onp.ToJustInt = 1
    ) -> onp.Array2D[np.float64]: ...

    #
    def fill_space(self, /) -> onp.Array2D[np.float64]: ...

class MultivariateNormalQMC(_QMCDistribution[_EngineT_co], Generic[_EngineT_co]):
    @overload
    def __init__(
        self: MultivariateNormalQMC[Sobol],
        /,
        mean: onp.ToFloat1D,
        cov: onp.ToFloat2D | None = None,
        *,
        cov_root: onp.ToFloat2D | None = None,
        inv_transform: bool = True,
        engine: None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        mean: onp.ToFloat1D,
        cov: onp.ToFloat2D | None = None,
        *,
        cov_root: onp.ToFloat2D | None = None,
        inv_transform: bool = True,
        engine: _EngineT_co,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self: MultivariateNormalQMC[Sobol],
        /,
        mean: onp.ToFloat1D,
        cov: onp.ToFloat2D | None = None,
        *,
        cov_root: onp.ToFloat2D | None = None,
        inv_transform: bool = True,
        engine: None = None,
        seed: onp.random.ToRNG | None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self,
        /,
        mean: onp.ToFloat1D,
        cov: onp.ToFloat2D | None = None,
        *,
        cov_root: onp.ToFloat2D | None = None,
        inv_transform: bool = True,
        engine: _EngineT_co,
        seed: onp.random.ToRNG | None,
    ) -> None: ...
    #
    def _correlate[InexactT: npc.inexact](
        self: MultivariateNormalQMC[QMCEngine[InexactT]], /, base_samples: onp.Array2D[InexactT]
    ) -> onp.Array2D[InexactT]: ...
    def _standard_normal_samples[InexactT: npc.inexact](
        self: MultivariateNormalQMC[QMCEngine[InexactT]], /, n: onp.ToJustInt = 1
    ) -> onp.Array2D[InexactT]: ...

class MultinomialQMC(_QMCDistribution[_EngineT_co], Generic[_EngineT_co]):
    pvals: Final[onp.Array1D[np.float32 | np.float64]]
    n_trials: Final[int | npc.integer]

    @overload
    def __init__(
        self: MultinomialQMC[Sobol],
        /,
        pvals: onp.ToJustFloat | onp.ToJustFloat1D,
        n_trials: onp.ToJustInt,
        *,
        engine: None = None,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        pvals: onp.ToJustFloat | onp.ToJustFloat1D,
        n_trials: onp.ToJustInt,
        *,
        engine: _EngineT_co,
        rng: onp.random.ToRNG | None = None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self: MultinomialQMC[Sobol],
        /,
        pvals: onp.ToJustFloat | onp.ToJustFloat1D,
        n_trials: onp.ToJustInt,
        *,
        engine: None = None,
        seed: onp.random.ToRNG | None,
    ) -> None: ...
    @overload  # will be deprecated in the future
    def __init__(
        self,
        /,
        pvals: onp.ToJustFloat | onp.ToJustFloat1D,
        n_trials: onp.ToJustInt,
        *,
        engine: _EngineT_co,
        seed: onp.random.ToRNG | None,
    ) -> None: ...

#
@overload
def check_random_state(seed: int | npc.integer | numbers.Integral | None = None) -> np.random.Generator: ...
@overload
def check_random_state[AnyRNG: (np.random.Generator, np.random.RandomState)](seed: AnyRNG) -> AnyRNG: ...

#
def scale(
    sample: onp.ToFloat2D, l_bounds: onp.ToFloat1D | onp.ToFloat, u_bounds: onp.ToFloat1D | onp.ToFloat, *, reverse: bool = False
) -> onp.Array2D[np.float64]: ...

#
def discrepancy(
    sample: onp.ToFloat2D, *, iterative: bool = False, method: _MethodDisc = "CD", workers: onp.ToJustInt = 1
) -> float | np.float64: ...

#
def geometric_discrepancy(
    sample: onp.ToJustFloat2D, method: _MethodDist = "mindist", metric: _Metric = "euclidean"
) -> float | np.float64: ...

#
def update_discrepancy(x_new: onp.ToJustFloat1D, sample: onp.ToJustFloat2D, initial_disc: opt.AnyFloat) -> float: ...
def primes_from_2_to(n: onp.ToInt) -> onp.Array1D[np.int_]: ...
def n_primes(n: onp.ToInt) -> list[int] | onp.Array1D[np.int_]: ...

#
def _select_optimizer(optimization: _MethodQMC | None, config: Mapping[str, object]) -> _Optimizer | None: ...
def _random_cd[FloatArrayT: onp.ArrayND[npc.floating]](
    best_sample: FloatArrayT, n_iters: onp.ToInt, n_nochange: onp.ToInt, rng: onp.random.RNG
) -> FloatArrayT: ...
def _l1_norm(sample: onp.ToJustFloat2D) -> float | np.float64: ...
def _lloyd_iteration[FloatArrayT: onp.ArrayND[npc.floating]](
    sample: FloatArrayT, decay: onp.ToFloat, qhull_options: str | None
) -> FloatArrayT: ...
def _lloyd_centroidal_voronoi_tessellation(
    sample: onp.ToJustFloat2D, *, tol: onp.ToFloat = 1e-5, maxiter: onp.ToJustInt = 10, qhull_options: str | None = None
) -> onp.Array2D[np.float64]: ...
def _ensure_in_unit_hypercube(sample: onp.ToJustFloat2D) -> onp.Array2D[np.float64]: ...

#
@overload
def _perturb_discrepancy(
    sample: onp.Array2D[npc.integer | np.bool], i1: SupportsIndex, i2: SupportsIndex, k: SupportsIndex, disc: onp.ToFloat
) -> float | np.float64: ...
@overload
def _perturb_discrepancy[InexactT: npc.inexact](
    sample: onp.Array2D[InexactT], i1: SupportsIndex, i2: SupportsIndex, k: SupportsIndex, disc: onp.ToFloat
) -> InexactT: ...

#
def van_der_corput(
    n: op.CanInt,
    base: onp.ToJustInt = 2,
    *,
    start_index: onp.ToJustInt = 0,
    scramble: bool = False,
    permutations: onp.ToInt | onp.ToIntND | None = None,
    rng: onp.random.ToRNG | None = None,
    workers: onp.ToJustInt = 1,
) -> onp.Array1D: ...

#
@overload
def _van_der_corput_permutation(base: SupportsIndex, *, random_state: onp.random.ToRNG | None = None) -> onp.Array2D[np.int_]: ...
@overload
def _van_der_corput_permutation(
    base: op.CanFloat, *, random_state: onp.random.ToRNG | None = None
) -> onp.Array2D[np.float64]: ...

#
def _validate_workers(workers: onp.ToJustInt = 1) -> int: ...
def _validate_bounds(
    l_bounds: onp.ToFloat1D, u_bounds: onp.ToFloat1D, d: onp.ToJustInt
) -> tuple[onp.Array1D[_Real], onp.Array1D[_Real]]: ...
