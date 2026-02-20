from typing import ClassVar, Final, Literal, final

from ._biasedurn import _PyFishersNCHypergeometric, _PyWalleniusNCHypergeometric
from ._distn_infrastructure import rv_discrete

__all__ = [
    "bernoulli",
    "bernoulli_gen",
    "betabinom",
    "betabinom_gen",
    "betanbinom",
    "betanbinom_gen",
    "binom",
    "binom_gen",
    "boltzmann",
    "boltzmann_gen",
    "dlaplace",
    "dlaplace_gen",
    "geom",
    "geom_gen",
    "hypergeom",
    "hypergeom_gen",
    "logser",
    "logser_gen",
    "nbinom",
    "nbinom_gen",
    "nchypergeom_fisher",
    "nchypergeom_fisher_gen",
    "nchypergeom_wallenius",
    "nchypergeom_wallenius_gen",
    "nhypergeom",
    "nhypergeom_gen",
    "planck",
    "planck_gen",
    "poisson",
    "poisson_binom",
    "poisson_binom_gen",
    "poisson_gen",
    "randint",
    "randint_gen",
    "skellam",
    "skellam_gen",
    "yulesimon",
    "yulesimon_gen",
    "zipf",
    "zipf_gen",
    "zipfian",
    "zipfian_gen",
]

class binom_gen(rv_discrete): ...
class bernoulli_gen(binom_gen): ...
class betabinom_gen(rv_discrete): ...
class nbinom_gen(rv_discrete): ...
class betanbinom_gen(rv_discrete): ...
class geom_gen(rv_discrete): ...
class hypergeom_gen(rv_discrete): ...
class nhypergeom_gen(rv_discrete): ...
class logser_gen(rv_discrete): ...
class poisson_gen(rv_discrete): ...
class planck_gen(rv_discrete): ...
class boltzmann_gen(rv_discrete): ...
class randint_gen(rv_discrete): ...
class zipf_gen(rv_discrete): ...
class zipfian_gen(rv_discrete): ...
class dlaplace_gen(rv_discrete): ...
class poisson_binom_gen(rv_discrete): ...
class skellam_gen(rv_discrete): ...
class yulesimon_gen(rv_discrete): ...

class _nchypergeom_gen(rv_discrete):
    rvs_name: ClassVar[Literal["rvs_fisher", "rvs_wallenius"] | None] = None
    dist: ClassVar[type[_PyFishersNCHypergeometric | _PyWalleniusNCHypergeometric] | None] = None

@final
class nchypergeom_fisher_gen(_nchypergeom_gen):
    rvs_name: ClassVar[Literal["rvs_fisher"]] = "rvs_fisher"  # pyright: ignore[reportIncompatibleVariableOverride]
    dist: ClassVar[type[_PyFishersNCHypergeometric]] = ...  # pyright: ignore[reportIncompatibleVariableOverride]

@final
class nchypergeom_wallenius_gen(_nchypergeom_gen):
    rvs_name: ClassVar[Literal["rvs_wallenius"]] = "rvs_wallenius"  # pyright: ignore[reportIncompatibleVariableOverride]
    dist: ClassVar[type[_PyWalleniusNCHypergeometric]] = ...  # pyright: ignore[reportIncompatibleVariableOverride]

binom: Final[binom_gen]
bernoulli: Final[bernoulli_gen]
betabinom: Final[betabinom_gen]
nbinom: Final[nbinom_gen]
betanbinom: Final[betanbinom_gen]
geom: Final[geom_gen]
hypergeom: Final[hypergeom_gen]
nhypergeom: Final[nhypergeom_gen]
logser: Final[logser_gen]
poisson: Final[poisson_gen]
planck: Final[planck_gen]
boltzmann: Final[boltzmann_gen]
randint: Final[randint_gen]
zipf: Final[zipf_gen]
zipfian: Final[zipfian_gen]
dlaplace: Final[dlaplace_gen]
poisson_binom: Final[poisson_binom_gen]
skellam: Final[skellam_gen]
yulesimon: Final[yulesimon_gen]
nchypergeom_fisher: Final[nchypergeom_fisher_gen]
nchypergeom_wallenius: Final[nchypergeom_wallenius_gen]
