# see `scipy/stats/_biasedurn.pyx`

from typing import Never, Self, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp

@type_check_only
class _PyNCHypergeometric:
    def __new__(cls, /, n: op.CanInt, m: op.CanInt, N: op.CanInt, odds: op.CanFloat, accuracy: op.CanFloat) -> Self: ...
    def mode(self, /) -> int: ...
    def mean(self, /) -> float: ...
    def variance(self, /) -> float: ...
    def probability(self, /, x: op.CanInt) -> float: ...
    def moments(self, /) -> tuple[float, float]: ...

###

# NOTE: These apprear to be broken, and will always raise `TypeError: no default __reduce__ due to non-trivial __cinit__`
def __setstate_cython__(self: Never, pyx_state: Never, /) -> None: ...  # undocumented
def __reduce_cython__(self: Never, /) -> Never: ...  # undocumented

class _PyFishersNCHypergeometric(_PyNCHypergeometric): ...  # undocumented
class _PyWalleniusNCHypergeometric(_PyNCHypergeometric): ...  # undocumented

# undocumented
class _PyStochasticLib3:
    def Random(self, /) -> float: ...
    def SetAccuracy(self, /, accur: op.CanFloat) -> None: ...
    def FishersNCHyp(self, /, n: op.CanInt, m: op.CanInt, N: op.CanInt, odds: op.CanFloat) -> _PyFishersNCHypergeometric: ...
    def WalleniusNCHyp(self, /, n: op.CanInt, m: op.CanInt, N: op.CanInt, odds: op.CanFloat) -> _PyWalleniusNCHypergeometric: ...
    #
    def rvs_fisher(
        self,
        /,
        n: op.CanInt,
        m: op.CanInt,
        N: op.CanInt,
        odds: op.CanFloat,
        size: op.CanInt,
        random_state: onp.random.RNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
    def rvs_wallenius(
        self,
        /,
        n: op.CanInt,
        m: op.CanInt,
        N: op.CanInt,
        odds: op.CanFloat,
        size: op.CanInt,
        random_state: onp.random.RNG | None = None,
    ) -> onp.Array1D[np.float64]: ...
