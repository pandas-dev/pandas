from _typeshed import Incomplete
from typing import Literal, TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_NumberT = TypeVar("_NumberT", bound=npc.number)
_ImageScalarT = TypeVar("_ImageScalarT", bound=np.uint8 | np.float32 | np.float64)
_Mode: TypeAlias = Literal[0, 1, 2]

###

# defined in scipy/signal/_correlate_nd.cc
def _correlateND(
    x: onp.ArrayND[_NumberT], y: onp.ArrayND[_NumberT], out: onp.ArrayND[_NumberT], mode: _Mode = 2
) -> onp.ArrayND[_NumberT]: ...

# defined in scipy/signal/_sigtoolsmodule.cc
def _convolve2d(
    in1: onp.ArrayND[_NumberT],
    in2: onp.ArrayND[_NumberT],
    flip: int = 1,
    mode: _Mode = 2,
    boundary: int = 0,
    fillvalue: Incomplete | None = None,
) -> Incomplete: ...

# defined in scipy/signal/_lfilter.cc
def _linear_filter(
    b: onp.ArrayND[_NumberT],
    a: onp.ArrayND[_NumberT],
    X: onp.ArrayND[_NumberT],
    axis: int = -1,
    Vi: onp.ArrayND[_NumberT] | None = None,
) -> onp.ArrayND[_NumberT]: ...

# defined in scipy/signal/_sigtoolsmodule.cc
def _remez(
    numtaps: int,
    bands: onp.ArrayND[np.float64],
    des: onp.ArrayND[np.float64],
    weight: onp.ArrayND[np.float64],
    type: Literal[1, 2, 3] = 1,
    fs: float = 1.0,
    maxiter: int = 25,
    grid_density: int = 16,
) -> onp.ArrayND[np.float64]: ...

#
# defined in scipy/signal/_sigtoolsmodule.cc
def _medfilt2d(image: onp.Array2D[_ImageScalarT], size: tuple[int, int]) -> onp.Array2D[_ImageScalarT]: ...
