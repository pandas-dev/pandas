from _typeshed import Incomplete
from typing import Literal

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _Mode = Literal[0, 1, 2]

###

# defined in scipy/signal/_correlate_nd.cc
def _correlateND[NumberT: npc.number](
    x: onp.ArrayND[NumberT], y: onp.ArrayND[NumberT], out: onp.ArrayND[NumberT], mode: _Mode = 2
) -> onp.ArrayND[NumberT]: ...

# defined in scipy/signal/_sigtoolsmodule.cc
def _convolve2d[NumberT: npc.number](
    in1: onp.ArrayND[NumberT],
    in2: onp.ArrayND[NumberT],
    flip: int = 1,
    mode: _Mode = 2,
    boundary: int = 0,
    fillvalue: Incomplete | None = None,
) -> Incomplete: ...

# defined in scipy/signal/_lfilter.cc
def _linear_filter[NumberT: npc.number](
    b: onp.ArrayND[NumberT],
    a: onp.ArrayND[NumberT],
    X: onp.ArrayND[NumberT],
    axis: int = -1,
    Vi: onp.ArrayND[NumberT] | None = None,
) -> onp.ArrayND[NumberT]: ...

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
def _medfilt2d[ST: np.uint8 | np.float32 | np.float64](image: onp.Array2D[ST], size: tuple[int, int]) -> onp.Array2D[ST]: ...
