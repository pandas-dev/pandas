# defined in scipy/spatial/src/distance_pybind.cpp

from typing import Any

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

def cdist_braycurtis(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_canberra(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_chebyshev(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_cityblock(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_dice(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_euclidean(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_hamming(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_jaccard(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_kulczynski1(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_minkowski(
    x: onp.ToComplex2D,
    y: onp.ToComplex2D,
    w: onp.ToFloat1D | None = None,
    out: onp.ArrayND[npc.inexact] | None = None,
    p: float = ...,
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_rogerstanimoto(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_russellrao(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_sokalmichener(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_sokalsneath(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_sqeuclidean(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def cdist_yule(
    x: onp.ToComplex2D, y: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...

#
def pdist_braycurtis(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_canberra(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_chebyshev(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_cityblock(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_dice(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_euclidean(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_hamming(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_jaccard(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_kulczynski1(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_minkowski(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None, p: float = ...
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_rogerstanimoto(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_russellrao(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_sokalmichener(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_sokalsneath(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_sqeuclidean(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
def pdist_yule(
    x: onp.ToComplex2D, w: onp.ToFloat1D | None = None, out: onp.ArrayND[npc.inexact] | None = None
) -> onp.ArrayND[np.float64 | Any]: ...
