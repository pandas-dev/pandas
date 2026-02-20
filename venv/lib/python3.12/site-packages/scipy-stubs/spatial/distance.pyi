# NOTE: Scipy already has a `distance.pyi` stub, but it has several errors, which are fixed here

from collections.abc import Callable, Sequence as Seq
from typing import Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "braycurtis",
    "canberra",
    "cdist",
    "chebyshev",
    "cityblock",
    "correlation",
    "cosine",
    "dice",
    "directed_hausdorff",
    "euclidean",
    "hamming",
    "is_valid_dm",
    "is_valid_y",
    "jaccard",
    "jensenshannon",
    "mahalanobis",
    "minkowski",
    "num_obs_dm",
    "num_obs_y",
    "pdist",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalsneath",
    "sqeuclidean",
    "squareform",
    "yule",
]

_MetricName: TypeAlias = Literal[
    "braycurtis",
    "canberra",
    "chebychev",
    "chebyshev",
    "cheby",
    "cheb",
    "ch",
    "cityblock",
    "cblock",
    "cb",
    "c",
    "correlation",
    "co",
    "cosine",
    "cos",
    "dice",
    "euclidean",
    "euclid",
    "eu",
    "e",
    "hamming",
    "hamm",
    "ha",
    "h",
    "minkowski",
    "mi",
    "m",
    "pnorm",
    "jaccard",
    "jacc",
    "ja",
    "j",
    "jensenshannon",
    "js",
    "mahalanobis",
    "mahal",
    "mah",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "se",
    "s",
    "sokalsneath",
    "sqeuclidean",
    "sqe",
    "sqeuclid",
    "yule",
]
_MetricFunc: TypeAlias = Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat | None]
_Metric: TypeAlias = _MetricName | _MetricFunc  # noqa: PYI047

_Force: TypeAlias = Literal["NO", "No", "no", "TOMATRIX", "ToMatrix", "tomatrix", "TOVECTOR", "ToVector", "tovector"]

_Seq2D: TypeAlias = Seq[Seq[_T]]
_FloatingND: TypeAlias = onp.ArrayND[npc.floating]
_InexactND: TypeAlias = onp.ArrayND[npc.inexact]

_T = TypeVar("_T")
_NumberT = TypeVar("_NumberT", bound=npc.number)
_ArrayT = TypeVar("_ArrayT", bound=onp.ArrayND[npc.number])

###

#
@overload
def cdist(
    XA: onp.ToFloat2D,
    XB: onp.ToFloat2D,
    metric: _MetricName = "euclidean",
    *,
    out: None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _FloatingND: ...
@overload
def cdist(
    XA: onp.ToComplex2D,
    XB: onp.ToComplex2D,
    metric: _MetricName = "euclidean",
    *,
    out: None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _InexactND: ...
@overload
def cdist(
    XA: onp.ToComplex2D,
    XB: onp.ToComplex2D,
    metric: _MetricName = "euclidean",
    *,
    out: _ArrayT,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _ArrayT: ...
@overload
def cdist(XA: onp.ToFloat2D, XB: onp.ToFloat2D, metric: _MetricFunc, *, out: None = None, **kwds: object) -> _FloatingND: ...
@overload
def cdist(XA: onp.ToComplex2D, XB: onp.ToComplex2D, metric: _MetricFunc, *, out: None = None, **kwds: object) -> _InexactND: ...
@overload
def cdist(XA: onp.ToComplex2D, XB: onp.ToComplex2D, metric: _MetricFunc, *, out: _ArrayT, **kwds: object) -> _ArrayT: ...

#
@overload
def pdist(
    X: onp.ToFloat2D,
    metric: _MetricName = "euclidean",
    *,
    out: None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _FloatingND: ...
@overload
def pdist(
    X: onp.ToComplex2D,
    metric: _MetricName = "euclidean",
    *,
    out: None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _InexactND: ...
@overload
def pdist(
    X: onp.ToComplex2D,
    metric: _MetricName = "euclidean",
    *,
    out: _ArrayT,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> _ArrayT: ...
@overload
def pdist(X: onp.ToFloat2D, metric: _MetricFunc, *, out: None = None, **kwargs: object) -> _FloatingND: ...
@overload
def pdist(X: onp.ToComplex2D, metric: _MetricFunc, *, out: None = None, **kwargs: object) -> _InexactND: ...
@overload
def pdist(X: onp.ToComplex2D, metric: _MetricFunc, *, out: _ArrayT, **kwargs: object) -> _ArrayT: ...

#
@overload  # 1-d int
def squareform(X: onp.ToJustIntStrict1D, force: _Force = "no", checks: bool = True) -> onp.Array2D[np.int_]: ...
@overload  # 1-d float
def squareform(X: onp.ToJustFloatStrict1D, force: _Force = "no", checks: bool = True) -> onp.Array2D[np.float64]: ...
@overload  # 1-d complex
def squareform(X: onp.ToJustComplexStrict1D, force: _Force = "no", checks: bool = True) -> onp.Array2D[np.complex128]: ...
@overload  # 1-d array-like
def squareform(
    X: Seq[_NumberT] | onp.CanArray1D[_NumberT], force: _Force = "no", checks: bool = True
) -> onp.Array2D[_NumberT]: ...
@overload  # 2-d int
def squareform(X: onp.ToJustIntStrict2D, force: _Force = "no", checks: bool = True) -> onp.Array1D[np.int_]: ...
@overload  # 2-d float
def squareform(X: onp.ToJustFloatStrict2D, force: _Force = "no", checks: bool = True) -> onp.Array1D[np.float64]: ...
@overload  # 2-d complex
def squareform(X: onp.ToJustComplexStrict2D, force: _Force = "no", checks: bool = True) -> onp.Array1D[np.complex128]: ...
@overload  # 2-d array-like
def squareform(
    X: _Seq2D[_NumberT] | Seq[onp.Array1D[_NumberT]] | onp.CanArray2D[_NumberT], force: _Force = "no", checks: bool = True
) -> onp.Array1D[_NumberT]: ...
@overload  # ?-d array-like
def squareform(
    X: Seq[onp.CanArrayND[_NumberT]] | onp.CanArrayND[_NumberT], force: _Force = "no", checks: bool = True
) -> onp.Array1D[_NumberT] | onp.Array2D[_NumberT]: ...

#
def correlation(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None, centered: bool = True) -> np.float64: ...

#
def cosine(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...

#
@overload
def mahalanobis(u: onp.ToFloat1D, v: onp.ToFloat1D, VI: onp.ToFloat2D) -> np.float64: ...
@overload
def mahalanobis(u: onp.ToComplex1D, v: onp.ToComplex1D, VI: onp.ToComplex1D) -> np.float64 | np.complex128: ...

#
@overload
def sokalsneath(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
@overload
def sokalsneath(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.float64 | np.complex128: ...

#
@overload
def sqeuclidean(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> npc.floating: ...
@overload
def sqeuclidean(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> npc.inexact: ...

#
@overload
def jensenshannon(
    p: onp.ToFloatStrict1D,
    q: onp.ToFloatStrict1D,
    base: onp.ToFloat | None = None,
    *,
    axis: int = 0,
    keepdims: onp.ToFalse = False,
) -> np.float32 | np.float64: ...
@overload
def jensenshannon(
    p: onp.ToFloatStrict1D, q: onp.ToFloatStrict1D, base: onp.ToFloat | None = None, *, axis: int = 0, keepdims: onp.ToTrue
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload
def jensenshannon(
    p: onp.ToFloatND, q: onp.ToFloatND, base: onp.ToFloat | None = None, *, axis: int = 0, keepdims: bool = False
) -> np.float32 | np.float64 | onp.ArrayND[np.float32 | np.float64]: ...

# NOTE: The output of the following functions is always real, but complex input usually results in a `ComplexWarning` at runtime.
def braycurtis(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def canberra(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def chebyshev(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> npc.floating: ...
def cityblock(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> npc.floating: ...
def dice(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> float: ...
def directed_hausdorff(
    u: onp.ToComplex2D, v: onp.ToComplex2D, rng: onp.random.ToRNG | None = 0, *, seed: onp.random.ToRNG | None = None
) -> tuple[float, int, int]: ...
def hamming(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def euclidean(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> float: ...
def jaccard(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToBool1D | None = None) -> np.float64: ...
def minkowski(u: onp.ToComplex1D, v: onp.ToComplex1D, p: float = 2, w: onp.ToFloat1D | None = None) -> float: ...
def rogerstanimoto(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> float: ...
def russellrao(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> float: ...
def seuclidean(u: onp.ToComplex1D, v: onp.ToComplex1D, V: onp.ToFloat1D) -> float: ...
def yule(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> float: ...

#
def num_obs_dm(d: onp.ToArray2D) -> int: ...
def num_obs_y(Y: onp.ToArray1D) -> int: ...
def is_valid_dm(D: onp.ToArray2D, tol: float = 0.0, throw: bool = False, name: str = "D", warning: bool = False) -> bool: ...
def is_valid_y(y: onp.ToArray1D, warning: bool = False, throw: bool = False, name: str | None = None) -> bool: ...
