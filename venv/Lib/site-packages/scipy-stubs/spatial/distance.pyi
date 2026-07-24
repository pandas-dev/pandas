from collections.abc import Callable, Sequence
from typing import Any, Literal, Never, overload

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

type _MetricName = Literal[
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

###

type _MetricFunc = Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat | None]
type _Metric = _MetricName | _MetricFunc  # noqa: PYI047

type _Force = Literal["NO", "No", "no", "TOMATRIX", "ToMatrix", "tomatrix", "TOVECTOR", "ToVector", "tovector"]

# workaround for mypy & pyright's failure to conform to the overload typing specification
type _JustAnyShape = tuple[Never, Never, Never, Never]

type _ToFloatStrictND = onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]

###

# NOTE: on numpy<2.1, both mypy and pyright reports false positive overload-overlap errors for overload 2 of `jensenshannon`
# pyright: reportOverlappingOverload=false
# mypy: disable-error-code=overload-overlap

# TODO(@jorenham): metric-specific overloads
# https://github.com/scipy/scipy-stubs/issues/404
@overload
def cdist(
    XA: onp.ToFloat2D,
    XB: onp.ToFloat2D,
    metric: _MetricName = "euclidean",
    *,
    out: onp.Array2D[np.float64] | None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> onp.Array2D[np.float64]: ...
@overload
def cdist(
    XA: onp.ToFloat2D, XB: onp.ToFloat2D, metric: _MetricFunc, *, out: onp.Array2D[np.float64] | None = None, **kwds: object
) -> onp.Array2D[np.float64]: ...

# TODO(@jorenham): metric-specific overloads
# https://github.com/scipy/scipy-stubs/issues/404
@overload
def pdist(
    X: onp.ToFloat2D,
    metric: _MetricName = "euclidean",
    *,
    out: onp.Array1D[np.float64] | None = None,
    p: float = 2,
    w: onp.ToFloat1D | None = None,
    V: onp.ToFloat2D | None = None,
    VI: onp.ToFloat2D | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def pdist(
    X: onp.ToFloat2D, metric: _MetricFunc, *, out: onp.Array1D[np.float64] | None = None, **kwargs: object
) -> onp.Array1D[np.float64]: ...

#
@overload  # ?d T@number
def squareform[NumberT: npc.number](
    X: onp.ArrayND[NumberT, _JustAnyShape], force: _Force = "no", checks: bool = True
) -> onp.ArrayND[NumberT]: ...
@overload  # 1d +int
def squareform(X: Sequence[int], force: _Force = "no", checks: bool = True) -> onp.Array2D[np.int_]: ...
@overload  # 1d ~float
def squareform(X: list[float], force: _Force = "no", checks: bool = True) -> onp.Array2D[np.float64]: ...
@overload  # 1d ~complex
def squareform(X: list[complex], force: _Force = "no", checks: bool = True) -> onp.Array2D[np.complex128]: ...
@overload  # 1d T@number
def squareform[NumberT: npc.number](
    X: onp.ToArrayStrict1D[NumberT, NumberT], force: _Force = "no", checks: bool = True
) -> onp.Array2D[NumberT]: ...
@overload  # 2d +int
def squareform(X: Sequence[Sequence[int]], force: _Force = "no", checks: bool = True) -> onp.Array1D[np.int_]: ...
@overload  # 2d ~float
def squareform(X: Sequence[list[float]], force: _Force = "no", checks: bool = True) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def squareform(X: Sequence[list[complex]], force: _Force = "no", checks: bool = True) -> onp.Array1D[np.complex128]: ...
@overload  # 2d T@number
def squareform[NumberT: npc.number](
    X: onp.ToArrayStrict2D[NumberT, NumberT], force: _Force = "no", checks: bool = True
) -> onp.Array1D[NumberT]: ...
@overload  # fallback
def squareform(X: onp.ToComplexND, force: _Force = "no", checks: bool = True) -> onp.ArrayND: ...

#
def correlation(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None, centered: bool = True) -> np.float64: ...

#
def cosine(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...

#
@overload
def mahalanobis(u: onp.ToFloat1D, v: onp.ToFloat1D, VI: onp.ToFloat2D) -> np.float64: ...
@overload
def mahalanobis(u: onp.ToJustComplex1D, v: onp.ToComplex1D, VI: onp.ToComplex2D) -> np.complex128: ...
@overload
def mahalanobis(u: onp.ToComplex1D, v: onp.ToJustComplex1D, VI: onp.ToComplex2D) -> np.complex128: ...
@overload
def mahalanobis(u: onp.ToComplex1D, v: onp.ToComplex1D, VI: onp.ToJustComplex2D) -> np.complex128: ...
@overload
def mahalanobis(u: onp.ToComplex1D, v: onp.ToComplex1D, VI: onp.ToComplex2D) -> np.float64 | np.complex128: ...

#
@overload
def sokalsneath(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
@overload
def sokalsneath(u: onp.ToJustComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.complex128: ...
@overload
def sokalsneath(u: onp.ToComplex1D, v: onp.ToJustComplex1D, w: onp.ToFloat1D | None = None) -> np.complex128: ...
@overload
def sokalsneath(u: onp.ToComplex1D, v: onp.ToComplex1D, w: onp.ToFloat1D | None = None) -> np.float64 | np.complex128: ...

#
@overload  # ?d, keepdims=False
def jensenshannon(
    p: _ToFloatStrictND, q: _ToFloatStrictND, base: float | None = None, *, axis: int = 0, keepdims: Literal[False] = False
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d, keepdims=True
def jensenshannon[ShapeT: tuple[int, ...]](
    p: onp.ArrayND[npc.floating | npc.integer, ShapeT],
    q: onp.ArrayND[npc.floating | npc.integer, ShapeT],
    base: float | None = None,
    *,
    axis: int = 0,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # 1d, keepdims=False
def jensenshannon(
    p: onp.ToFloatStrict1D, q: onp.ToFloatStrict1D, base: float | None = None, *, axis: int = 0, keepdims: Literal[False] = False
) -> np.float64: ...
@overload  # 1d, keepdims=True
def jensenshannon(
    p: onp.ToFloatStrict1D, q: onp.ToFloatStrict1D, base: float | None = None, *, axis: int = 0, keepdims: Literal[True]
) -> onp.Array1D[np.float64]: ...
@overload  # 2d, keepdims=False
def jensenshannon(
    p: onp.ToFloatStrict2D, q: onp.ToFloatStrict2D, base: float | None = None, *, axis: int = 0, keepdims: Literal[False] = False
) -> onp.Array1D[np.float64]: ...
@overload  # 2d, keepdims=True
def jensenshannon(
    p: onp.ToFloatStrict2D, q: onp.ToFloatStrict2D, base: float | None = None, *, axis: int = 0, keepdims: Literal[True]
) -> onp.Array2D[np.float64]: ...
@overload  # nd, keepdims=False
def jensenshannon(
    p: onp.ToFloatND, q: onp.ToFloatND, base: float | None = None, *, axis: int = 0, keepdims: Literal[False] = False
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # nd, keepdims=True
def jensenshannon(
    p: onp.ToFloatND, q: onp.ToFloatND, base: float | None = None, *, axis: int = 0, keepdims: Literal[True]
) -> onp.ArrayND[np.float64]: ...

# NOTE: These technically also accept complex inputs, but will unsafely downcast to float and emit a warning, so we disallow it.
def braycurtis(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def canberra(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def chebyshev(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def cityblock(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def dice(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> float: ...
def directed_hausdorff(
    u: onp.ToFloat2D, v: onp.ToFloat2D, rng: onp.random.ToRNG | None = 0, *, seed: onp.random.ToRNG | None = None
) -> tuple[float, int, int]: ...
def hamming(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def jaccard(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
def rogerstanimoto(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> float: ...
def russellrao(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> float: ...
def seuclidean(u: onp.ToFloat1D, v: onp.ToFloat1D, V: onp.ToFloat1D) -> float: ...
def yule(u: onp.ToFloat1D, v: onp.ToFloat1D, w: onp.ToFloat1D | None = None) -> float: ...

type _ToFloatMax2D = onp.ToFloatStrict2D | onp.ToFloatStrict1D
type _ToFloatMax3D = onp.ToFloatStrict3D | _ToFloatMax2D
type _ToComplexMax2D = onp.ToComplexStrict2D | onp.ToComplexStrict1D
type _ToJustComplexMax2D = onp.ToJustComplexStrict2D | onp.ToJustComplexStrict1D

#
@overload  # 1d, 1d
def minkowski(
    u: onp.ToFloatStrict1D, v: onp.ToFloatStrict1D, p: float = 2, w: onp.ToFloatStrict1D | None = None
) -> np.float64: ...
@overload  # 1d, 2d
def minkowski(
    u: onp.ToFloatStrict1D, v: onp.ToFloatStrict2D, p: float = 2, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # 2d, {1,2}d
def minkowski(
    u: onp.ToFloatStrict2D, v: _ToFloatMax2D, p: float = 2, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # {1,2}d, 3d
def minkowski(
    u: _ToFloatMax2D, v: onp.ToFloatStrict3D, p: float = 2, w: onp.ToFloat1D | None = None
) -> onp.Array2D[np.float64]: ...
@overload  # 3d, {1,2,3}d
def minkowski(
    u: onp.ToFloatStrict3D, v: _ToFloatMax3D, p: float = 2, w: onp.ToFloat1D | None = None
) -> onp.Array2D[np.float64]: ...
@overload  # Nd, Nd
def minkowski(
    u: onp.ToFloatND, v: onp.ToFloatND, p: float = 2, w: onp.ToFloat1D | None = None
) -> onp.ArrayND[np.float64] | Any: ...

#
@overload  # 1d, 1d
def euclidean(u: onp.ToFloatStrict1D, v: onp.ToFloatStrict1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
@overload  # 1d, 2d
def euclidean(u: onp.ToFloatStrict1D, v: onp.ToFloatStrict2D, w: onp.ToFloat1D | None = None) -> onp.Array1D[np.float64]: ...
@overload  # 2d, {1,2}d
def euclidean(u: onp.ToFloatStrict2D, v: _ToFloatMax2D, w: onp.ToFloat1D | None = None) -> onp.Array1D[np.float64]: ...
@overload  # {1,2}d, 3d
def euclidean(u: _ToFloatMax2D, v: onp.ToFloatStrict3D, w: onp.ToFloat1D | None = None) -> onp.Array2D[np.float64]: ...
@overload  # 3d, {1,2,3}d
def euclidean(u: onp.ToFloatStrict3D, v: _ToFloatMax3D, w: onp.ToFloat1D | None = None) -> onp.Array2D[np.float64]: ...
@overload  # Nd, Nd
def euclidean(u: onp.ToFloatND, v: onp.ToFloatND, w: onp.ToFloat1D | None = None) -> onp.ArrayND[np.float64] | Any: ...

#
@overload  # 1d float, 1d float
def sqeuclidean(u: onp.ToFloatStrict1D, v: onp.ToFloatStrict1D, w: onp.ToFloat1D | None = None) -> np.float64: ...
@overload  # 1d complex, 1d +complex
def sqeuclidean(u: onp.ToJustComplexStrict1D, v: onp.ToComplexStrict1D, w: onp.ToFloat1D | None = None) -> np.complex128: ...
@overload  # 1d +complex, 1d complex
def sqeuclidean(u: onp.ToComplexStrict1D, v: onp.ToJustComplexStrict1D, w: onp.ToFloat1D | None = None) -> np.complex128: ...
@overload  # 1d +complex, 1d +complex
def sqeuclidean(
    u: onp.ToComplexStrict1D, v: onp.ToComplexStrict1D, w: onp.ToFloat1D | None = None
) -> np.float64 | np.complex128: ...
@overload  # 1d float, 2d float
def sqeuclidean(u: onp.ToFloatStrict1D, v: onp.ToFloatStrict2D, w: onp.ToFloat1D | None = None) -> onp.Array1D[np.float64]: ...
@overload  # 1d complex, 2d +complex
def sqeuclidean(
    u: onp.ToJustComplexStrict1D, v: onp.ToComplexStrict2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # 1d +complex, 2d complex
def sqeuclidean(
    u: onp.ToComplexStrict1D, v: onp.ToJustComplexStrict2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # 1d +complex, 2d +complex
def sqeuclidean(
    u: onp.ToComplexStrict1D, v: onp.ToComplexStrict2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.float64 | np.complex128]: ...
@overload  # 2d float, {1,2}d float
def sqeuclidean(u: onp.ToFloatStrict2D, v: _ToFloatMax2D, w: onp.ToFloat1D | None = None) -> onp.Array1D[np.float64]: ...
@overload  # 2d complex, {1,2}d +complex
def sqeuclidean(
    u: onp.ToJustComplexStrict2D, v: _ToComplexMax2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex, {1,2}d complex
def sqeuclidean(
    u: onp.ToComplexStrict2D, v: _ToJustComplexMax2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex, {1,2}d +complex
def sqeuclidean(
    u: onp.ToComplexStrict2D, v: _ToComplexMax2D, w: onp.ToFloat1D | None = None
) -> onp.Array1D[np.float64 | np.complex128]: ...
@overload  # Nd float, Nd float
def sqeuclidean(u: onp.ToFloatND, v: onp.ToFloatND, w: onp.ToFloat1D | None = None) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd complex, Nd +complex
def sqeuclidean(
    u: onp.ToJustComplexND, v: onp.ToComplexND, w: onp.ToFloat1D | None = None
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # Nd +complex, Nd complex
def sqeuclidean(
    u: onp.ToComplexND, v: onp.ToJustComplexND, w: onp.ToFloat1D | None = None
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # Nd +complex, Nd +complex
def sqeuclidean(
    u: onp.ToComplexND, v: onp.ToComplexND, w: onp.ToFloat1D | None = None
) -> onp.ArrayND[np.float64 | np.complex128] | Any: ...

#
def num_obs_dm(d: onp.ToArray2D) -> int: ...
def num_obs_y(Y: onp.ToArray1D) -> int: ...
def is_valid_dm(D: onp.ToArray2D, tol: float = 0.0, throw: bool = False, name: str = "D", warning: bool = False) -> bool: ...
def is_valid_y(y: onp.ToArray1D, warning: bool = False, throw: bool = False, name: str | None = None) -> bool: ...
