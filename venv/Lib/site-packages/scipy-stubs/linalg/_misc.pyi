from typing import Any, Literal, Never, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy.linalg import LinAlgError

__all__ = ["LinAlgError", "LinAlgWarning", "bandwidth", "norm"]

###

type _Inf = float
type _Order = Literal["fro", "nuc", 0, 1, -1, 2, -2] | _Inf
type _Axis = SupportsIndex | tuple[SupportsIndex, SupportsIndex]
type _SubScalar = npc.inexact64 | npc.integer | np.bool

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
type _WorkaroundForPyright = tuple[int] | tuple[Any, ...]

###

# On numpy<2.1, pyright reports 6 false positive incompatible overload errors here.
# pyright: reportOverlappingOverload=false

class LinAlgWarning(RuntimeWarning): ...

# NOTE: mypy reports two false positive `overload-overlap` error with numpy<2.5
# mypy: disable-error-code=overload-overlap

@overload  # scalar, axis: None = ...
def norm(
    a: complex | _SubScalar, ord: _Order | None = None, axis: None = None, keepdims: bool = False, check_finite: bool = True
) -> np.float64: ...
@overload  # inexact32, axis: None = ...
def norm(
    a: npc.inexact32, ord: _Order | None = None, axis: None = None, keepdims: bool = False, check_finite: bool = True
) -> np.float32: ...
@overload  # longdouble, axis: None = ...
def norm(
    a: npc.inexact80, ord: _Order | None = None, axis: None = None, keepdims: bool = False, check_finite: bool = True
) -> np.longdouble: ...
@overload  # scalar array, axis: None = ..., keepdims: False = ...
def norm(
    a: onp.ToArrayND[_SubScalar, _SubScalar],
    ord: _Order | None = None,
    axis: None = None,
    keepdims: Literal[False] = False,
    check_finite: bool = True,
) -> np.float64: ...
@overload  # float64-coercible array, keepdims: True (positional)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[_SubScalar, ShapeT], ord: _Order | None, axis: _Axis | None, keepdims: Literal[True], check_finite: bool = True
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # float64-coercible array, keepdims: True (keyword)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[_SubScalar, ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # float64-coercible array-like, keepdims: True (positional)
def norm(
    a: onp.ToArrayND[complex, _SubScalar],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # float64-coercible array-like, keepdims: True (keyword)
def norm(
    a: onp.ToArrayND[complex, _SubScalar],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # shaped inexact32 array, keepdims: True (positional)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.inexact32, ShapeT],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload  # shaped longdouble array, keepdims: True (positional)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.inexact80, ShapeT],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload  # shaped inexact32 array, keepdims: True (keyword)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.inexact32, ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload  # shaped longdouble array, keepdims: True (keyword)
def norm[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.inexact80, ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload  # scalar array-like, keepdims: True (positional)
def norm(
    a: onp.ToArrayND[npc.inexact32, npc.inexact32],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # scalar array-like, keepdims: True (positional)
def norm(
    a: onp.ToArrayND[npc.inexact80, npc.inexact80],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.longdouble]: ...
@overload  # scalar array-like, keepdims: True (keyword)
def norm(
    a: onp.ToArrayND[npc.inexact32, npc.inexact32],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # scalar array-like, keepdims: True (keyword)
def norm(
    a: onp.ToArrayND[npc.inexact80, npc.inexact80],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[np.longdouble]: ...
@overload  # array-like, axis: None = ..., keepdims: False = ...
def norm(
    a: onp.ToComplexND, ord: _Order | None = None, axis: None = None, keepdims: Literal[False] = False, check_finite: bool = True
) -> np.float64: ...
@overload  # array-like, keepdims: True (positional)
def norm(
    a: onp.ToComplexND, ord: _Order | None, axis: _Axis | None, keepdims: Literal[True], check_finite: bool = True
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...
@overload  # array-like, keepdims: True (keyword)
def norm(
    a: onp.ToComplexND,
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: Literal[True],
    check_finite: bool = True,
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...
@overload  # catch-all
def norm(
    a: onp.ToArrayND, ord: _Order | None = None, axis: _Axis | None = None, keepdims: bool = False, check_finite: bool = True
) -> npc.floating | onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...

#
def _datacopied(arr: onp.ArrayND[Any], original: onp.CanArrayND[Any]) -> bool: ...  # undocumented

#
@overload  # pyright workaround
def bandwidth(a: onp.ArrayND[npc.number, tuple[Never, Never, Never, Never]]) -> tuple[np.int64 | Any, np.int64 | Any]: ...
@overload
def bandwidth(a: onp.ToComplexStrict2D) -> tuple[np.int64, np.int64]: ...
@overload
def bandwidth(a: onp.ToComplexStrict3D) -> tuple[onp.Array1D[np.int64], onp.Array1D[np.int64]]: ...
@overload
def bandwidth(a: onp.ToComplexND) -> tuple[onp.ArrayND[np.int64] | Any, onp.ArrayND[np.int64] | Any]: ...
