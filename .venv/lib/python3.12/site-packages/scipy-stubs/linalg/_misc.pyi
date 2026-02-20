from typing import Any, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy.linalg import LinAlgError

__all__ = ["LinAlgError", "LinAlgWarning", "norm"]

_Inf: TypeAlias = float
_Order: TypeAlias = Literal["fro", "nuc", 0, 1, -1, 2, -2] | _Inf
_Axis: TypeAlias = op.CanIndex | tuple[op.CanIndex, op.CanIndex]
_SubScalar: TypeAlias = npc.inexact64 | npc.integer | np.bool_

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
_WorkaroundForPyright: TypeAlias = tuple[int] | tuple[Any, ...]

###

class LinAlgWarning(RuntimeWarning): ...

# NOTE: the mypy errors are false positives (join vs union)

@overload  # scalar, axis: None = ...
def norm(
    a: complex | _SubScalar,
    ord: _Order | None = None,
    axis: None = None,
    keepdims: op.CanBool = False,
    check_finite: onp.ToBool = True,
) -> np.float64: ...
@overload  # inexact32, axis: None = ...
def norm(
    a: npc.inexact32, ord: _Order | None = None, axis: None = None, keepdims: op.CanBool = False, check_finite: onp.ToBool = True
) -> np.float32: ...
@overload  # longdouble, axis: None = ...
def norm(
    a: npc.inexact80, ord: _Order | None = None, axis: None = None, keepdims: op.CanBool = False, check_finite: onp.ToBool = True
) -> np.longdouble: ...
@overload  # scalar array, axis: None = ..., keepdims: False = ...
def norm(
    a: onp.CanArrayND[_SubScalar] | onp.SequenceND[onp.CanArrayND[_SubScalar]] | onp.SequenceND[_SubScalar],
    ord: _Order | None = None,
    axis: None = None,
    keepdims: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> np.float64: ...
@overload  # float64-coercible array, keepdims: True (positional)
def norm(
    a: onp.ArrayND[_SubScalar, _ShapeT],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # float64-coercible array, keepdims: True (keyword)
def norm(
    a: onp.ArrayND[_SubScalar, _ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # float64-coercible array-like, keepdims: True (positional)
def norm(
    a: onp.SequenceND[onp.CanArrayND[_SubScalar]] | onp.SequenceND[complex | _SubScalar],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # float64-coercible array-like, keepdims: True (keyword)
def norm(
    a: onp.SequenceND[onp.CanArrayND[_SubScalar]] | onp.SequenceND[complex | _SubScalar],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # shaped inexact32 array, keepdims: True (positional)
def norm(
    a: onp.ArrayND[npc.inexact32, _ShapeT],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # shaped longdouble array, keepdims: True (positional)
def norm(
    a: onp.ArrayND[npc.inexact80, _ShapeT],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.longdouble, _ShapeT]: ...
@overload  # shaped inexact32 array, keepdims: True (keyword)
def norm(
    a: onp.ArrayND[npc.inexact32, _ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # shaped longdouble array, keepdims: True (keyword)
def norm(
    a: onp.ArrayND[npc.inexact80, _ShapeT],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.longdouble, _ShapeT]: ...
@overload  # scalar array-like, keepdims: True (positional)
def norm(
    a: onp.SequenceND[onp.CanArrayND[npc.inexact32]] | onp.SequenceND[npc.inexact32],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # scalar array-like, keepdims: True (positional)
def norm(
    a: onp.SequenceND[onp.CanArrayND[npc.inexact80]] | onp.SequenceND[npc.inexact80],
    ord: _Order | None,
    axis: _Axis | None,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.longdouble]: ...
@overload  # scalar array-like, keepdims: True (keyword)
def norm(
    a: onp.SequenceND[onp.CanArrayND[npc.inexact32]] | onp.SequenceND[npc.inexact32],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # scalar array-like, keepdims: True (keyword)
def norm(
    a: onp.SequenceND[onp.CanArrayND[npc.inexact80]] | onp.SequenceND[npc.inexact80],
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[np.longdouble]: ...
@overload  # array-like, axis: None = ..., keepdims: False = ...
def norm(
    a: onp.ToComplexND,
    ord: _Order | None = None,
    axis: None = None,
    keepdims: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> np.float64: ...
@overload  # array-like, keepdims: True (positional)
def norm(
    a: onp.ToComplexND, ord: _Order | None, axis: _Axis | None, keepdims: onp.ToTrue, check_finite: onp.ToBool = True
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...
@overload  # array-like, keepdims: True (keyword)
def norm(
    a: onp.ToComplexND,
    ord: _Order | None = None,
    axis: _Axis | None = None,
    *,
    keepdims: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...
@overload  # catch-all
def norm(
    a: onp.ToArrayND,
    ord: _Order | None = None,
    axis: _Axis | None = None,
    keepdims: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> npc.floating | onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...

#
def _datacopied(arr: onp.ArrayND[Any], original: onp.CanArrayND[Any]) -> bool: ...  # undocumented
