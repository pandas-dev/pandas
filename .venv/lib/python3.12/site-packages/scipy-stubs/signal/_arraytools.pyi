from typing import TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp

_SCT = TypeVar("_SCT", bound=np.generic)

###

def axis_slice(
    a: onp.ArrayND[_SCT],
    start: op.CanIndex | None = None,
    stop: op.CanIndex | None = None,
    step: op.CanIndex | None = None,
    axis: op.CanIndex = -1,
) -> onp.ArrayND[_SCT]: ...
def axis_reverse(a: onp.ArrayND[_SCT], axis: op.CanIndex = -1) -> onp.ArrayND[_SCT]: ...

#
def odd_ext(x: onp.ArrayND[_SCT], n: onp.ToInt, axis: op.CanIndex = -1) -> onp.ArrayND[_SCT]: ...
def even_ext(x: onp.ArrayND[_SCT], n: onp.ToInt, axis: op.CanIndex = -1) -> onp.ArrayND[_SCT]: ...
def const_ext(x: onp.ArrayND[_SCT], n: onp.ToInt, axis: op.CanIndex = -1) -> onp.ArrayND[_SCT]: ...
def zero_ext(x: onp.ArrayND[_SCT], n: onp.ToInt, axis: op.CanIndex = -1) -> onp.ArrayND[_SCT]: ...

#
@overload
def _validate_fs(fs: None, allow_none: onp.ToTrue = True) -> None: ...
@overload
def _validate_fs(fs: onp.ToFloat, allow_none: onp.ToBool = True) -> float: ...
