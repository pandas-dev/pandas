import datetime as dt
from typing import Any, Final, Never, TypeVar, overload
from typing_extensions import CapsuleType

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_ArrayT = TypeVar("_ArrayT", bound=onp.ArrayND[Any, onp.AtLeast1D])

###

__pyx_capi__: Final[dict[str, CapsuleType]] = ...

#
@overload
def squeeze_element(arr: onp.ArrayND[Any, tuple[Never]]) -> Any: ...
@overload
def squeeze_element(arr: _ArrayT) -> _ArrayT: ...

# TODO: use generic item type of `np.generic` when we require `numpy>=2.2`
@overload
def squeeze_element(arr: onp.Array0D[np.bool_]) -> bool: ...
@overload
def squeeze_element(arr: onp.Array0D[npc.integer]) -> int: ...
@overload
def squeeze_element(arr: onp.Array0D[np.float16 | np.float32 | np.float64]) -> float: ...  # type: ignore[overload-overlap]
@overload
def squeeze_element(arr: onp.Array0D[npc.floating80]) -> np.longdouble: ...
@overload
def squeeze_element(arr: onp.Array0D[np.complex64 | np.complex128]) -> complex: ...  # type: ignore[overload-overlap]
@overload
def squeeze_element(arr: onp.Array0D[npc.complexfloating160]) -> np.clongdouble: ...
@overload
def squeeze_element(arr: onp.Array0D[np.bytes_]) -> bytes: ...
@overload
def squeeze_element(arr: onp.Array0D[np.str_]) -> str: ...
@overload
def squeeze_element(arr: onp.Array0D[np.datetime64]) -> dt.date | int | None: ...
@overload
def squeeze_element(arr: onp.Array0D[np.timedelta64]) -> dt.timedelta | int | None: ...
@overload
def squeeze_element(arr: onp.Array0D[np.object_]) -> Any: ...

#
def chars_to_strings(in_arr: onp.ArrayND[np.str_]) -> onp.ArrayND[np.str_]: ...
