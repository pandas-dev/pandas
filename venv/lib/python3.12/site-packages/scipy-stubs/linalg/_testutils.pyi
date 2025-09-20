from collections.abc import Callable, Iterable
from typing import Any, Final, Generic, TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp

_T = TypeVar("_T")
_ScalarT = TypeVar("_ScalarT", bound=np.generic)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

_Ignored: TypeAlias = object

class _FakeMatrix(Generic[_ScalarT, _ShapeT]):  # undocumented
    _data: onp.ArrayND[_ScalarT, _ShapeT]

    def __init__(self, /, data: onp.ArrayND[_ScalarT, _ShapeT]) -> None: ...
    __array_interface__: Final[Callable[[], dict[str, Any]]] = ...

class _FakeMatrix2(Generic[_ScalarT, _ShapeT]):  # undocumented
    _data: onp.ArrayND[_ScalarT, _ShapeT]
    def __init__(self, /, data: onp.ArrayND[_ScalarT, _ShapeT]) -> None: ...
    def __array__(self, /, dtype: np.dtype[Any] | None = None, copy: bool | None = None) -> onp.ArrayND[_ScalarT, _ShapeT]: ...

def _get_array(shape: _ShapeT, dtype: type[_ScalarT]) -> onp.ArrayND[_ScalarT, _ShapeT]: ...  # undocumented
def _id(x: _T) -> _T: ...  # undocumented

#
def assert_no_overwrite(
    call: Callable[..., _Ignored], shapes: Iterable[tuple[int, ...]], dtypes: Iterable[type[np.generic]] | None = None
) -> None: ...  # undocumented
