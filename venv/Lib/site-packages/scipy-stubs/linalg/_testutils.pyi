from _typeshed import Unused
from collections.abc import Callable, Iterable
from typing import Any, Final

import numpy as np
import optype.numpy as onp

class _FakeMatrix[ScalarT: np.generic, ShapeT: tuple[int, ...]]:  # undocumented
    _data: onp.ArrayND[ScalarT, ShapeT]

    def __init__(self, /, data: onp.ArrayND[ScalarT, ShapeT]) -> None: ...
    __array_interface__: Final[Callable[[], dict[str, Any]]] = ...

    @property
    def shape(self) -> ShapeT: ...

class _FakeMatrix2[ScalarT: np.generic, ShapeT: tuple[int, ...]]:  # undocumented
    _data: onp.ArrayND[ScalarT, ShapeT]
    def __init__(self, /, data: onp.ArrayND[ScalarT, ShapeT]) -> None: ...
    def __array__(self, /, dtype: np.dtype[Any] | None = None, copy: bool | None = None) -> onp.ArrayND[ScalarT, ShapeT]: ...
    @property
    def shape(self) -> ShapeT: ...

def _get_array[ScalarT: np.generic, ShapeT: tuple[int, ...]](
    shape: ShapeT, dtype: type[ScalarT]
) -> onp.ArrayND[ScalarT, ShapeT]: ...  # undocumented
def _id[T](x: T) -> T: ...  # undocumented

#
def assert_no_overwrite(
    call: Callable[..., Unused], shapes: Iterable[tuple[int, ...]], dtypes: Iterable[type[np.generic]] | None = None
) -> None: ...  # undocumented
