from collections.abc import Callable, Mapping, Sequence
from typing import Any, ClassVar, Final, Generic, Literal, final, override
from typing_extensions import TypeVar

import optype as op

from scipy._lib._uarray._backend import _Backend

_T_co = TypeVar("_T_co", covariant=True, default=Any)

type _ToBackend[T] = type[_Backend[T]] | Literal["scipy"]

###

@final
class _ScipyBackend(_Backend[_T_co], Generic[_T_co]):
    __ua_domain__: ClassVar = "numpy.scipy.fft"

    @override
    @staticmethod
    def __ua_function__(method: Callable[..., _T_co], args: Sequence[object], kwargs: Mapping[str, object]) -> _T_co: ...

_named_backends: Final[dict[str, type[_Backend[Any]]]] = ...

def _backend_from_arg[T](backend: _ToBackend[T]) -> type[_Backend[T]]: ...
def set_global_backend(backend: _ToBackend[Any], coerce: bool = False, only: bool = False, try_last: bool = False) -> None: ...
def register_backend(backend: _ToBackend[Any]) -> None: ...
def set_backend(backend: _ToBackend[Any], coerce: bool = False, only: bool = False) -> op.CanWith[None, None]: ...
def skip_backend(backend: _ToBackend[Any]) -> op.CanWith[None, None]: ...
