from collections.abc import Callable, Mapping, Sequence
from typing import Any, ClassVar, Final, Generic, Literal, TypeAlias, final
from typing_extensions import TypeVar, override

import optype as op

from scipy._lib._uarray._backend import _Backend

_T = TypeVar("_T", default=Any)
_T_co = TypeVar("_T_co", covariant=True, default=Any)

_ToBackend: TypeAlias = type[_Backend[_T]] | Literal["scipy"]

###

@final
class _ScipyBackend(_Backend[_T_co], Generic[_T_co]):
    __ua_domain__: ClassVar = "numpy.scipy.fft"

    @override
    @staticmethod
    def __ua_function__(method: Callable[..., _T_co], args: Sequence[object], kwargs: Mapping[str, object]) -> _T_co: ...

_named_backends: Final[dict[str, type[_Backend[Any]]]] = ...

def _backend_from_arg(backend: _ToBackend[_T]) -> type[_Backend[_T]]: ...  # undocumented
def set_global_backend(backend: _ToBackend, coerce: bool = False, only: bool = False, try_last: bool = False) -> None: ...
def register_backend(backend: _ToBackend) -> None: ...
def set_backend(backend: _ToBackend, coerce: bool = False, only: bool = False) -> op.CanWith[None, None]: ...
def skip_backend(backend: _ToBackend) -> op.CanWith[None, None]: ...
