from collections.abc import Callable, Iterator
from contextlib import _GeneratorContextManager
from typing import Concatenate, Final, Generic, ParamSpec, overload
from typing_extensions import TypeVar, deprecated

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True, default=object)
_Tss = ParamSpec("_Tss")
_FT = TypeVar("_FT", bound=Callable[..., object])

__version__: Final[str] = ...

class ContextManager(_GeneratorContextManager[_T_co], Generic[_T_co]):
    def __init__(self, /, g: Callable[_Tss, Iterator[_T_co]], *a: _Tss.args, **k: _Tss.kwargs) -> None: ...

def contextmanager(g: Callable[_Tss, Iterator[_T]], *a: _Tss.args, **k: _Tss.kwargs) -> Callable[[], ContextManager[_T]]: ...

#
def decorate(func: _FT, caller: Callable[Concatenate[_FT, _Tss], _T]) -> Callable[_Tss, _T]: ...
@overload
def decorator(caller: Callable[Concatenate[_FT, _Tss], _T], _func: None = None) -> Callable[[_FT], Callable[_Tss, _T]]: ...
@overload
@deprecated("this is obsolete behavior; you should use decorate instead")
def decorator(caller: Callable[Concatenate[_FT, _Tss], _T], _func: _FT) -> Callable[_Tss, _T]: ...
