from collections.abc import Callable
from contextlib import _GeneratorContextManager
from typing import Final, ParamSpec, TypeVar

__all__ = ["assert_deallocated", "gc_state", "set_gc_state"]

###

_T = TypeVar("_T")
_Tss = ParamSpec("_Tss")

###

IS_PYPY: Final[bool] = ...

class ReferenceError(AssertionError): ...

def set_gc_state(state: bool) -> None: ...
def gc_state(state: bool) -> _GeneratorContextManager[None]: ...
def assert_deallocated(func: Callable[_Tss, _T], *args: _Tss.args, **kwargs: _Tss.kwargs) -> _GeneratorContextManager[_T]: ...
