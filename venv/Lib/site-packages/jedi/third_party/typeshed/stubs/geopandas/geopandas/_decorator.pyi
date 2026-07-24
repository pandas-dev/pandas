from collections.abc import Callable
from typing import TypeVar, overload
from typing_extensions import TypeAlias

_AnyCallable: TypeAlias = Callable[..., object]
_Func = TypeVar("_Func", bound=_AnyCallable)

# We (ab)use this decorator to also copy the signature of the source function (overload 1)
# The advantages are:
# - avoid copying all parameters and types manually while conserving type safety
# - signature properly handeled in IDEs (at least with Pylance)
# - docstring from the original function properly displayed (at least with Pylance)
# Using the other overloads returns the signature of the decorated function instead
@overload
def doc(func: _Func, /, **params: object) -> Callable[[_AnyCallable], _Func]: ...
@overload
def doc(docstring: str, /, *docstrings: str | _AnyCallable, **params: object) -> Callable[[_Func], _Func]: ...
@overload
def doc(
    docstring1: str | _AnyCallable, docstring2: str | _AnyCallable, /, *docstrings: str | _AnyCallable, **params: object
) -> Callable[[_Func], _Func]: ...
