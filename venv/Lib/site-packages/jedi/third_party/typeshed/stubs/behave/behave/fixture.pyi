from _typeshed import Incomplete
from collections.abc import Callable
from typing import Any, Concatenate, ParamSpec, TypeVar

from behave.runner import Context

_T = TypeVar("_T")
_F = TypeVar("_F", bound=Callable[..., Any])
_P = ParamSpec("_P")

def use_fixture(
    fixture_func: Callable[Concatenate[Context, _P], _T], context: Context, *fixture_args: _P.args, **fixture_kwargs: _P.kwargs
) -> _T: ...
def fixture(func: _F | None = None, name: str | None = None, pattern: str | None = None) -> _F: ...
def __getattr__(name: str) -> Incomplete: ...
