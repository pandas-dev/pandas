from _typeshed import Incomplete
from collections.abc import Callable
from typing import Concatenate, TypeVar

from behave.runner import Context

_F = TypeVar("_F", bound=Callable[Concatenate[Context, ...], None])

def given(step_text: str, **kwargs) -> Callable[[_F], _F]: ...
def when(step_text: str, **kwargs) -> Callable[[_F], _F]: ...
def then(step_text: str, **kwargs) -> Callable[[_F], _F]: ...
def step(step_text: str, **kwargs) -> Callable[[_F], _F]: ...

# Title-case aliases
Given = given
When = when
Then = then
Step = step

def __getattr__(name: str) -> Incomplete: ...
