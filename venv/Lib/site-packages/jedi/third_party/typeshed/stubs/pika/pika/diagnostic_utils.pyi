from collections.abc import Callable
from logging import Logger
from typing import Any, TypeVar

_F = TypeVar("_F", bound=Callable[..., Any])

def create_log_exception_decorator(logger: Logger) -> Callable[[_F], _F]: ...
