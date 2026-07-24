import logging
from collections.abc import Callable
from typing import Any, TypeVar
from typing_extensions import TypeAlias

import click

_AnyCallable: TypeAlias = Callable[..., Any]
_FC = TypeVar("_FC", bound=_AnyCallable | click.Command)

def simple_verbosity_option(logger: logging.Logger | str | None = None, *names: str, **kwargs: Any) -> Callable[[_FC], _FC]: ...
