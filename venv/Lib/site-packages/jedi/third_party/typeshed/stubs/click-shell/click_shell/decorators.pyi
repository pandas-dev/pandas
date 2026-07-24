from collections.abc import Callable
from typing import Any

from click.decorators import _AnyCallable

from .core import Shell

def shell(name: str | None = None, **attrs: Any) -> Callable[[_AnyCallable], Shell]: ...
