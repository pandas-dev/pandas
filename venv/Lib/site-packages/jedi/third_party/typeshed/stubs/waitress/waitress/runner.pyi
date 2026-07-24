from _typeshed import Unused
from collections.abc import Callable, Sequence
from io import TextIOWrapper
from typing import Final

HELP: Final[str]

def show_help(stream: TextIOWrapper, name: str, error: str | None = None) -> None: ...
def show_exception(stream: TextIOWrapper) -> None: ...
def run(argv: Sequence[str] = ..., _serve: Callable[..., Unused] = ...) -> None: ...
