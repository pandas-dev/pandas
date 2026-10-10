from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # import modules that have public classes/functions
    from pandas.io import (
        formats,
        json,
        stata,
    )

    # mark only those modules as public
    __all__ = ["formats", "json", "stata"]


def __getattr__(name: str) -> object:
    if name in {"json", "stata"}:
        return import_module(f"pandas.io.{name}")
    raise AttributeError(f"module 'pandas.io' has no attribute '{name}'")


def __dir__() -> list[str]:
    return [*list(globals().keys()), "json", "stata"]
