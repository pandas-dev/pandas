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


def _lazy_submodules() -> set[str]:
    # submodules pandas.io.api defers under PEP 810 and no longer imports
    # eagerly, so pd.io.<name> would otherwise disappear on Python 3.15+
    api = import_module("pandas.io.api")
    return {name.rsplit(".", 1)[1] for name in api.__lazy_modules__}


def __getattr__(name: str) -> object:
    if name in _lazy_submodules():
        return import_module(f"pandas.io.{name}")
    raise AttributeError(f"module 'pandas.io' has no attribute '{name}'")


def __dir__() -> list[str]:
    return [*globals(), *_lazy_submodules()]
