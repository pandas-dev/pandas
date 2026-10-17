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
    # submodules that PEP 810 lazy imports keep out of a plain `import pandas`
    # on Python 3.15+, so pd.io.<name> would otherwise disappear: the backends
    # pandas.io.api defers, and pandas.io.common (deferred by core.frame and
    # io.formats.format)
    api = import_module("pandas.io.api")
    return {"common", *(name.rsplit(".", 1)[1] for name in api.__lazy_modules__)}


def __getattr__(name: str) -> object:
    if name in _lazy_submodules():
        return import_module(f"pandas.io.{name}")
    raise AttributeError(f"module 'pandas.io' has no attribute '{name}'")


def __dir__() -> list[str]:
    return sorted({*globals(), *_lazy_submodules()})
