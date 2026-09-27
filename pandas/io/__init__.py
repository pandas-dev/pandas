# ruff: noqa: TC004
from typing import TYPE_CHECKING

__lazy_modules__ = (
    "pandas.io.json",
    "pandas.io.stata",
)

from pandas.io import (
    json,
    stata,
)

if TYPE_CHECKING:
    # import modules that have public classes/functions
    from pandas.io import (
        formats,
        json,
        stata,
    )

    # mark only those modules as public
    __all__ = ["formats", "json", "stata"]
