from __future__ import annotations

from .application import Application
from .config import Config, ConfigINI
from .switches import (
    CSV,
    CountOf,
    ExistingDirectory,
    ExistingFile,
    Flag,
    NonexistentPath,
    Predicate,
    Range,
    Set,
    SwitchAttr,
    SwitchError,
    autoswitch,
    positional,
    switch,
)

__all__ = (
    "CSV",
    "Application",
    "Config",
    "ConfigINI",
    "CountOf",
    "ExistingDirectory",
    "ExistingFile",
    "Flag",
    "NonexistentPath",
    "Predicate",
    "Range",
    "Set",
    "SwitchAttr",
    "SwitchError",
    "autoswitch",
    "positional",
    "switch",
)
