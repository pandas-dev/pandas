from __future__ import annotations

__lazy_modules__ = {
    f"{__spec__.parent}.application",
    f"{__spec__.parent}.config",
    f"{__spec__.parent}.switches",
}

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
