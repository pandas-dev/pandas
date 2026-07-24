import dataclasses
import importlib.abc
from _typeshed import Unused
from collections.abc import Iterator
from importlib.machinery import ModuleSpec
from pathlib import Path
from types import ModuleType
from typing import Any, Final

from .globals import Indirect

__all__ = [
    "COMPAT_PACKAGE_NAME",
    "PACKAGE_NAME",
    "PluginSpec",
    "directories",
    "load_all_plugins",
    "load_plugins",
    "register_plugin_spec",
]

PACKAGE_NAME: Final = "yt_dlp_plugins"
COMPAT_PACKAGE_NAME: Final = "ytdlp_plugins"

@dataclasses.dataclass
class PluginSpec:
    module_name: str
    suffix: str
    destination: Indirect[dict[str, Any]]
    plugin_destination: Indirect[dict[str, Any]]

class PluginLoader(importlib.abc.Loader):
    def exec_module(self, module: ModuleType) -> None: ...

class PluginFinder(importlib.abc.MetaPathFinder):
    packages: set[str]
    def __init__(self, *packages: str) -> None: ...
    def search_locations(self, fullname: str) -> Iterator[Path]: ...
    def find_spec(self, fullname: str, path: Unused | None = None, target: Unused | None = None) -> ModuleSpec | None: ...
    def invalidate_caches(self) -> None: ...

def directories() -> list[str]: ...
def load_plugins(plugin_spec: PluginSpec) -> dict[str, Any]: ...
def load_all_plugins() -> None: ...
def register_plugin_spec(plugin_spec: PluginSpec) -> None: ...
