import configparser
import importlib.metadata
from collections.abc import Generator
from logging import Logger
from typing import Any, Final, NamedTuple

LOG: Logger
FLAKE8_GROUPS: Final[frozenset[str]]
BANNED_PLUGINS: dict[str, str]

class Plugin(NamedTuple):
    package: str
    version: str
    entry_point: importlib.metadata.EntryPoint

class LoadedPlugin(NamedTuple):
    plugin: Plugin
    obj: Any
    parameters: dict[str, bool]
    @property
    def entry_name(self) -> str: ...
    @property
    def display_name(self) -> str: ...

class Checkers(NamedTuple):
    tree: list[LoadedPlugin]
    logical_line: list[LoadedPlugin]
    physical_line: list[LoadedPlugin]

class Plugins(NamedTuple):
    checkers: Checkers
    reporters: dict[str, LoadedPlugin]
    disabled: list[LoadedPlugin]
    def all_plugins(self) -> Generator[LoadedPlugin]: ...
    def versions_str(self) -> str: ...

class PluginOptions(NamedTuple):
    local_plugin_paths: tuple[str, ...]
    enable_extensions: frozenset[str]
    require_plugins: frozenset[str]
    @classmethod
    def blank(cls) -> PluginOptions: ...

def parse_plugin_options(
    cfg: configparser.RawConfigParser, cfg_dir: str, *, enable_extensions: str | None, require_plugins: str | None
) -> PluginOptions: ...
def find_plugins(cfg: configparser.RawConfigParser, opts: PluginOptions) -> list[Plugin]: ...
def load_plugins(plugins: list[Plugin], opts: PluginOptions) -> Plugins: ...
