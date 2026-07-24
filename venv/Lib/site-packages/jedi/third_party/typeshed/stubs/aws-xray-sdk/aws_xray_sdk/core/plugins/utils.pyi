from collections.abc import Iterable
from types import ModuleType
from typing import Final

module_prefix: Final[str]
PLUGIN_MAPPING: Final[dict[str, str]]

def get_plugin_modules(plugins: Iterable[str]) -> tuple[ModuleType, ...]: ...
