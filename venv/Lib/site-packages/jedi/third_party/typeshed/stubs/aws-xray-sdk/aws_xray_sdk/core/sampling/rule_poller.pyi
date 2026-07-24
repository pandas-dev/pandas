from logging import Logger
from typing import Final

from .connector import ServiceConnector
from .rule_cache import RuleCache

log: Logger
DEFAULT_INTERVAL: Final = 300

class RulePoller:
    def __init__(self, cache: RuleCache, connector: ServiceConnector) -> None: ...
    def start(self) -> None: ...
    def wake_up(self) -> None: ...
