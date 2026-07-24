from logging import Logger

from .connector import ServiceConnector
from .rule_cache import RuleCache
from .rule_poller import RulePoller

log: Logger

class TargetPoller:
    def __init__(self, cache: RuleCache, rule_poller: RulePoller, connector: ServiceConnector) -> None: ...
    def start(self) -> None: ...
