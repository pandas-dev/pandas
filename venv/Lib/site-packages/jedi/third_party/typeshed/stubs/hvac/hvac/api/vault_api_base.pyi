from _typeshed import Incomplete
from abc import ABCMeta
from logging import Logger

from hvac.adapters import Adapter

logger: Logger

class VaultApiBase(metaclass=ABCMeta):
    def __init__(self, adapter: Adapter[Incomplete]) -> None: ...
