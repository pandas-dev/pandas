from abc import ABCMeta, abstractmethod
from logging import Logger
from typing import Any

from hvac.adapters import Adapter
from hvac.api.vault_api_base import VaultApiBase

logger: Logger

class VaultApiCategory(VaultApiBase, metaclass=ABCMeta):
    implemented_class_names: list[str]
    def __init__(self, adapter: Adapter[Any]) -> None: ...
    def __getattr__(self, item): ...
    @property
    def adapter(self) -> Adapter[Any]: ...
    @adapter.setter
    def adapter(self, adapter: Adapter[Any]) -> None: ...
    @property
    @abstractmethod
    def implemented_classes(self): ...
    @property
    def unimplemented_classes(self) -> list[str]: ...
    @staticmethod
    def get_private_attr_name(class_name): ...
