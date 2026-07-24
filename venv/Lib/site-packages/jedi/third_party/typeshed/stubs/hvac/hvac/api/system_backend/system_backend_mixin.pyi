import logging
from abc import ABCMeta

from hvac.api.vault_api_base import VaultApiBase

logger: logging.Logger

class SystemBackendMixin(VaultApiBase, metaclass=ABCMeta): ...
