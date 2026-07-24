from __future__ import annotations

import logging
from tempfile import mkdtemp
from typing import TYPE_CHECKING

from virtualenv.util.path import safe_delete

from .via_disk_folder import AppDataDiskFolder

if TYPE_CHECKING:
    from typing import NoReturn

LOGGER = logging.getLogger(__name__)


class TempAppData(AppDataDiskFolder):
    transient = True
    can_update = False

    def __init__(self) -> None:
        super().__init__(folder=mkdtemp())
        LOGGER.debug("created temporary app data folder %s", self.lock.path)

    def reset(self) -> None:
        """This is a temporary folder, is already empty to start with."""

    def close(self) -> None:
        LOGGER.debug("remove temporary app data folder %s", self.lock.path)
        safe_delete(self.lock.path)

    def embed_update_log(self, distribution: str, for_py_version: str) -> NoReturn:
        raise NotImplementedError


__all__ = [
    "TempAppData",
]
