from logging import Logger
from typing import ClassVar

log: Logger

class XRayConfig:
    name: ClassVar[str]
    def ready(self) -> None: ...
