from logging import Logger
from typing import ClassVar

log: Logger

class SDKConfig:
    XRAY_ENABLED_KEY: ClassVar[str]
    DISABLED_ENTITY_NAME: ClassVar[str]
    @classmethod
    def sdk_enabled(cls) -> bool: ...
    @classmethod
    def set_sdk_enabled(cls, value: bool | None) -> None: ...
