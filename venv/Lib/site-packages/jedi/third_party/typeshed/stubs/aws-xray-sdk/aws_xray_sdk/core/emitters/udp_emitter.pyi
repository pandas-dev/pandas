from logging import Logger
from typing import Final

from aws_xray_sdk.core.models.entity import Entity

log: Logger
PROTOCOL_HEADER: Final[str]
PROTOCOL_DELIMITER: Final[str]
DEFAULT_DAEMON_ADDRESS: Final[str]

class UDPEmitter:
    def __init__(self, daemon_address: str = "127.0.0.1:2000") -> None: ...
    def send_entity(self, entity: Entity) -> None: ...
    def set_daemon_address(self, address: str | None) -> None: ...
    @property
    def ip(self) -> str: ...
    @property
    def port(self) -> int: ...
