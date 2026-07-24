from typing import Any

from google._upb._message import MessageMeta

MESSAGE_CLASS_CACHE: dict[str, Any]

class GeneratedProtocolMessageType(MessageMeta):
    def __new__(  # noqa: Y034
        cls, name: str, bases: tuple[type, ...], dictionary: dict[str, Any]
    ) -> GeneratedProtocolMessageType: ...
