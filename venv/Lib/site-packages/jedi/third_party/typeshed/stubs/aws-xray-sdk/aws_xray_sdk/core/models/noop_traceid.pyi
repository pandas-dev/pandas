from typing import ClassVar

class NoOpTraceId:
    VERSION: ClassVar[str]
    DELIMITER: ClassVar[str]
    start_time: str
    def __init__(self) -> None: ...
    def to_id(self) -> str: ...
