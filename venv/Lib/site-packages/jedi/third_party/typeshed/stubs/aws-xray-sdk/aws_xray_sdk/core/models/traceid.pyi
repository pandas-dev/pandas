from typing import ClassVar

class TraceId:
    VERSION: ClassVar[str]
    DELIMITER: ClassVar[str]
    start_time: int
    def __init__(self) -> None: ...
    def to_id(self) -> str: ...
