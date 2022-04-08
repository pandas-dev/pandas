from typing import Any

class NDFrameIndexerBase:
    name: str
    obj: Any

    def __init__(self, name: str, obj: Any) -> None: ...
    @property
    def ndim(self) -> int: ...
