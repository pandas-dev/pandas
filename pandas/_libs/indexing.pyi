# flake8: noqa

from typing import Any

class _NDFrameIndexerBase:
    obj: Any
    name: Any
    def __init__(self, name, obj): ...
    @property
    def ndim(self): ...
