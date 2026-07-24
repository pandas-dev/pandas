from typing import Literal

from gdb import Progspace

class MissingFileHandler:
    @property
    def name(self) -> str: ...
    enabled: bool
    def __init__(self, name: str) -> None: ...

def register_handler(
    handler_type: Literal["debug", "objfile"], locus: Progspace | None, handler: MissingFileHandler, replace: bool = False
) -> None: ...
