from typing import NamedTuple

__all__ = ["RL_Codecs"]

class StdCodecData(NamedTuple):
    exceptions: dict[int, int | None] | None
    rexceptions: dict[int, int | None] | None

class ExtCodecData(NamedTuple):
    baseName: str
    exceptions: dict[int, int | None] | None
    rexceptions: dict[int, int | None] | None

class RL_Codecs:
    def __init__(self) -> None: ...
    @staticmethod
    def register() -> None: ...
    @staticmethod
    def add_dynamic_codec(name: str, exceptions, rexceptions) -> None: ...
    @staticmethod
    def remove_dynamic_codec(name: str) -> None: ...
    @staticmethod
    def reset_dynamic_codecs() -> None: ...
