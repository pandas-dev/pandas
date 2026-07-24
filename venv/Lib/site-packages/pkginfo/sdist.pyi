from typing import Sequence
from .distribution import Distribution as Distribution

class NoSuchFile(ValueError):
    fqp: str
    def __init__(self, fqp: str) -> None: ...

class UnknownArchiveFormat(ValueError):
    fqp: str
    def __init__(self, fqp: str) -> None: ...

class InvalidPkgInfo(ValueError):
    fqp: str
    candidates: Sequence[Sequence[str]]
    def __init__(self, fqp: str, candidates: Sequence[Sequence[str]]) -> None:
        ...

class NoPkgInfo(ValueError):
    fqp: str
    def __init__(self, fqp: str) -> None: ...

class InvalidUnpackedSDist(ValueError):
    fqp: str
    def __init__(self, fqp: str, raised: Exception) -> None: ...

class SDist(Distribution):
    filename: str
    metadata_version: str | None
    def __init__(self, filename: str, metadata_version: str | None = ...) -> None: ...
    def read(self) -> bytes: ...

class UnpackedSDist(SDist):
    def __init__(self, filename: str, metadata_version: str | None = ...) -> None: ...
    def read(self) -> bytes: ...
