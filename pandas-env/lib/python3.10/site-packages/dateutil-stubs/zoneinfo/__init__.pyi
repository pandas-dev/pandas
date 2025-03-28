from _typeshed import Incomplete
from typing import IO
from typing_extensions import TypeAlias

__all__ = ["get_zonefile_instance", "gettz", "gettz_db_metadata"]

_MetadataType: TypeAlias = dict[str, Incomplete]

class ZoneInfoFile:
    zones: dict[Incomplete, Incomplete]
    metadata: _MetadataType | None
    def __init__(self, zonefile_stream: IO[bytes] | None = None) -> None: ...
    def get(self, name, default: Incomplete | None = None): ...

def get_zonefile_instance(new_instance: bool = False) -> ZoneInfoFile: ...
def gettz(name): ...
def gettz_db_metadata() -> _MetadataType: ...
