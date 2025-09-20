from _typeshed import Incomplete
from collections.abc import Callable
from io import BytesIO
from tarfile import _Fileobj
from typing import Final, TypeVar, overload
from typing_extensions import Self, TypeAlias, deprecated

from dateutil.tz import tzfile as _tzfile

_T = TypeVar("_T")
_MetadataType: TypeAlias = dict[str, Incomplete]

__all__ = ["get_zonefile_instance", "gettz", "gettz_db_metadata"]

ZONEFILENAME: Final[str]
METADATA_FN: Final[str]

class tzfile(_tzfile):
    def __reduce__(self) -> tuple[Callable[[str], Self], tuple[str, ...]]: ...

def getzoneinfofile_stream() -> BytesIO | None: ...

class ZoneInfoFile:
    zones: dict[str, _tzfile]
    metadata: _MetadataType | None
    def __init__(self, zonefile_stream: _Fileobj | None = None) -> None: ...
    @overload
    def get(self, name: str, default: None = None) -> _tzfile | None: ...
    @overload
    def get(self, name: str, default: _tzfile) -> _tzfile: ...
    @overload
    def get(self, name: str, default: _T) -> _tzfile | _T: ...

def get_zonefile_instance(new_instance: bool = False) -> ZoneInfoFile: ...
@deprecated(
    "zoneinfo.gettz() will be removed in future versions, to use the dateutil-provided "
    "zoneinfo files, instantiate a ZoneInfoFile object and use ZoneInfoFile.zones.get() instead."
)
def gettz(name: str) -> _tzfile: ...
@deprecated(
    "zoneinfo.gettz_db_metadata() will be removed in future versions, to use the "
    "dateutil-provided zoneinfo files, ZoneInfoFile object and query the 'metadata' attribute instead."
)
def gettz_db_metadata() -> _MetadataType: ...
