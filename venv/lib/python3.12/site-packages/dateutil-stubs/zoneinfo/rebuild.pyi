from _typeshed import StrOrBytesPath, Unused
from collections.abc import Iterable

from ..zoneinfo import _MetadataType

def rebuild(
    filename: StrOrBytesPath,
    tag: Unused | None = None,
    format: str = "gz",
    zonegroups: Iterable[str] = [],
    metadata: _MetadataType | None = None,
) -> None: ...
