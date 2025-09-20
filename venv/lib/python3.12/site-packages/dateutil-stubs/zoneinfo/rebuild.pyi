from _typeshed import StrOrBytesPath
from collections.abc import Iterable

def rebuild(filename: StrOrBytesPath, tag=None, format: str = "gz", zonegroups: Iterable[str] = [], metadata=None) -> None: ...
