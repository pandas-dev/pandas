from _typeshed import SupportsItems
from typing import TypeVar

_KT = TypeVar("_KT")
_VT = TypeVar("_VT")

BOOL_CONFIG: set[str]

def asbool(obj: object) -> bool: ...
def convert_config(config: SupportsItems[_KT, _VT]) -> dict[_KT, _VT | bool]: ...
