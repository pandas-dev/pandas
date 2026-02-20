import enum
from typing import Any, Final

from ._array_api import _CapabilitiesTable  # type-check-only

###

ALIASES: Final[dict[str, set[str]]] = ...
BACKEND_NAMES_MAP: Final[dict[str, str]] = ...

class BackendSupportStatus(enum.Enum):
    YES = 1
    NO = 2
    OUT_OF_SCOPE = 3
    UNKNOWN = 4

def _process_capabilities_table_entry(entry: dict[str, Any] | None) -> dict[str, dict[str, bool]]: ...
def is_named_function_like_object(obj: object) -> bool: ...
def make_flat_capabilities_table(
    modules: str | list[str], backend_type: str, /, *, capabilities_table: _CapabilitiesTable | None = None
) -> list[dict[str, str]]: ...
def calculate_table_statistics(flat_table: list[dict[str, str]]) -> dict[str, tuple[dict[str, str], bool]]: ...
