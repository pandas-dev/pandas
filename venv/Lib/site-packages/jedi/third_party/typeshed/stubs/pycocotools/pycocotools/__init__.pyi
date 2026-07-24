from typing import TypedDict, type_check_only

# Unused in this module, but imported in multiple submodules.
@type_check_only
class _EncodedRLE(TypedDict):  # noqa: Y049
    size: list[int]
    counts: str | bytes
