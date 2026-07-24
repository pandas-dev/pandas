from typing import TypedDict, type_check_only

@type_check_only
class _Versions(TypedDict):
    dirty: bool
    error: None
    full_revisionid: str
    version: str

version_json: str

def get_versions() -> _Versions: ...
