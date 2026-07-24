from typing import TypedDict

version_json: str
_Versions = TypedDict("_Versions", {"date": str, "dirty": bool, "error": None, "full-revisionid": str, "version": str})

def get_versions() -> _Versions: ...
