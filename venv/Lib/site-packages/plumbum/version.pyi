from typing_extensions import TypeAlias

__all__ = [
    "__commit_id__",
    "__version__",
    "__version_tuple__",
    "commit_id",
    "version",
    "version_tuple",
]

VERSION_TUPLE: TypeAlias = tuple[int | str, ...]
COMMIT_ID: TypeAlias = str | None

version: str
__version__: str
__version_tuple__: VERSION_TUPLE
version_tuple: VERSION_TUPLE
commit_id: COMMIT_ID
__commit_id__: COMMIT_ID
