from typing import Any, Literal

from .YoutubeDL import YoutubeDL

class Cache:
    def __init__(self, ydl: YoutubeDL) -> None: ...
    @property
    def enabled(self) -> bool: ...
    def store(
        self, section: str, key: str, data: Any, dtype: Literal["json"] = "json"  # data is anything JSON serializable.
    ) -> None: ...
    def load(
        self,
        section: str,
        key: str,
        dtype: Literal["json"] = "json",
        default: Any = None,  # Returned if not enabled or if the cache entry is not found.
        *,
        min_ver: str | None = None,
    ) -> Any: ...  # Anything JSON serializable.
    def remove(self) -> None: ...
