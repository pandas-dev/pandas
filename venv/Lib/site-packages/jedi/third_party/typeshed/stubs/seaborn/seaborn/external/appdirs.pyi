from typing import Literal

__version__: str
__version_info__: tuple[int, int, int]
system: str

def user_cache_dir(
    appname: str | None = None, appauthor: Literal[False] | str | None = None, version: str | None = None, opinion: bool = True
) -> str: ...
