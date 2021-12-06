import sys
from _typeshed import AnyPath, StrPath
from collections.abc import Callable
from configparser import RawConfigParser
from threading import Thread
from typing import IO, Any, Optional, Pattern, Union

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

if sys.version_info >= (3, 7):
    _Path = AnyPath
else:
    _Path = StrPath

DEFAULT_LOGGING_CONFIG_PORT: int
RESET_ERROR: int  # undocumented
IDENTIFIER: Pattern[str]  # undocumented

def dictConfig(config: dict[str, Any]) -> None: ...

if sys.version_info >= (3, 10):
    def fileConfig(
        fname: Union[_Path, IO[str], RawConfigParser],
        defaults: Optional[dict[str, str]] = ...,
        disable_existing_loggers: bool = ...,
        encoding: Optional[str] = ...,
    ) -> None: ...

else:
    def fileConfig(
        fname: Union[_Path, IO[str], RawConfigParser],
        defaults: Optional[dict[str, str]] = ...,
        disable_existing_loggers: bool = ...,
    ) -> None: ...

def valid_ident(s: str) -> Literal[True]: ...  # undocumented
def listen(port: int = ..., verify: Optional[Callable[[bytes], Optional[bytes]]] = ...) -> Thread: ...
def stopListening() -> None: ...
