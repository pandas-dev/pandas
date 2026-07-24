from typing import Final, Literal
from typing_extensions import TypeAlias

_Mode: TypeAlias = Literal[
    "GAE-LOCAL", "GAE-SERVER", "MOD_WSGI", "INTERACTIVE", "FCGI-FLUP", "UWSGI", "FLASK", "CHERRYPY", "BOTTLE", "UNKNOWN"
]
__version__: Final[str]

def detect_mode() -> _Mode: ...

MODE: Final[_Mode]
MAIN_FILE: Final[str | None]
MAIN_DIR: Final[str | None]
PONY_DIR: Final[str]
