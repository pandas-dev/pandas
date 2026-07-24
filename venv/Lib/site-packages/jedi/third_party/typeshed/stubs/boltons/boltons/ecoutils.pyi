from typing import Any

ECO_VERSION: str
HAVE_URANDOM: bool
INSTANCE_ID: str
IS_64BIT: bool
HAVE_UCS4: bool
HAVE_READLINE: bool
SQLITE_VERSION: str
OPENSSL_VERSION: str
TKINTER_VERSION: str
ZLIB_VERSION: str
EXPAT_VERSION: str
CPU_COUNT: int
HAVE_THREADING: bool
HAVE_IPV6: bool
RLIMIT_FDS_SOFT: int
RLIMIT_FDS_HARD: int
START_TIME_INFO: dict[str, str | float]

def getrandbits(k: int) -> int: ...
def get_python_info() -> dict[str, Any]: ...
def get_profile(**kwargs) -> dict[str, Any]: ...
def get_profile_json(indent: bool = False) -> str: ...
def main() -> None: ...
def dumps(val: Any, indent: int) -> str: ...
