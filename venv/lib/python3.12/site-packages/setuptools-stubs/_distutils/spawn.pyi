from _typeshed import StrPath
from collections.abc import MutableSequence
from subprocess import _ENV

def spawn(
    cmd: MutableSequence[bytes | StrPath],
    search_path: bool = True,
    verbose: bool = False,
    dry_run: bool = False,
    env: _ENV | None = None,
) -> None: ...
def find_executable(executable: str, path: str | None = None) -> str | None: ...
