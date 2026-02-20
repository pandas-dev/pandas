from _typeshed import StrOrBytesPath, StrPath, Unused
from collections.abc import MutableSequence, Sequence
from subprocess import _ENV
from typing import Literal, overload

@overload
def spawn(
    cmd: Sequence[StrOrBytesPath], search_path: Literal[False], verbose: Unused = False, env: _ENV | None = None
) -> None: ...
@overload
def spawn(
    cmd: MutableSequence[bytes | StrPath], search_path: Literal[True] = True, verbose: Unused = False, env: _ENV | None = None
) -> None: ...
def find_executable(executable: str, path: str | None = None) -> str | None: ...
