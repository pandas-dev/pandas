from _typeshed import GenericPath, StrOrBytesPath, StrPath
from collections.abc import Iterator
from typing import AnyStr

VCS_NAMES: list[str]
IGNORED_EXTENSIONS: list[str]

def list_directory(path: GenericPath[AnyStr], include_directories: bool = True) -> Iterator[AnyStr]: ...
def mtime(path: StrOrBytesPath) -> str: ...
def md5(path: StrPath) -> str: ...
