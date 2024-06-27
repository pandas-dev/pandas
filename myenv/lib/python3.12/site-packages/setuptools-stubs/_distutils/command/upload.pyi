from _typeshed import Incomplete, StrOrBytesPath
from collections.abc import Iterable, Iterator
from typing import AnyStr, ClassVar, TypeVar, overload

from ..config import PyPIRCCommand

_T = TypeVar("_T")

class upload(PyPIRCCommand):
    description: ClassVar[str]
    username: str
    password: str
    show_response: int
    sign: bool
    identity: Incomplete
    def initialize_options(self) -> None: ...
    repository: Incomplete
    realm: Incomplete
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
    def upload_file(self, command: str, pyversion: str, filename: StrOrBytesPath) -> None: ...

@overload
def make_iterable(values: None) -> list[None]: ...
@overload
def make_iterable(values: AnyStr) -> Iterator[AnyStr]: ...
@overload
def make_iterable(values: Iterable[_T]) -> Iterator[_T]: ...
@overload
def make_iterable(values: _T) -> Iterator[_T]: ...
