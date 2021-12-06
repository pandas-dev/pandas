import sys
from _typeshed import SupportsGetItem, SupportsItemAccess
from builtins import type as _type
from collections.abc import Iterable, Iterator, Mapping
from typing import IO, Any, Optional, Protocol, TypeVar, Union

_T = TypeVar("_T", bound=FieldStorage)

def parse(
    fp: Optional[IO[Any]] = ...,
    environ: SupportsItemAccess[str, str] = ...,
    keep_blank_values: bool = ...,
    strict_parsing: bool = ...,
    separator: str = ...,
) -> dict[str, list[str]]: ...

if sys.version_info < (3, 8):
    def parse_qs(qs: str, keep_blank_values: bool = ..., strict_parsing: bool = ...) -> dict[str, list[str]]: ...
    def parse_qsl(qs: str, keep_blank_values: bool = ..., strict_parsing: bool = ...) -> list[tuple[str, str]]: ...

if sys.version_info >= (3, 7):
    def parse_multipart(
        fp: IO[Any], pdict: SupportsGetItem[str, bytes], encoding: str = ..., errors: str = ..., separator: str = ...
    ) -> dict[str, list[Any]]: ...

else:
    def parse_multipart(fp: IO[Any], pdict: SupportsGetItem[str, bytes]) -> dict[str, list[bytes]]: ...

class _Environ(Protocol):
    def __getitem__(self, __k: str) -> str: ...
    def keys(self) -> Iterable[str]: ...

def parse_header(line: str) -> tuple[str, dict[str, str]]: ...
def test(environ: _Environ = ...) -> None: ...
def print_environ(environ: _Environ = ...) -> None: ...
def print_form(form: dict[str, Any]) -> None: ...
def print_directory() -> None: ...
def print_environ_usage() -> None: ...

if sys.version_info < (3, 8):
    def escape(s: str, quote: Optional[bool] = ...) -> str: ...

class MiniFieldStorage:
    # The first five "Any" attributes here are always None, but mypy doesn't support that
    filename: Any
    list: Any
    type: Any
    file: Optional[IO[bytes]]
    type_options: dict[Any, Any]
    disposition: Any
    disposition_options: dict[Any, Any]
    headers: dict[Any, Any]
    name: Any
    value: Any
    def __init__(self, name: Any, value: Any) -> None: ...
    def __repr__(self) -> str: ...

_list = list

class FieldStorage(object):
    FieldStorageClass: Optional[_type]
    keep_blank_values: int
    strict_parsing: int
    qs_on_post: Optional[str]
    headers: Mapping[str, str]
    fp: IO[bytes]
    encoding: str
    errors: str
    outerboundary: bytes
    bytes_read: int
    limit: Optional[int]
    disposition: str
    disposition_options: dict[str, str]
    filename: Optional[str]
    file: Optional[IO[bytes]]
    type: str
    type_options: dict[str, str]
    innerboundary: bytes
    length: int
    done: int
    list: Optional[_list[Any]]
    value: Union[None, bytes, _list[Any]]
    def __init__(
        self,
        fp: Optional[IO[Any]] = ...,
        headers: Optional[Mapping[str, str]] = ...,
        outerboundary: bytes = ...,
        environ: SupportsGetItem[str, str] = ...,
        keep_blank_values: int = ...,
        strict_parsing: int = ...,
        limit: Optional[int] = ...,
        encoding: str = ...,
        errors: str = ...,
        max_num_fields: Optional[int] = ...,
        separator: str = ...,
    ) -> None: ...
    def __enter__(self: _T) -> _T: ...
    def __exit__(self, *args: Any) -> None: ...
    def __repr__(self) -> str: ...
    def __iter__(self) -> Iterator[str]: ...
    def __getitem__(self, key: str) -> Any: ...
    def getvalue(self, key: str, default: Any = ...) -> Any: ...
    def getfirst(self, key: str, default: Any = ...) -> Any: ...
    def getlist(self, key: str) -> _list[Any]: ...
    def keys(self) -> _list[str]: ...
    def __contains__(self, key: str) -> bool: ...
    def __len__(self) -> int: ...
    def __bool__(self) -> bool: ...
    # In Python 3 it returns bytes or str IO depending on an internal flag
    def make_file(self) -> IO[Any]: ...
