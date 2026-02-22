from collections.abc import Callable, Iterable, Iterator, Mapping
from re import Pattern
from typing import Any, TypeVar, overload

from . import resolver as resolver  # Help mypy a bit; this is implied by loader and dumper
from .constructor import BaseConstructor
from .cyaml import *
from .cyaml import _CLoader
from .dumper import *
from .dumper import _Inf
from .emitter import _WriteStream
from .error import *
from .events import *
from .loader import *
from .loader import _Loader
from .nodes import *
from .reader import _ReadStream
from .representer import BaseRepresenter
from .resolver import BaseResolver
from .tokens import *

_T = TypeVar("_T")
_Constructor = TypeVar("_Constructor", bound=BaseConstructor)
_Representer = TypeVar("_Representer", bound=BaseRepresenter)

__with_libyaml__: bool
__version__: str

def warnings(settings=None): ...
def scan(stream, Loader: type[_Loader | _CLoader] = ...): ...
def parse(stream, Loader: type[_Loader | _CLoader] = ...): ...
def compose(stream, Loader: type[_Loader | _CLoader] = ...): ...
def compose_all(stream, Loader: type[_Loader | _CLoader] = ...): ...
def load(stream: _ReadStream, Loader: type[_Loader | _CLoader]) -> Any: ...
def load_all(stream: _ReadStream, Loader: type[_Loader | _CLoader]) -> Iterator[Any]: ...
def full_load(stream: _ReadStream) -> Any: ...
def full_load_all(stream: _ReadStream) -> Iterator[Any]: ...
def safe_load(stream: _ReadStream) -> Any: ...
def safe_load_all(stream: _ReadStream) -> Iterator[Any]: ...
def unsafe_load(stream: _ReadStream) -> Any: ...
def unsafe_load_all(stream: _ReadStream) -> Iterator[Any]: ...
def emit(
    events,
    stream: _WriteStream[Any] | None = None,
    Dumper=...,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
): ...
@overload
def serialize_all(
    nodes,
    stream: _WriteStream[Any],
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> None: ...
@overload
def serialize_all(
    nodes,
    stream: None = None,
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> str: ...
@overload
def serialize_all(
    nodes,
    stream: None = None,
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> bytes: ...
@overload
def serialize(
    node,
    stream: _WriteStream[Any],
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> None: ...
@overload
def serialize(
    node,
    stream: None = None,
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> str: ...
@overload
def serialize(
    node,
    stream: None = None,
    Dumper=...,
    *,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
) -> bytes: ...
@overload
def dump_all(
    documents: Iterable[Any],
    stream: _WriteStream[Any],
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> None: ...
@overload
def dump_all(
    documents: Iterable[Any],
    stream: None = None,
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> str: ...
@overload
def dump_all(
    documents: Iterable[Any],
    stream: None = None,
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> bytes: ...
@overload
def dump(
    data: Any,
    stream: _WriteStream[Any],
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> None: ...
@overload
def dump(
    data: Any,
    stream: None = None,
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> str: ...
@overload
def dump(
    data: Any,
    stream: None = None,
    Dumper=...,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> bytes: ...
@overload
def safe_dump_all(
    documents: Iterable[Any],
    stream: _WriteStream[Any],
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> None: ...
@overload
def safe_dump_all(
    documents: Iterable[Any],
    stream: None = None,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> str: ...
@overload
def safe_dump_all(
    documents: Iterable[Any],
    stream: None = None,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> bytes: ...
@overload
def safe_dump(
    data: Any,
    stream: _WriteStream[Any],
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str | None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> None: ...
@overload
def safe_dump(
    data: Any,
    stream: None = None,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: None = None,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> str: ...
@overload
def safe_dump(
    data: Any,
    stream: None = None,
    *,
    default_style: str | None = None,
    default_flow_style: bool | None = False,
    canonical: bool | None = None,
    indent: int | None = None,
    width: int | _Inf | None = None,
    allow_unicode: bool | None = None,
    line_break: str | None = None,
    encoding: str,
    explicit_start: bool | None = None,
    explicit_end: bool | None = None,
    version: tuple[int, int] | None = None,
    tags: Mapping[str, str] | None = None,
    sort_keys: bool = True,
) -> bytes: ...
def add_implicit_resolver(
    tag: str,
    regexp: Pattern[str],
    first: Iterable[Any] | None = None,
    Loader: type[BaseResolver] | None = None,
    Dumper: type[BaseResolver] = ...,
) -> None: ...
def add_path_resolver(
    tag: str,
    path: Iterable[Any],
    kind: type[Any] | None = None,
    Loader: type[BaseResolver] | None = None,
    Dumper: type[BaseResolver] = ...,
) -> None: ...
@overload
def add_constructor(
    tag: str, constructor: Callable[[Loader | FullLoader | UnsafeLoader, Node], Any], Loader: None = None
) -> None: ...
@overload
def add_constructor(tag: str, constructor: Callable[[_Constructor, Node], Any], Loader: type[_Constructor]) -> None: ...
@overload
def add_multi_constructor(
    tag_prefix: str, multi_constructor: Callable[[Loader | FullLoader | UnsafeLoader, str, Node], Any], Loader: None = None
) -> None: ...
@overload
def add_multi_constructor(
    tag_prefix: str, multi_constructor: Callable[[_Constructor, str, Node], Any], Loader: type[_Constructor]
) -> None: ...
@overload
def add_representer(data_type: type[_T], representer: Callable[[Dumper, _T], Node]) -> None: ...
@overload
def add_representer(data_type: type[_T], representer: Callable[[_Representer, _T], Node], Dumper: type[_Representer]) -> None: ...
@overload
def add_multi_representer(data_type: type[_T], multi_representer: Callable[[Dumper, _T], Node]) -> None: ...
@overload
def add_multi_representer(
    data_type: type[_T], multi_representer: Callable[[_Representer, _T], Node], Dumper: type[_Representer]
) -> None: ...

class YAMLObjectMetaclass(type):
    def __init__(cls, name, bases, kwds) -> None: ...

class YAMLObject(metaclass=YAMLObjectMetaclass):
    __slots__ = ()
    yaml_loader: Any
    yaml_dumper: Any
    yaml_tag: Any
    yaml_flow_style: Any
    @classmethod
    def from_yaml(cls, loader, node): ...
    @classmethod
    def to_yaml(cls, dumper, data): ...
