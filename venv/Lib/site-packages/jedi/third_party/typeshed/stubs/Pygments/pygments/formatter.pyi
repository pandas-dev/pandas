import types
from _typeshed import SupportsWrite
from collections.abc import Iterable, Sequence
from typing import Any, ClassVar, Generic, TypeVar, overload

from pygments.style import Style
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["Formatter"]

class Formatter(Generic[_T]):
    name: ClassVar[str]  # Set to None, but always overridden with a non-None value in subclasses.
    aliases: ClassVar[Sequence[str]]  # Not intended to be mutable
    filenames: ClassVar[Sequence[str]]  # Not intended to be mutable
    unicodeoutput: ClassVar[bool]
    style: type[Style]
    full: bool
    title: str
    encoding: str | None
    options: dict[str, Any]  # arbitrary values used by subclasses
    @overload
    def __init__(
        self: Formatter[str],
        *,
        style: type[Style] | str = "default",
        full: bool = False,
        title: str = "",
        encoding: None = None,
        outencoding: None = None,
        **options: Any,  # arbitrary values used by subclasses
    ) -> None: ...
    @overload
    def __init__(
        self: Formatter[bytes],
        *,
        style: type[Style] | str = "default",
        full: bool = False,
        title: str = "",
        encoding: str,
        outencoding: None = None,
        **options: Any,  # arbitrary values used by subclasses
    ) -> None: ...
    @overload
    def __init__(
        self: Formatter[bytes],
        *,
        style: type[Style] | str = "default",
        full: bool = False,
        title: str = "",
        encoding: None = None,
        outencoding: str,
        **options: Any,  # arbitrary values used by subclasses
    ) -> None: ...
    def __class_getitem__(cls, name: Any) -> types.GenericAlias: ...
    def get_style_defs(self, arg: str = "") -> str: ...
    def format(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[_T]) -> None: ...
