from collections.abc import Callable, Container, Iterable, Iterator
from re import Pattern
from typing import Final, Protocol, type_check_only
from typing_extensions import TypeAlias

from html5lib.filters.base import Filter
from html5lib.filters.sanitizer import Filter as SanitizerFilter
from html5lib.treewalkers.base import TreeWalker

from . import _HTMLAttrKey
from .css_sanitizer import CSSSanitizer
from .html5lib_shim import BleachHTMLParser, BleachHTMLSerializer
from .linkifier import _Token

ALLOWED_TAGS: Final[frozenset[str]]
ALLOWED_ATTRIBUTES: Final[dict[str, list[str]]]
ALLOWED_PROTOCOLS: Final[frozenset[str]]

INVISIBLE_CHARACTERS: Final[str]
INVISIBLE_CHARACTERS_RE: Final[Pattern[str]]
INVISIBLE_REPLACEMENT_CHAR: Final = "?"

class NoCssSanitizerWarning(UserWarning): ...

@type_check_only
class _FilterConstructor(Protocol):
    def __call__(self, *, source: BleachSanitizerFilter) -> Filter: ...

# _FilterConstructor used to be called _Filter
# this alias is obsolete and can potentially be removed in the future
_Filter: TypeAlias = _FilterConstructor  # noqa: Y047

_AttributeFilter: TypeAlias = Callable[[str, str, str], bool]
_AttributeDict: TypeAlias = dict[str, list[str] | _AttributeFilter] | dict[str, list[str]] | dict[str, _AttributeFilter]
_Attributes: TypeAlias = _AttributeFilter | _AttributeDict | list[str]

class Cleaner:
    tags: Iterable[str]
    attributes: _Attributes
    protocols: Iterable[str]
    strip: bool
    strip_comments: bool
    filters: Iterable[_FilterConstructor]
    css_sanitizer: CSSSanitizer | None
    parser: BleachHTMLParser
    walker: TreeWalker
    serializer: BleachHTMLSerializer
    def __init__(
        self,
        tags: Iterable[str] = ...,
        attributes: _Attributes = ...,
        protocols: Iterable[str] = ...,
        strip: bool = False,
        strip_comments: bool = True,
        filters: Iterable[_FilterConstructor] | None = None,
        css_sanitizer: CSSSanitizer | None = None,
    ) -> None: ...
    def clean(self, text: str) -> str: ...

def attribute_filter_factory(attributes: _Attributes) -> _AttributeFilter: ...

class BleachSanitizerFilter(SanitizerFilter):
    allowed_tags: frozenset[str]
    allowed_protocols: frozenset[str]
    attr_filter: _AttributeFilter
    strip_disallowed_tags: bool
    strip_html_comments: bool
    attr_val_is_uri: frozenset[_HTMLAttrKey]
    svg_attr_val_allows_ref: frozenset[_HTMLAttrKey]
    svg_allow_local_href: frozenset[_HTMLAttrKey]
    css_sanitizer: CSSSanitizer | None
    def __init__(
        self,
        source: TreeWalker,
        allowed_tags: Iterable[str] = ...,
        attributes: _Attributes = ...,
        allowed_protocols: Iterable[str] = ...,
        attr_val_is_uri: frozenset[_HTMLAttrKey] = ...,
        svg_attr_val_allows_ref: frozenset[_HTMLAttrKey] = ...,
        svg_allow_local_href: frozenset[_HTMLAttrKey] = ...,
        strip_disallowed_tags: bool = False,
        strip_html_comments: bool = True,
        css_sanitizer: CSSSanitizer | None = None,
    ) -> None: ...
    def sanitize_stream(self, token_iterator: Iterable[_Token]) -> Iterator[_Token]: ...
    def merge_characters(self, token_iterator: Iterable[_Token]) -> Iterator[_Token]: ...
    def __iter__(self) -> Iterator[_Token]: ...
    def sanitize_token(self, token: _Token) -> _Token | list[_Token] | None: ...  # type: ignore[override]
    def sanitize_characters(self, token: _Token) -> _Token | list[_Token]: ...
    def sanitize_uri_value(self, value: str, allowed_protocols: Container[str]) -> str | None: ...
    def allow_token(self, token: _Token) -> _Token: ...
    def disallowed_token(self, token: _Token) -> _Token: ...
