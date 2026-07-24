from _typeshed import SupportsItems
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any, Literal, NamedTuple, Protocol, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

from webob._types import AsymmetricPropertyWithDelete

_T = TypeVar("_T")
_ListOrTuple: TypeAlias = list[_T] | tuple[_T, ...]
_ParsedAccept: TypeAlias = tuple[str, float, list[tuple[str, str]], list[str | tuple[str, str]]]

@type_check_only
class _SupportsStr(Protocol):
    def __str__(self) -> str: ...  # noqa: Y029

_AnyAcceptHeader: TypeAlias = AcceptValidHeader | AcceptInvalidHeader | AcceptNoHeader
_AnyAcceptCharsetHeader: TypeAlias = AcceptCharsetValidHeader | AcceptCharsetInvalidHeader | AcceptCharsetNoHeader
_AnyAcceptEncodingHeader: TypeAlias = AcceptEncodingValidHeader | AcceptEncodingInvalidHeader | AcceptEncodingNoHeader
_AnyAcceptLanguageHeader: TypeAlias = AcceptLanguageValidHeader | AcceptLanguageInvalidHeader | AcceptLanguageNoHeader

_AcceptProperty: TypeAlias = AsymmetricPropertyWithDelete[
    _AnyAcceptHeader,
    (
        _AnyAcceptHeader
        | SupportsItems[str, float | tuple[float, str]]
        | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
        | _SupportsStr
        | str
        | None
    ),
]
_AcceptCharsetProperty: TypeAlias = AsymmetricPropertyWithDelete[
    _AnyAcceptCharsetHeader,
    (
        _AnyAcceptCharsetHeader
        | SupportsItems[str, float]
        | _ListOrTuple[str | tuple[str, float] | list[Any]]
        | _SupportsStr
        | str
        | None
    ),
]
_AcceptEncodingProperty: TypeAlias = AsymmetricPropertyWithDelete[
    _AnyAcceptEncodingHeader,
    (
        _AnyAcceptEncodingHeader
        | SupportsItems[str, float]
        | _ListOrTuple[str | tuple[str, float] | list[Any]]
        | _SupportsStr
        | str
        | None
    ),
]
_AcceptLanguageProperty: TypeAlias = AsymmetricPropertyWithDelete[
    _AnyAcceptLanguageHeader,
    (
        _AnyAcceptLanguageHeader
        | SupportsItems[str, float]
        | _ListOrTuple[str | tuple[str, float] | list[Any]]
        | _SupportsStr
        | str
        | None
    ),
]

@type_check_only
class _AcceptOffer(NamedTuple):
    type: str
    subtype: str
    params: tuple[tuple[str, str], ...]

class AcceptOffer(_AcceptOffer):
    __slots__ = ()

class Accept:
    @classmethod
    def parse(cls, value: str) -> Iterator[_ParsedAccept]: ...
    @classmethod
    def parse_offer(cls, offer: str | AcceptOffer) -> AcceptOffer: ...

class AcceptValidHeader(Accept):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> list[_ParsedAccept]: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    def __add__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def __bool__(self) -> Literal[True]: ...
    def __contains__(self, offer: str) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __radd__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def accept_html(self) -> bool: ...
    @property
    def accepts_html(self) -> bool: ...
    def acceptable_offers(self, offers: Sequence[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float | None: ...

class _AcceptInvalidOrNoHeader(Accept):
    def __bool__(self) -> Literal[False]: ...
    def __contains__(self, offer: str) -> Literal[True]: ...
    def __iter__(self) -> Iterator[str]: ...
    def accept_html(self) -> bool: ...
    @property
    def accepts_html(self) -> bool: ...
    def acceptable_offers(self, offers: Sequence[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float: ...

class AcceptNoHeader(_AcceptInvalidOrNoHeader):
    @property
    def header_value(self) -> None: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptValidHeader | Literal[""]) -> AcceptValidHeader: ...
    @overload
    def __add__(self, other: AcceptNoHeader | AcceptInvalidHeader | None) -> Self: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptValidHeader: ...
    @overload
    def __radd__(self, other: AcceptValidHeader | Literal[""]) -> AcceptValidHeader: ...
    @overload
    def __radd__(self, other: AcceptNoHeader | AcceptInvalidHeader | None) -> Self: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptValidHeader: ...

class AcceptInvalidHeader(_AcceptInvalidOrNoHeader):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptValidHeader | Literal[""]) -> AcceptValidHeader: ...
    @overload
    def __add__(self, other: AcceptInvalidHeader | AcceptNoHeader | None) -> AcceptNoHeader: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptValidHeader | AcceptNoHeader: ...
    @overload
    def __radd__(self, other: AcceptValidHeader | Literal[""]) -> AcceptValidHeader: ...
    @overload
    def __radd__(self, other: AcceptInvalidHeader | AcceptNoHeader | None) -> AcceptNoHeader: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptHeader
            | SupportsItems[str, float | tuple[float, str]]
            | _ListOrTuple[str | tuple[str, float, str] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptValidHeader | AcceptNoHeader: ...

@overload
def create_accept_header(header_value: AcceptValidHeader | Literal[""]) -> AcceptValidHeader: ...
@overload
def create_accept_header(header_value: AcceptInvalidHeader) -> AcceptInvalidHeader: ...
@overload
def create_accept_header(header_value: None | AcceptNoHeader) -> AcceptNoHeader: ...
@overload
def create_accept_header(header_value: str) -> AcceptValidHeader | AcceptInvalidHeader: ...
@overload
def create_accept_header(header_value: _AnyAcceptHeader | str | None) -> _AnyAcceptHeader: ...
def accept_property() -> _AcceptProperty: ...

class AcceptCharset:
    @classmethod
    def parse(cls, value: str) -> Iterator[tuple[str, float]]: ...

class AcceptCharsetValidHeader(AcceptCharset):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> list[tuple[str, float]]: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    def __add__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def __bool__(self) -> Literal[True]: ...
    def __contains__(self, offer: str) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __radd__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def acceptable_offers(self, offers: Sequence[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float | None: ...

class _AcceptCharsetInvalidOrNoHeader(AcceptCharset):
    def __bool__(self) -> Literal[False]: ...
    def __contains__(self, offer: str) -> Literal[True]: ...
    def __iter__(self) -> Iterator[str]: ...
    def acceptable_offers(self, offers: Iterable[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float | None: ...

class AcceptCharsetNoHeader(_AcceptCharsetInvalidOrNoHeader):
    @property
    def header_value(self) -> None: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptCharsetValidHeader) -> AcceptCharsetValidHeader: ...
    @overload
    def __add__(self, other: AcceptCharsetInvalidHeader | AcceptCharsetNoHeader | Literal[""] | None) -> Self: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptCharsetValidHeader: ...
    @overload
    def __radd__(self, other: AcceptCharsetValidHeader) -> AcceptCharsetValidHeader: ...
    @overload
    def __radd__(self, other: AcceptCharsetInvalidHeader | AcceptCharsetNoHeader | Literal[""] | None) -> Self: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptCharsetValidHeader: ...

class AcceptCharsetInvalidHeader(_AcceptCharsetInvalidOrNoHeader):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptCharsetValidHeader) -> AcceptCharsetValidHeader: ...
    @overload
    def __add__(
        self, other: AcceptCharsetInvalidHeader | AcceptCharsetNoHeader | Literal[""] | None
    ) -> AcceptCharsetNoHeader: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptCharsetValidHeader | AcceptCharsetNoHeader: ...
    @overload
    def __radd__(self, other: AcceptCharsetValidHeader) -> AcceptCharsetValidHeader: ...
    @overload
    def __radd__(
        self, other: AcceptCharsetInvalidHeader | AcceptCharsetNoHeader | Literal[""] | None
    ) -> AcceptCharsetNoHeader: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptCharsetHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptCharsetValidHeader | AcceptCharsetNoHeader: ...

@overload
def create_accept_charset_header(header_value: AcceptCharsetValidHeader | Literal[""]) -> AcceptCharsetValidHeader: ...
@overload
def create_accept_charset_header(header_value: AcceptCharsetInvalidHeader) -> AcceptCharsetInvalidHeader: ...
@overload
def create_accept_charset_header(header_value: AcceptCharsetNoHeader | None) -> AcceptCharsetNoHeader: ...
@overload
def create_accept_charset_header(header_value: str) -> AcceptCharsetValidHeader | AcceptCharsetInvalidHeader: ...
@overload
def create_accept_charset_header(header_value: _AnyAcceptCharsetHeader | str | None) -> _AnyAcceptCharsetHeader: ...
def accept_charset_property() -> _AcceptCharsetProperty: ...

class AcceptEncoding:
    @classmethod
    def parse(cls, value: str) -> Iterator[tuple[str, float]]: ...

class AcceptEncodingValidHeader(AcceptEncoding):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> list[tuple[str, float]]: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    def __add__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def __bool__(self) -> Literal[True]: ...
    def __contains__(self, offer: str) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __radd__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def acceptable_offers(self, offers: Sequence[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float | None: ...

class _AcceptEncodingInvalidOrNoHeader(AcceptEncoding):
    def __bool__(self) -> Literal[False]: ...
    def __contains__(self, offer: str) -> Literal[True]: ...
    def __iter__(self) -> Iterator[str]: ...
    def acceptable_offers(self, offers: Iterable[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    def quality(self, offer: str) -> float | None: ...

class AcceptEncodingNoHeader(_AcceptEncodingInvalidOrNoHeader):
    @property
    def header_value(self) -> None: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptEncodingValidHeader | Literal[""]) -> AcceptEncodingValidHeader: ...
    @overload
    def __add__(self, other: AcceptEncodingInvalidHeader | AcceptEncodingNoHeader | None) -> Self: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptEncodingValidHeader: ...
    @overload
    def __radd__(self, other: AcceptEncodingValidHeader | Literal[""]) -> AcceptEncodingValidHeader: ...
    @overload
    def __radd__(self, other: AcceptEncodingInvalidHeader | AcceptEncodingNoHeader | None) -> Self: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptEncodingValidHeader: ...

class AcceptEncodingInvalidHeader(_AcceptEncodingInvalidOrNoHeader):
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> None: ...
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    @overload
    def __add__(self, other: AcceptEncodingValidHeader | Literal[""]) -> AcceptEncodingValidHeader: ...
    @overload
    def __add__(self, other: AcceptEncodingInvalidHeader | AcceptEncodingNoHeader | None) -> AcceptEncodingNoHeader: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptEncodingValidHeader | AcceptEncodingNoHeader: ...
    @overload
    def __radd__(self, other: AcceptEncodingValidHeader | Literal[""]) -> AcceptEncodingValidHeader: ...
    @overload
    def __radd__(self, other: AcceptEncodingInvalidHeader | AcceptEncodingNoHeader | None) -> AcceptEncodingNoHeader: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptEncodingHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptEncodingValidHeader | AcceptEncodingNoHeader: ...

@overload
def create_accept_encoding_header(header_value: AcceptEncodingValidHeader | Literal[""]) -> AcceptEncodingValidHeader: ...
@overload
def create_accept_encoding_header(header_value: AcceptEncodingInvalidHeader) -> AcceptEncodingInvalidHeader: ...
@overload
def create_accept_encoding_header(header_value: AcceptEncodingNoHeader | None) -> AcceptEncodingNoHeader: ...
@overload
def create_accept_encoding_header(header_value: str) -> AcceptEncodingValidHeader | AcceptEncodingInvalidHeader: ...
@overload
def create_accept_encoding_header(header_value: _AnyAcceptEncodingHeader | str | None) -> _AnyAcceptEncodingHeader: ...
def accept_encoding_property() -> _AcceptEncodingProperty: ...

class AcceptLanguage:
    @classmethod
    def parse(cls, value: str) -> Iterator[tuple[str, float]]: ...

class AcceptLanguageValidHeader(AcceptLanguage):
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> list[tuple[str, float]]: ...
    def __add__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def __bool__(self) -> Literal[True]: ...
    def __contains__(self, offer: str) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __radd__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self: ...
    def basic_filtering(self, language_tags: Sequence[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    @overload
    def lookup(
        self, language_tags: Sequence[str], default_range: str | None, default_tag: str, default: None = None
    ) -> str | None: ...
    @overload
    def lookup(
        self, language_tags: Sequence[str], *, default_range: str | None = None, default_tag: str, default: None = None
    ) -> str | None: ...
    @overload
    def lookup(
        self, language_tags: Sequence[str], default_range: str | None, default_tag: None, default: _T | Callable[[], _T]
    ) -> _T | str | None: ...
    @overload
    def lookup(
        self, language_tags: Sequence[str], default_range: str | None, default_tag: str, default: _T | Callable[[], _T]
    ) -> _T | str: ...
    @overload
    def lookup(
        self,
        language_tags: Sequence[str],
        *,
        default_range: str | None = None,
        default_tag: None = None,
        default: _T | Callable[[], _T],
    ) -> _T | str | None: ...
    @overload
    def lookup(
        self, language_tags: Sequence[str], *, default_range: str | None = None, default_tag: str, default: _T | Callable[[], _T]
    ) -> _T | str: ...
    def quality(self, offer: str) -> float | None: ...

class _AcceptLanguageInvalidOrNoHeader(AcceptLanguage):
    def __bool__(self) -> Literal[False]: ...
    def __contains__(self, offer: str) -> Literal[True]: ...
    def __iter__(self) -> Iterator[str]: ...
    def basic_filtering(self, language_tags: Iterable[str]) -> list[tuple[str, float]]: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: None = None) -> str | None: ...
    @overload
    def best_match(self, offers: Iterable[str | tuple[str, float] | list[Any]], default_match: str) -> str: ...
    @overload
    def lookup(self, language_tags: object, default_range: object, default_tag: str, default: object = None) -> str: ...
    @overload
    def lookup(
        self, language_tags: object = None, *, default_range: object = None, default_tag: str, default: object = None
    ) -> str: ...
    @overload
    def lookup(self, language_tags: object, default_range: object, default_tag: None, default: _T | Callable[[], _T]) -> _T: ...
    @overload
    def lookup(
        self,
        language_tags: object = None,
        *,
        default_range: object = None,
        default_tag: None = None,
        default: _T | Callable[[], _T],
    ) -> _T: ...
    def quality(self, offer: str) -> float | None: ...

class AcceptLanguageNoHeader(_AcceptLanguageInvalidOrNoHeader):
    def __init__(self) -> None: ...
    def copy(self) -> Self: ...
    @property
    def header_value(self) -> None: ...
    @property
    def parsed(self) -> None: ...
    @overload
    def __add__(self, other: AcceptLanguageValidHeader) -> AcceptLanguageValidHeader: ...
    @overload
    def __add__(self, other: AcceptLanguageInvalidHeader | AcceptLanguageNoHeader | Literal[""] | None) -> Self: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptLanguageValidHeader: ...
    @overload
    def __radd__(self, other: AcceptLanguageValidHeader) -> AcceptLanguageValidHeader: ...
    @overload
    def __radd__(self, other: AcceptLanguageInvalidHeader | AcceptLanguageNoHeader | Literal[""] | None) -> Self: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> Self | AcceptLanguageValidHeader: ...

class AcceptLanguageInvalidHeader(_AcceptLanguageInvalidOrNoHeader):
    def __init__(self, header_value: str) -> None: ...
    def copy(self) -> Self: ...
    @property
    def header_value(self) -> str: ...
    @property
    def parsed(self) -> None: ...
    @overload
    def __add__(self, other: AcceptLanguageValidHeader) -> AcceptLanguageValidHeader: ...
    @overload
    def __add__(
        self, other: AcceptLanguageInvalidHeader | AcceptLanguageNoHeader | Literal[""] | None
    ) -> AcceptLanguageNoHeader: ...
    @overload
    def __add__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptLanguageValidHeader | AcceptLanguageNoHeader: ...
    @overload
    def __radd__(self, other: AcceptLanguageValidHeader) -> AcceptLanguageValidHeader: ...
    @overload
    def __radd__(
        self, other: AcceptLanguageInvalidHeader | AcceptLanguageNoHeader | Literal[""] | None
    ) -> AcceptLanguageNoHeader: ...
    @overload
    def __radd__(
        self,
        other: (
            _AnyAcceptLanguageHeader
            | SupportsItems[str, float]
            | _ListOrTuple[str | tuple[str, float] | list[Any]]
            | _SupportsStr
            | str
            | None
        ),
    ) -> AcceptLanguageValidHeader | AcceptLanguageNoHeader: ...

@overload
def create_accept_language_header(header_value: AcceptLanguageValidHeader | Literal[""]) -> AcceptLanguageValidHeader: ...
@overload
def create_accept_language_header(header_value: AcceptLanguageNoHeader | None) -> AcceptLanguageNoHeader: ...
@overload
def create_accept_language_header(header_value: AcceptLanguageInvalidHeader) -> AcceptLanguageInvalidHeader: ...
@overload
def create_accept_language_header(header_value: str) -> AcceptLanguageValidHeader | AcceptLanguageInvalidHeader: ...
@overload
def create_accept_language_header(header_value: _AnyAcceptLanguageHeader | str | None) -> _AnyAcceptLanguageHeader: ...
def accept_language_property() -> _AcceptLanguageProperty: ...
