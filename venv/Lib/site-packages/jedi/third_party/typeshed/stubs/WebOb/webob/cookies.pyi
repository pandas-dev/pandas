from _typeshed import sentinel
from _typeshed.wsgi import WSGIEnvironment
from collections.abc import Collection, ItemsView, Iterator, KeysView, MutableMapping, ValuesView
from datetime import date, datetime, timedelta
from time import _TimeTuple, struct_time
from typing import Any, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from webob._types import AsymmetricProperty
from webob.request import BaseRequest
from webob.response import Response

__all__ = [
    "Cookie",
    "CookieProfile",
    "SignedCookieProfile",
    "SignedSerializer",
    "JSONSerializer",
    "Base64Serializer",
    "make_cookie",
]

_T = TypeVar("_T")
# we accept both the official spelling and the one used in the WebOb docs
# the implementation compares after lower() so technically there are more
# valid spellings, but it seems more natural to support these two spellings
_SameSitePolicy: TypeAlias = Literal["Strict", "Lax", "None", "strict", "lax", "none"]

@type_check_only
class _Serializer(Protocol):
    def dumps(self, appstruct: Any, /) -> bytes: ...
    def loads(self, bstruct: bytes, /) -> Any: ...

class RequestCookies(MutableMapping[str, str]):
    def __init__(self, environ: WSGIEnvironment) -> None: ...
    def __setitem__(self, name: str, value: str) -> None: ...
    def __getitem__(self, name: str) -> str: ...
    @overload
    def get(self, name: str, default: None = None) -> str | None: ...
    @overload
    def get(self, name: str, default: str) -> str: ...
    @overload
    def get(self, name: str, default: _T) -> str | _T: ...
    def __delitem__(self, name: str) -> None: ...
    def keys(self) -> KeysView[str]: ...
    def values(self) -> ValuesView[str]: ...
    def items(self) -> ItemsView[str, str]: ...
    def __contains__(self, name: object) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...
    def clear(self) -> None: ...

class Cookie(dict[bytes, Morsel]):
    def __init__(self, input: str | None = None) -> None: ...
    def load(self, data: str) -> None: ...
    def add(self, key: str | bytes, val: str | bytes) -> Morsel | dict[bytes, bytes]: ...
    def __setitem__(self, key: str | bytes, val: str | bytes) -> Morsel | dict[bytes, bytes]: ...  # type: ignore[override]
    def serialize(self, full: bool = True) -> str: ...
    def values(self) -> list[Morsel]: ...  # type: ignore[override]
    def __str__(self, full: bool = True) -> str: ...

class Morsel(dict[bytes, bytes | bool | None]):
    __slots__ = ("name", "value")
    name: bytes
    value: bytes
    def __init__(self, name: str | bytes, value: str | bytes) -> None: ...
    @property
    def path(self) -> bytes | None: ...
    @path.setter
    def path(self, v: bytes | None) -> None: ...
    @property
    def domain(self) -> bytes | None: ...
    @domain.setter
    def domain(self, v: bytes | None) -> None: ...
    @property
    def comment(self) -> bytes | None: ...
    @comment.setter
    def comment(self, v: bytes | None) -> None: ...
    expires: AsymmetricProperty[bytes | None, datetime | date | timedelta | _TimeTuple | struct_time | int | str | bytes | None]
    max_age: AsymmetricProperty[bytes | None, timedelta | int | str | bytes | None]
    httponly: AsymmetricProperty[bool, bool | None]
    secure: AsymmetricProperty[bool, bool | None]
    samesite: AsymmetricProperty[bytes, _SameSitePolicy | bytes]
    def __setitem__(self, k: str | bytes, v: bytes | bool | None) -> None: ...
    def serialize(self, full: bool = True) -> str: ...
    def __str__(self, full: bool = True) -> str: ...

def make_cookie(
    name: str | bytes,
    value: str | bytes | None,
    max_age: int | timedelta | None = None,
    path: str = "/",
    domain: str | None = None,
    secure: bool | None = False,
    httponly: bool | None = False,
    comment: str | None = None,
    samesite: _SameSitePolicy | None = None,
) -> str: ...

class JSONSerializer:
    def dumps(self, appstruct: Any) -> bytes: ...
    def loads(self, bstruct: bytes | str) -> Any: ...

class Base64Serializer:
    serializer: _Serializer
    def __init__(self, serializer: _Serializer | None = None) -> None: ...
    def dumps(self, appstruct: Any) -> bytes: ...
    def loads(self, bstruct: bytes) -> Any: ...

class SignedSerializer:
    salt: str | bytes
    secret: str | bytes
    hashalg: str
    salted_secret: bytes
    digest_size: int
    serializer: _Serializer
    def __init__(
        self, secret: str | bytes, salt: str | bytes, hashalg: str = "sha512", serializer: _Serializer | None = None
    ) -> None: ...
    def dumps(self, appstruct: Any) -> bytes: ...
    def loads(self, bstruct: bytes) -> Any: ...

class CookieProfile:
    cookie_name: str
    secure: bool
    max_age: int | timedelta | None
    httponly: bool | None
    samesite: _SameSitePolicy | None
    path: str
    domains: Collection[str] | None
    serializer: _Serializer
    request: BaseRequest | None
    def __init__(
        self,
        cookie_name: str,
        secure: bool = False,
        max_age: int | timedelta | None = None,
        httponly: bool | None = None,
        samesite: _SameSitePolicy | None = None,
        path: str = "/",
        # even though the docs claim any iterable is fine, that is
        # clearly not the case judging by the implementation
        domains: Collection[str] | None = None,
        serializer: _Serializer | None = None,
    ) -> None: ...
    def __call__(self, request: BaseRequest) -> CookieProfile: ...
    def bind(self, request: BaseRequest) -> CookieProfile: ...
    def get_value(self) -> Any | None: ...
    def set_cookies(
        self,
        response: Response,
        value: Any,
        domains: Collection[str] = sentinel,
        max_age: int | timedelta | None = sentinel,
        path: str = sentinel,
        secure: bool = sentinel,
        httponly: bool = sentinel,
        samesite: _SameSitePolicy | None = sentinel,
    ) -> Response: ...
    def get_headers(
        self,
        value: Any,
        domains: Collection[str] = sentinel,
        max_age: int | timedelta | None = sentinel,
        path: str = sentinel,
        secure: bool = sentinel,
        httponly: bool = sentinel,
        samesite: _SameSitePolicy | None = sentinel,
    ) -> list[tuple[str, str]]: ...

class SignedCookieProfile(CookieProfile):
    secret: str | bytes
    salt: str | bytes
    hashalg: str
    original_serializer: _Serializer
    def __init__(
        self,
        secret: str,
        salt: str,
        cookie_name: str,
        secure: bool = False,
        max_age: int | timedelta | None = None,
        httponly: bool | None = False,
        samesite: _SameSitePolicy | None = None,
        path: str = "/",
        domains: Collection[str] | None = None,
        hashalg: str = "sha512",
        serializer: _Serializer | None = None,
    ) -> None: ...
    def __call__(self, request: BaseRequest) -> SignedCookieProfile: ...
    def bind(self, request: BaseRequest) -> SignedCookieProfile: ...
