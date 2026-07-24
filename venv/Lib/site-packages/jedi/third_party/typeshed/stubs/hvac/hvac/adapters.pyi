from abc import ABCMeta, abstractmethod
from collections.abc import Mapping
from typing import Any, Generic, TypeVar, type_check_only
from typing_extensions import Self

from requests import Response, Session

_R = TypeVar("_R")

class Adapter(Generic[_R], metaclass=ABCMeta):
    @classmethod
    def from_adapter(cls, adapter: Adapter[_R]) -> Self: ...
    base_uri: str
    token: str | None
    namespace: str | None
    session: Session
    allow_redirects: bool
    ignore_exceptions: bool
    strict_http: bool
    request_header: bool
    def __init__(
        self,
        base_uri: str = "http://localhost:8200",
        token: str | None = None,
        cert: tuple[str, str] | None = None,
        verify: bool = True,
        timeout: int = 30,
        proxies: Mapping[str, str] | None = None,
        allow_redirects: bool = True,
        session: Session | None = None,
        namespace: str | None = None,
        ignore_exceptions: bool = False,
        strict_http: bool = False,
        request_header: bool = True,
    ) -> None: ...
    @staticmethod
    def urljoin(*args: object) -> str: ...
    def close(self) -> None: ...
    def get(self, url: str, **kwargs: Any) -> _R: ...
    def post(self, url: str, **kwargs: Any) -> _R: ...
    def put(self, url: str, **kwargs: Any) -> _R: ...
    def delete(self, url: str, **kwargs: Any) -> _R: ...
    def list(self, url: str, **kwargs: Any) -> _R: ...
    def head(self, url: str, **kwargs: Any) -> _R: ...
    def login(self, url: str, use_token: bool = True, **kwargs: Any) -> Response: ...
    @abstractmethod
    def get_login_token(self, response: _R) -> str: ...
    @abstractmethod
    def request(
        self, method, url: str, headers: Mapping[str, str] | None = None, raise_exception: bool = True, **kwargs: Any
    ) -> _R: ...

@type_check_only
class _GenericRawAdapter(Adapter[_R]):
    def get_login_token(self, response: _R) -> str: ...
    def request(
        self, method: str, url: str, headers: Mapping[str, str] | None = None, raise_exception: bool = True, **kwargs: Any
    ) -> _R: ...

class RawAdapter(_GenericRawAdapter[Response]): ...

class JSONAdapter(_GenericRawAdapter[Response | dict[Any, Any]]):
    def get_login_token(self, response: Response | dict[Any, Any]) -> str: ...
    def request(self, *args: Any, **kwargs: Any) -> Response | dict[Any, Any]: ...

Request = RawAdapter
