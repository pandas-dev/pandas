import ssl
from _socket import _Address
from _typeshed import StrOrBytesPath
from collections.abc import Callable, Iterable, Mapping
from socket import AddressFamily, SocketKind
from typing import Any
from typing_extensions import TypeAlias

from ..socks import sockssocket
from ..utils.networking import HTTPHeaderDict
from .common import Request, RequestHandler, Response

def ssl_load_certs(context: ssl.SSLContext, use_certifi: bool = True) -> None: ...
def ssl_load_windows_store_certs(ssl_context: ssl.SSLContext, storename: str) -> None: ...
def make_socks_proxy_opts(socks_proxy: str) -> dict[str, Any]: ...
def get_redirect_method(method: str, status: int) -> str: ...
def make_ssl_context(
    verify: bool = True,
    client_certificate: StrOrBytesPath | None = None,
    client_certificate_key: StrOrBytesPath | None = None,
    client_certificate_password: Callable[[], str | bytes | bytearray] | str | bytes | bytearray | None = None,
    legacy_support: bool = False,
    use_certifi: bool = True,
) -> ssl.SSLContext: ...

class InstanceStoreMixin:
    def __init__(self, **kwargs: Any) -> None: ...  # Passed to non-existent parent so MRO works.

def add_accept_encoding_header(headers: HTTPHeaderDict, supported_encodings: Iterable[str]) -> None: ...
def wrap_request_errors(
    func: Callable[[RequestHandler, Request], Response | None],
) -> Callable[[RequestHandler, Request], None]: ...

_IPAddress: TypeAlias = tuple[
    AddressFamily, SocketKind, int, str, tuple[str, int] | tuple[str, int, int, int] | tuple[int, bytes]
]

def create_socks_proxy_socket(
    dest_addr: _Address, proxy_args: Mapping[str, Any], proxy_ip_addr: _IPAddress, timeout: float | None, source_address: _Address
) -> sockssocket: ...
def create_connection(
    address: tuple[str, int],
    timeout: int = ...,
    source_address: _Address | None = None,
    *,
    _create_socket_func: Callable[[_IPAddress, int, _Address], sockssocket] = ...,
) -> sockssocket: ...
