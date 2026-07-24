import http.client
import ssl
import sys
import urllib.request
from _typeshed import Incomplete, SupportsKeysAndGetItem
from typing import Any, TypeVar

import socks

_K = TypeVar("_K")
_V = TypeVar("_V")

def merge_dict(a: dict[_K, _V], b: SupportsKeysAndGetItem[_K, _V]) -> dict[_K, _V]: ...  # undocumented
def is_ip(s: str) -> bool: ...  # undocumented

socks4_no_rdns: set[str]  # undocumented

class SocksiPyConnection(http.client.HTTPConnection):  # undocumented
    proxyargs: tuple[int, str, int | None, bool, str | None, str | None]
    sock: socks.socksocket
    def __init__(
        self,
        proxytype: int,
        proxyaddr: str,
        proxyport: int | None = None,
        rdns: bool = True,
        username: str | None = None,
        password: str | None = None,
        host: str | None = None,
        port: int | None = None,
        timeout: float | None = ...,
        source_address: tuple[str, int] | None = None,
        blocksize: int = 8192,
    ) -> None: ...
    def connect(self) -> None: ...

class SocksiPyConnectionS(http.client.HTTPSConnection):  # undocumented
    proxyargs: tuple[int, str, int | None, bool, str | None, str | None]
    sock: socks.socksocket
    if sys.version_info >= (3, 12):
        def __init__(
            self,
            proxytype: int,
            proxyaddr: str,
            proxyport: int | None = None,
            rdns: bool = True,
            username: str | None = None,
            password: str | None = None,
            host: str | None = None,
            port: int | None = None,
            *,
            timeout: float | None = ...,
            source_address: tuple[str, int] | None = None,
            context: ssl.SSLContext | None = None,
            blocksize: int = 8192,
        ) -> None: ...
    else:
        def __init__(
            self,
            proxytype: int,
            proxyaddr: str,
            proxyport: int | None = None,
            rdns: bool = True,
            username: str | None = None,
            password: str | None = None,
            host: str | None = None,
            port: int | None = None,
            key_file: str | None = None,
            cert_file: str | None = None,
            timeout: float | None = ...,
            source_address: tuple[str, int] | None = None,
            *,
            context: ssl.SSLContext | None = None,
            check_hostname: bool | None = None,
            blocksize: int = 8192,
        ) -> None: ...

    def connect(self) -> None: ...

class SocksiPyHandler(urllib.request.HTTPHandler, urllib.request.HTTPSHandler):
    args: tuple[Incomplete, ...]  # undocumented
    kw: dict[str, Incomplete]  # undocumented
    def __init__(
        self,
        proxytype: int,
        proxyaddr: str,
        proxyport: int | None = None,
        rdns: bool = True,
        username: str | None = None,
        password: str | None = None,
        *,
        source_address: tuple[str, int] | None = None,
        blocksize: int = 8192,
        **kwargs: Any,  # any additional arguments to `SocksiPyConnection` or `SocksiPyConnectionS`
    ) -> None: ...
    def http_open(self, req: urllib.request.Request) -> http.client.HTTPResponse: ...  # undocumented
    def https_open(self, req: urllib.request.Request) -> http.client.HTTPResponse: ...  # undocumented
