from _typeshed.wsgi import WSGIApplication
from collections.abc import Mapping, Sequence
from logging import Logger
from typing import Final, NamedTuple, TypeVar

_T = TypeVar("_T")

PROXY_HEADERS: Final[frozenset[str]]

class Forwarded(NamedTuple):
    by: str
    for_: str
    host: str
    proto: str

class MalformedProxyHeader(Exception):
    header: str
    reason: str
    value: str
    def __init__(self, header: str, reason: str, value: str) -> None: ...

def proxy_headers_middleware(
    app: WSGIApplication,
    trusted_proxy: str | None = None,
    trusted_proxy_count: int = 1,
    trusted_proxy_headers: set[str] | None = None,
    clear_untrusted: bool = True,
    log_untrusted: bool = False,
    logger: Logger = ...,
) -> WSGIApplication: ...
def parse_proxy_headers(
    environ: Mapping[str, str], trusted_proxy_count: int, trusted_proxy_headers: set[str], logger: Logger = ...
) -> set[str]: ...
def strip_brackets(addr: Sequence[_T]) -> Sequence[_T]: ...
def clear_untrusted_headers(
    environ: Mapping[str, str], untrusted_headers: Sequence[str], log_warning: bool = False, logger: Logger = ...
) -> None: ...
