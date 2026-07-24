from collections.abc import Sequence
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

@type_check_only
class _ProxyConfigDict(TypedDict):
    http: NotRequired[str]
    https: NotRequired[str]
    ftpProxy: NotRequired[str]
    noProxy: NotRequired[str]

@type_check_only
class _Environment(TypedDict):
    http_proxy: NotRequired[str]
    HTTP_PROXY: NotRequired[str]
    https_proxy: NotRequired[str]
    HTTPS_PROXY: NotRequired[str]
    ftp_proxy: NotRequired[str]
    FTP_PROXY: NotRequired[str]
    no_proxy: NotRequired[str]
    NO_PROXY: NotRequired[str]

class ProxyConfig(dict[str, str]):
    @property
    def http(self) -> str | None: ...
    @property
    def https(self) -> str | None: ...
    @property
    def ftp(self) -> str | None: ...
    @property
    def no_proxy(self) -> str | None: ...
    @staticmethod
    def from_dict(config: _ProxyConfigDict) -> ProxyConfig: ...
    def get_environment(self) -> _Environment: ...
    def inject_proxy_environment(self, environment: None | Sequence[str]) -> None | Sequence[str]: ...
