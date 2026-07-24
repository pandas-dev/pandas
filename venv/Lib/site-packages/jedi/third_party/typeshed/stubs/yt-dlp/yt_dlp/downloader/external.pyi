import enum
import functools
from _typeshed import Unused
from collections.abc import Iterable, Mapping
from typing import Any, Literal

from ..extractor.common import _InfoDict
from .fragment import FragmentFD

class Features(enum.Enum):
    TO_STDOUT = 1
    MULTIPLE_FORMATS = 2

class ExternalFD(FragmentFD):
    SUPPORTED_PROTOCOLS: tuple[str, ...]
    SUPPORTED_FEATURES: tuple[Features, ...]
    @classmethod
    def get_basename(cls) -> str: ...
    @property
    def EXE_NAME(cls) -> str: ...
    @functools.cached_property
    def exe(self) -> str: ...
    @classmethod
    def available(cls, path: str | None = None) -> str | Literal[False]: ...
    @classmethod
    def supports(cls, info_dict: _InfoDict) -> bool: ...
    @classmethod
    def can_download(cls, info_dict: _InfoDict, path: str | None = None) -> bool: ...

class CurlFD(ExternalFD):
    AVAILABLE_OPT: str

class AxelFD(ExternalFD):
    AVAILABLE_OPT: str

class WgetFD(ExternalFD):
    AVAILABLE_OPT: str

class Aria2cFD(ExternalFD):
    AVAILABLE_OPT: str
    SUPPORTED_PROTOCOLS: tuple[str, ...]
    @staticmethod
    def supports_manifest(manifest: str) -> bool: ...
    def aria2c_rpc(self, rpc_port: int, rpc_secret: str, method: str, params: Iterable[str] = ()) -> Any: ...

class HttpieFD(ExternalFD):
    AVAILABLE_OPT: str
    EXE_NAME: str

class FFmpegFD(ExternalFD):
    SUPPORTED_PROTOCOLS: tuple[str, ...]
    SUPPORTED_FEATURES: tuple[Features, ...]
    @classmethod
    def available(cls, path: str | None = None) -> bool: ...  # type: ignore[override]
    def on_process_started(self, proc: Unused, stdin: Unused) -> None: ...
    @classmethod
    def can_merge_formats(cls, info_dict: _InfoDict, params: Mapping[str, Any]) -> bool: ...

def list_external_downloaders() -> list[str]: ...
def get_external_downloader(external_downloader: str) -> type[ExternalFD]: ...
