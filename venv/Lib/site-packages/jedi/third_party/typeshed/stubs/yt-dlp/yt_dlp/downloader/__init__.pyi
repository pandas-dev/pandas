from collections.abc import Mapping
from typing import Any, Literal
from typing_extensions import TypeAlias

from ..extractor.common import _InfoDict
from ..utils._utils import NO_DEFAULT
from .common import FileDownloader

__all__ = ["FileDownloader", "get_suitable_downloader", "shorten_protocol_name"]

_Protocol: TypeAlias = Literal[
    "dash_frag_urls",
    "ftps",
    "http_dash_segments",
    "http_dash_segments_generator",
    "https",
    "m3u8",
    "m3u8_frag_urls",
    "m3u8_native",
    "rtmp_ffmpeg",
    "websocket_frag",
]

def get_suitable_downloader(
    info_dict: _InfoDict,
    params: Mapping[str, Any] = {},
    default: FileDownloader | type[NO_DEFAULT] = ...,
    protocol: _Protocol | None = None,
    to_stdout: bool = False,
) -> FileDownloader: ...
def shorten_protocol_name(proto: _Protocol, simplify: bool = False) -> str: ...
