import threading
from collections.abc import Mapping

from ..extractor.common import _InfoDict
from ..utils.networking import HTTPHeaderDict
from .common import FileDownloader

class BunnyCdnFD(FileDownloader):
    def real_download(self, filename: str, info_dict: _InfoDict) -> bool | None: ...
    def ping_thread(
        self,
        stop_event: threading.Event,
        url: str,
        headers: HTTPHeaderDict | Mapping[str, str] | None,
        secret: str,
        context_id: str,
    ) -> None: ...
