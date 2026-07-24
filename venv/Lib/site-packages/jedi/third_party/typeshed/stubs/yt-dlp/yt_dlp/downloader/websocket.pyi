from _typeshed import SupportsWrite

from ..extractor.common import _InfoDict
from .common import FileDownloader

class FFmpegSinkFD(FileDownloader):
    def real_download(self, filename: str, info_dict: _InfoDict) -> bool: ...
    async def real_connection(self, sink: SupportsWrite[bytes], info_dict: _InfoDict) -> None: ...

class WebSocketFragmentFD(FFmpegSinkFD):
    async def real_connection(self, sink: SupportsWrite[bytes], info_dict: _InfoDict) -> None: ...
