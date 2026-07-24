from ..extractor.common import _InfoDict
from .common import FileDownloader

class RtspFD(FileDownloader):
    def real_download(self, filename: str, info_dict: _InfoDict) -> bool: ...
