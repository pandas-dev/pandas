from typing import Literal

from ..extractor.common import _InfoDict
from .fragment import FragmentFD

class DashSegmentsFD(FragmentFD):
    FD_NAME: Literal["dashsegments"]
    def real_download(self, filename: str, info_dict: _InfoDict) -> bool: ...
