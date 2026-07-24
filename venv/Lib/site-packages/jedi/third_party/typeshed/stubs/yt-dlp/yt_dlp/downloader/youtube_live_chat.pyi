from typing import Any

from ..extractor.common import _InfoDict
from .fragment import FragmentFD

class YoutubeLiveChatFD(FragmentFD):
    def real_download(self, filename: str, info_dict: _InfoDict) -> bool: ...
    @staticmethod
    def parse_live_timestamp(action: dict[str, Any]) -> int | None: ...
