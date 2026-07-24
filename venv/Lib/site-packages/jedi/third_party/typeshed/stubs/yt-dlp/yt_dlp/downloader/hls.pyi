from ..extractor.common import _InfoDict
from .fragment import FragmentFD

class HlsFD(FragmentFD):
    FD_NAME: str
    @classmethod
    def can_download(cls, manifest: str, info_dict: _InfoDict, allow_unplayable_formats: bool = False) -> bool: ...
