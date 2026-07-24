from typing import ClassVar

from .common import InfoExtractor

class RtmpIE(InfoExtractor):
    IE_DESC: ClassVar[bool]

class MmsIE(InfoExtractor):
    IE_DESC: ClassVar[bool]

class ViewSourceIE(InfoExtractor):
    IE_DESC: ClassVar[bool]
