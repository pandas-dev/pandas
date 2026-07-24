from typing import ClassVar

from .common import InfoExtractor

class CommonMistakesIE(InfoExtractor):
    IE_DESC: ClassVar[bool]

class UnicodeBOMIE(InfoExtractor):
    IE_DESC: ClassVar[bool]

class BlobIE(InfoExtractor):
    IE_DESC: ClassVar[bool]
