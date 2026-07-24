from _typeshed import Incomplete
from typing import Final, Literal

__version__: Final[str]

class PDFImage:
    image: Incomplete
    x: Incomplete
    y: Incomplete
    width: Incomplete
    height: Incomplete
    filename: Incomplete
    imageCaching: bool | Literal[0, 1]
    colorSpace: str
    bitsPerComponent: int
    filters: Incomplete
    source: Incomplete
    def __init__(self, image, x, y, width=None, height=None, caching: bool | Literal[0, 1] = 0) -> None: ...
    def jpg_imagedata(self) -> tuple[list[str], Incomplete, Incomplete]: ...
    def cache_imagedata(self) -> list[str]: ...
    def PIL_imagedata(self) -> tuple[list[str], Incomplete, Incomplete]: ...
    def non_jpg_imagedata(self, image) -> tuple[list[str], int, int]: ...
    imageData: Incomplete
    imgwidth: Incomplete
    imgheight: Incomplete
    def getImageData(self, preserveAspectRatio: bool = False) -> None: ...
    def drawInlineImage(
        self,
        canvas,
        preserveAspectRatio: bool = False,
        anchor: str = "sw",
        anchorAtXY: bool = False,
        showBoundary: bool = False,
        extraReturn=None,
    ) -> bool: ...
    def format(self, document) -> bytes: ...
