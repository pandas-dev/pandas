from _typeshed import SupportsRead
from pathlib import Path
from types import ModuleType
from typing import Any, Literal
from typing_extensions import TypeAlias

from openpyxl.drawing.spreadsheet_drawing import _AnchorBase

# Is actually PIL.Image.Image
_PILImageImage: TypeAlias = Any
# same as first parameter of PIL.Image.open
_PILImageFilePath: TypeAlias = str | bytes | Path | SupportsRead[bytes]

PILImage: ModuleType | Literal[False]

class Image:
    anchor: str | _AnchorBase
    ref: _PILImageImage | _PILImageFilePath
    width: int
    height: int
    format: str
    def __init__(self, img: _PILImageImage | _PILImageFilePath) -> None: ...
    @property
    def path(self) -> str: ...
