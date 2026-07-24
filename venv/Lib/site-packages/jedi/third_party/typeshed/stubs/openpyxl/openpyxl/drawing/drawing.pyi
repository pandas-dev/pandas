from _typeshed import Incomplete

from .spreadsheet_drawing import AbsoluteAnchor, OneCellAnchor

class Drawing:
    count: int
    name: str
    description: str
    coordinates: Incomplete
    left: int
    top: int
    resize_proportional: bool
    rotation: int
    anchortype: str
    anchorcol: int
    anchorrow: int
    def __init__(self) -> None: ...
    @property
    def width(self) -> int: ...
    @width.setter
    def width(self, w: int) -> None: ...
    @property
    def height(self) -> int: ...
    @height.setter
    def height(self, h: int) -> None: ...
    def set_dimension(self, w: int = 0, h: int = 0) -> None: ...
    @property
    def anchor(self) -> AbsoluteAnchor | OneCellAnchor: ...
