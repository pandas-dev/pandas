import abc

from ...._types import Box
from ....main import ActiveWithNeighbors
from ...base import BaseImage

class QRModuleDrawer(abc.ABC, metaclass=abc.ABCMeta):
    needs_neighbors: bool = False
    img: BaseImage
    def initialize(self, img: BaseImage) -> None: ...
    @abc.abstractmethod
    def drawrect(self, box: Box, is_active: bool | ActiveWithNeighbors) -> None: ...
