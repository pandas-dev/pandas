from collections.abc import Collection, Sequence
from pathlib import Path
from typing import Generic, Literal, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

import numpy as np
import numpy.typing as npt

from . import _EncodedRLE

PYTHON_VERSION: int

@type_check_only
class _Image(TypedDict):
    id: int
    width: int
    height: int
    file_name: str

_TPolygonSegmentation: TypeAlias = list[list[float]]

@type_check_only
class _RLE(TypedDict):
    size: list[int]
    counts: list[int]

@type_check_only
class _Annotation(TypedDict):
    id: int
    image_id: int
    category_id: int
    segmentation: _TPolygonSegmentation | _RLE | _EncodedRLE
    area: float
    bbox: list[float]
    iscrowd: int

_TSeg = TypeVar("_TSeg", _TPolygonSegmentation, _RLE, _EncodedRLE)

@type_check_only
class _AnnotationG(TypedDict, Generic[_TSeg]):
    id: int
    image_id: int
    category_id: int
    segmentation: _TSeg
    area: float
    bbox: list[float]
    iscrowd: int

@type_check_only
class _Category(TypedDict):
    id: int
    name: str
    supercategory: str

@type_check_only
class _Dataset(TypedDict):
    images: list[_Image]
    annotations: list[_Annotation]
    categories: list[_Category]

class COCO:
    anns: dict[int, _Annotation]
    dataset: _Dataset
    cats: dict[int, _Category]
    imgs: dict[int, _Image]
    imgToAnns: dict[int, list[_Annotation]]
    catToImgs: dict[int, list[int]]
    def __init__(self, annotation_file: str | Path | None = None) -> None: ...
    def createIndex(self) -> None: ...
    def info(self) -> None: ...
    def getAnnIds(
        self,
        imgIds: Collection[int] | int = [],
        catIds: Collection[int] | int = [],
        areaRng: Sequence[float] = [],
        iscrowd: bool | None = None,
    ) -> list[int]: ...
    def getCatIds(
        self, catNms: Collection[str] | str = [], supNms: Collection[str] | str = [], catIds: Collection[int] | int = []
    ) -> list[int]: ...
    def getImgIds(self, imgIds: Collection[int] | int = [], catIds: list[int] | int = []) -> list[int]: ...
    def loadAnns(self, ids: Collection[int] | int = []) -> list[_Annotation]: ...
    def loadCats(self, ids: Collection[int] | int = []) -> list[_Category]: ...
    def loadImgs(self, ids: Collection[int] | int = []) -> list[_Image]: ...
    def showAnns(self, anns: Sequence[_Annotation], draw_bbox: bool = False) -> None: ...
    def loadRes(self, resFile: str) -> COCO: ...
    def download(self, tarDir: str | None = None, imgIds: Collection[int] = []) -> Literal[-1] | None: ...
    def loadNumpyAnnotations(self, data: npt.NDArray[np.float64]) -> list[_Annotation]: ...
    @overload
    def annToRLE(self, ann: _AnnotationG[_RLE]) -> _RLE: ...
    @overload
    def annToRLE(self, ann: _AnnotationG[_EncodedRLE]) -> _EncodedRLE: ...
    @overload
    def annToRLE(self, ann: _AnnotationG[_TPolygonSegmentation]) -> _EncodedRLE: ...
    def annToMask(self, ann: _Annotation) -> npt.NDArray[np.uint8]: ...
