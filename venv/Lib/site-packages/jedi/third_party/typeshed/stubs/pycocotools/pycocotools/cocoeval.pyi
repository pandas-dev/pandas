from typing import Literal, TypedDict, type_check_only
from typing_extensions import TypeAlias

import numpy as np
import numpy.typing as npt

from .coco import COCO

_NDFloatArray: TypeAlias = npt.NDArray[np.float64]
_TIOU: TypeAlias = Literal["segm", "bbox", "keypoints"]

@type_check_only
class _ImageEvaluationResult(TypedDict):
    image_id: int
    category_id: int
    aRng: list[int]
    maxDet: int
    dtIds: list[int]
    gtIds: list[int]
    dtMatches: _NDFloatArray
    gtMatches: _NDFloatArray
    dtScores: list[float]
    gtIgnore: _NDFloatArray
    dtIgnore: _NDFloatArray

@type_check_only
class _EvaluationResult(TypedDict):
    params: Params
    counts: list[int]
    date: str
    precision: _NDFloatArray
    recall: _NDFloatArray
    scores: _NDFloatArray

class COCOeval:
    cocoGt: COCO
    cocoDt: COCO
    evalImgs: list[_ImageEvaluationResult]
    eval: _EvaluationResult
    params: Params
    stats: _NDFloatArray
    ious: dict[tuple[int, int], list[float]]
    def __init__(self, cocoGt: COCO | None = None, cocoDt: COCO | None = None, iouType: _TIOU = "segm") -> None: ...
    def evaluate(self) -> None: ...
    def computeIoU(self, imgId: int, catId: int) -> list[float]: ...
    def computeOks(self, imgId: int, catId: int) -> _NDFloatArray: ...
    def evaluateImg(self, imgId: int, catId: int, aRng: list[int], maxDet: int) -> _ImageEvaluationResult: ...
    def accumulate(self, p: Params | None = None) -> None: ...
    def summarize(self) -> None: ...

class Params:
    imgIds: list[int]
    catIds: list[int]
    iouThrs: _NDFloatArray
    recThrs: _NDFloatArray
    maxDets: list[int]
    areaRng: list[list[float]]
    areaRngLbl: list[str]
    useCats: int
    kpt_oks_sigmas: _NDFloatArray
    iouType: _TIOU
    useSegm: int | None
    def __init__(self, iouType: _TIOU = "segm") -> None: ...
    def setDetParams(self) -> None: ...
    def setKpParams(self) -> None: ...
