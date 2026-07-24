from collections.abc import Callable
from typing import Literal as L, LiteralString, Protocol, overload, type_check_only

import numpy as np
import optype.numpy as onp

__all__ = "AscentDataset", "CanFetch", "Dataset", "ECGDataset", "Face2Dataset", "Face3Dataset", "Fetcher"

@type_check_only
class CanFetch(Protocol):
    def fetch(self, dataset_name: LiteralString, /) -> LiteralString: ...

type AscentDataset = onp.Array2D[np.uint8]
type ECGDataset = onp.Array1D[np.float64]
type Face2Dataset = onp.Array2D[np.uint8]
type Face3Dataset = onp.Array3D[np.uint8]
type _FaceDataset = Face2Dataset | Face3Dataset
type Dataset = AscentDataset | ECGDataset | _FaceDataset

type _AscentFetcher = Callable[[], AscentDataset]
type _ECGFetcher = Callable[[], ECGDataset]
type Fetcher = _AscentFetcher | _ECGFetcher | _FaceFetcher

@type_check_only
class _FaceFetcher(Protocol):
    @overload
    def __call__(self, /, gray: L[True]) -> Face2Dataset: ...
    @overload
    def __call__(self, /, gray: L[False] = False) -> Face3Dataset: ...
