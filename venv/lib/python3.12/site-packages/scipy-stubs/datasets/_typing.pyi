from collections.abc import Callable
from typing import Literal as L, LiteralString, Protocol, TypeAlias, overload, type_check_only

import numpy as np
import optype.numpy as onp

__all__ = "AscentDataset", "CanFetch", "Dataset", "ECGDataset", "Face2Dataset", "Face3Dataset", "Fetcher"

@type_check_only
class CanFetch(Protocol):
    def fetch(self, dataset_name: LiteralString, /) -> LiteralString: ...

AscentDataset: TypeAlias = onp.Array2D[np.uint8]
ECGDataset: TypeAlias = onp.Array1D[np.float64]
Face2Dataset: TypeAlias = onp.Array2D[np.uint8]
Face3Dataset: TypeAlias = onp.Array3D[np.uint8]
_FaceDataset: TypeAlias = Face2Dataset | Face3Dataset
Dataset: TypeAlias = AscentDataset | ECGDataset | _FaceDataset

_AscentFetcher: TypeAlias = Callable[[], AscentDataset]
_ECGFetcher: TypeAlias = Callable[[], ECGDataset]

@type_check_only
class _FaceFetcher(Protocol):
    @overload
    def __call__(self, /, gray: L[True, 1]) -> Face2Dataset: ...
    @overload
    def __call__(self, /, gray: L[False, 0] = False) -> Face3Dataset: ...

Fetcher: TypeAlias = _AscentFetcher | _ECGFetcher | _FaceFetcher
