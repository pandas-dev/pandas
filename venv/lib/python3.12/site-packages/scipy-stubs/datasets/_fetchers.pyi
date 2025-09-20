from typing import Final, Literal as L, LiteralString, overload

from ._typing import AscentDataset, CanFetch, ECGDataset, Face2Dataset, Face3Dataset

###

data_fetcher: Final[CanFetch | None] = ...  # undocumented

def fetch_data(
    dataset_name: L["ascent.dat", "ecg.dat", "face.dat"], data_fetcher: CanFetch | None = None
) -> LiteralString: ...  # undocumented

#
def ascent() -> AscentDataset: ...
def electrocardiogram() -> ECGDataset: ...
@overload
def face(gray: L[True, 1]) -> Face2Dataset: ...
@overload
def face(gray: L[False, 0] = False) -> Face3Dataset: ...
