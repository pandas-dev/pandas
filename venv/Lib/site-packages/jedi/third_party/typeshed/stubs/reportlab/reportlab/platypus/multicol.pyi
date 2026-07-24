from collections.abc import Sequence

from .flowables import Flowable, _Container, _FindSplitterMixin

class MultiCol(_Container, _FindSplitterMixin, Flowable):
    contents: Sequence[Flowable]
    widths: Sequence[float | str]
    minHeightNeeded: float
    def __init__(
        self,
        contents: Sequence[Flowable],
        widths: Sequence[float | str],
        minHeightNeeded: float = 36,
        spaceBefore: float | None = None,
        spaceAfter: float | None = None,
    ) -> None: ...
    def nWidths(self, aW: float) -> list[float]: ...

__all__ = ["MultiCol"]
