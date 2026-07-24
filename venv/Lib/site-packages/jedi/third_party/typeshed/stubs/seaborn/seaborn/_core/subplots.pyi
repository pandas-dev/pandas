from collections.abc import Iterator
from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from seaborn._core.plot import FacetSpec, PairSpec

class Subplots:
    subplot_spec: dict[str, Any]
    def __init__(self, subplot_spec: dict[str, Any], facet_spec: FacetSpec, pair_spec: PairSpec) -> None: ...
    def init_figure(
        self,
        pair_spec: PairSpec,
        pyplot: bool = False,
        figure_kws: dict[str, Any] | None = None,
        target: Axes | Figure | SubFigure | None = None,
    ) -> Figure: ...
    def __iter__(self) -> Iterator[dict[str, Any]]: ...
    def __len__(self) -> int: ...
