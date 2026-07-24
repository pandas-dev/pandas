from dataclasses import dataclass

from numpy.typing import ArrayLike
from seaborn._stats.base import Stat

@dataclass
class Count(Stat): ...

@dataclass
class Hist(Stat):
    stat: str = "count"
    bins: str | int | ArrayLike = "auto"
    binwidth: float | None = None
    binrange: tuple[float, float] | None = None
    common_norm: bool | list[str] = True
    common_bins: bool | list[str] = True
    cumulative: bool = False
    discrete: bool = False
    def __post_init__(self) -> None: ...
