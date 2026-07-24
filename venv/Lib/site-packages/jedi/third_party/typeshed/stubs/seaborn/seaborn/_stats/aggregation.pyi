from collections.abc import Callable
from dataclasses import dataclass

from seaborn._core.typing import Vector
from seaborn._stats.base import Stat

@dataclass
class Agg(Stat):
    func: str | Callable[[Vector], float] = "mean"

@dataclass
class Est(Stat):
    func: str | Callable[[Vector], float] = "mean"
    errorbar: str | tuple[str, float] = ("ci", 95)
    n_boot: int = 1000
    seed: int | None = None

@dataclass
class Rolling(Stat): ...
