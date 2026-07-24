from dataclasses import dataclass

from seaborn._stats.base import Stat

@dataclass
class Perc(Stat):
    k: int | list[float] = 5
    method: str = "linear"
