from dataclasses import dataclass

from seaborn._stats.base import Stat

@dataclass
class PolyFit(Stat):
    order: int = 2
    gridsize: int = 100

@dataclass
class OLSFit(Stat): ...
