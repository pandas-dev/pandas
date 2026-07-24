from dataclasses import dataclass
from typing import ClassVar

from pandas import DataFrame
from seaborn._core.groupby import GroupBy
from seaborn._core.scales import Scale

@dataclass
class Stat:
    group_by_orient: ClassVar[bool]
    def __call__(self, data: DataFrame, groupby: GroupBy, orient: str, scales: dict[str, Scale]) -> DataFrame: ...
