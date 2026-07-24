from collections.abc import Collection
from typing import Any

import pandas as pd
from shapely import Geometry
from shapely.geometry.base import BaseGeometry

from ..geoseries import GeoSeries

def collect(
    x: Collection[Geometry] | GeoSeries | pd.Series[Any] | Geometry, multi: bool = False  # Cannot use pd.Series[BaseGeometry]
) -> BaseGeometry: ...
