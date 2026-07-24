from typing import Final

from ._config import options as options
from ._exports import (
    gpd as gpd,
    list_layers as list_layers,
    np as np,
    pd as pd,
    read_feather as read_feather,
    read_file as read_file,
    read_parquet as read_parquet,
    read_postgis as read_postgis,
)
from .array import points_from_xy as points_from_xy
from .geodataframe import GeoDataFrame as GeoDataFrame
from .geoseries import GeoSeries as GeoSeries
from .tools import clip as clip, overlay as overlay, sjoin as sjoin, sjoin_nearest as sjoin_nearest
from .tools._show_versions import show_versions as show_versions

__version__: Final[str]
