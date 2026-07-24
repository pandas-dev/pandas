from typing import Final

from . import affinity as affinity
from ._coverage import *
from ._geometry import *
from .constructive import *
from .coordinates import *
from .creation import *
from .errors import setup_signal_checks as setup_signal_checks
from .geometry import (
    GeometryCollection as GeometryCollection,
    LinearRing as LinearRing,
    LineString as LineString,
    MultiLineString as MultiLineString,
    MultiPoint as MultiPoint,
    MultiPolygon as MultiPolygon,
    Point as Point,
    Polygon as Polygon,
)
from .io import *
from .lib import (
    Geometry as Geometry,
    GEOSException as GEOSException,
    geos_capi_version as geos_capi_version,
    geos_capi_version_string as geos_capi_version_string,
    geos_version as geos_version,
    geos_version_string as geos_version_string,
)
from .linear import *
from .measurement import *
from .predicates import *
from .set_operations import *
from .strtree import *

__version__: Final[str]
