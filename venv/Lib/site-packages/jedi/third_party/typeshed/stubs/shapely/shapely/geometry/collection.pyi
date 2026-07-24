from collections.abc import Collection
from typing import overload
from typing_extensions import Self

from .._typing import OptGeoArrayLike
from .base import BaseMultipartGeometry, GeometrySequence, _GeoT_co

class GeometryCollection(BaseMultipartGeometry[_GeoT_co]):
    # Overloads of __new__ are used because mypy is unable to narrow the typevar otherwise
    __slots__: list[str] = []
    @overload
    def __new__(
        self, geoms: BaseMultipartGeometry[_GeoT_co] | GeometrySequence[BaseMultipartGeometry[_GeoT_co]] | Collection[_GeoT_co]
    ) -> Self: ...
    @overload
    def __new__(self, geoms: OptGeoArrayLike = None) -> Self: ...
    # more precise base overrides
    @property
    def boundary(self) -> None: ...
