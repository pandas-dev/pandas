from collections.abc import Callable, Iterable
from typing import Protocol, type_check_only

from ..base import _ConvertibleToGeoSeries
from ..geodataframe import GeoDataFrame

@type_check_only
class _GeoCoder(Protocol):
    # Represents a geopy.geocoders.base.GeoCoder subclass without actually depending on geopy
    def geocode(self, query: str, /): ...
    def reverse(self, coords, /, exactly_one: bool = ...): ...

# TODO: Use something like `provider: Callable[P, _GeoCoder], **kwargs: P.kwargs` in the functions
# below if this ever becomes a thing
def geocode(strings: Iterable[str], provider: str | Callable[..., _GeoCoder] | None = None, **kwargs) -> GeoDataFrame: ...
def reverse_geocode(
    points: _ConvertibleToGeoSeries, provider: str | Callable[..., _GeoCoder] | None = None, **kwargs
) -> GeoDataFrame: ...
