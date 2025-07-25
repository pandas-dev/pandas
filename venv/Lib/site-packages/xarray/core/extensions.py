from __future__ import annotations

import warnings

from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.datatree import DataTree


class AccessorRegistrationWarning(Warning):
    """Warning for conflicts in accessor registration."""


class _CachedAccessor:
    """Custom property-like object (descriptor) for caching accessors."""

    def __init__(self, name, accessor):
        self._name = name
        self._accessor = accessor

    def __get__(self, obj, cls):
        if obj is None:
            # we're accessing the attribute of the class, i.e., Dataset.geo
            return self._accessor

        # Use the same dict as @pandas.util.cache_readonly.
        # It must be explicitly declared in obj.__slots__.
        try:
            cache = obj._cache
        except AttributeError:
            cache = obj._cache = {}

        try:
            return cache[self._name]
        except KeyError:
            pass

        try:
            accessor_obj = self._accessor(obj)
        except AttributeError as err:
            # __getattr__ on data object will swallow any AttributeErrors
            # raised when initializing the accessor, so we need to raise as
            # something else (GH933):
            raise RuntimeError(f"error initializing {self._name!r} accessor.") from err

        cache[self._name] = accessor_obj
        return accessor_obj


def _register_accessor(name, cls):
    def decorator(accessor):
        if hasattr(cls, name):
            warnings.warn(
                f"registration of accessor {accessor!r} under name {name!r} for type {cls!r} is "
                "overriding a preexisting attribute with the same name.",
                AccessorRegistrationWarning,
                stacklevel=2,
            )
        setattr(cls, name, _CachedAccessor(name, accessor))
        return accessor

    return decorator


def register_dataarray_accessor(name):
    """Register a custom accessor on xarray.DataArray objects.

    Parameters
    ----------
    name : str
        Name under which the accessor should be registered. A warning is issued
        if this name conflicts with a preexisting attribute.

    See Also
    --------
    register_dataset_accessor
    """
    return _register_accessor(name, DataArray)


def register_dataset_accessor(name):
    """Register a custom property on xarray.Dataset objects.

    Parameters
    ----------
    name : str
        Name under which the accessor should be registered. A warning is issued
        if this name conflicts with a preexisting attribute.

    Examples
    --------
    In your library code:

    >>> @xr.register_dataset_accessor("geo")
    ... class GeoAccessor:
    ...     def __init__(self, xarray_obj):
    ...         self._obj = xarray_obj
    ...
    ...     @property
    ...     def center(self):
    ...         # return the geographic center point of this dataset
    ...         lon = self._obj.latitude
    ...         lat = self._obj.longitude
    ...         return (float(lon.mean()), float(lat.mean()))
    ...
    ...     def plot(self):
    ...         # plot this array's data on a map, e.g., using Cartopy
    ...         pass
    ...

    Back in an interactive IPython session:

    >>> ds = xr.Dataset(
    ...     {"longitude": np.linspace(0, 10), "latitude": np.linspace(0, 20)}
    ... )
    >>> ds.geo.center
    (10.0, 5.0)
    >>> ds.geo.plot()  # plots data on a map

    See Also
    --------
    register_dataarray_accessor
    """
    return _register_accessor(name, Dataset)


def register_datatree_accessor(name):
    """Register a custom accessor on DataTree objects.

    Parameters
    ----------
    name : str
        Name under which the accessor should be registered. A warning is issued
        if this name conflicts with a preexisting attribute.

    See Also
    --------
    xarray.register_dataarray_accessor
    xarray.register_dataset_accessor
    """
    return _register_accessor(name, DataTree)
