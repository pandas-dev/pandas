from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

import numpy as np

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    AbstractDataStore,
    BackendArray,
    BackendEntrypoint,
    robust_getitem,
)
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core import indexing
from xarray.core.utils import (
    Frozen,
    FrozenDict,
    close_on_error,
    is_dict_like,
    is_remote_uri,
)
from xarray.core.variable import Variable
from xarray.namedarray.pycompat import integer_types

if TYPE_CHECKING:
    import os
    from io import BufferedIOBase

    from xarray.core.dataset import Dataset


class PydapArrayWrapper(BackendArray):
    def __init__(self, array):
        self.array = array

    @property
    def shape(self) -> tuple[int, ...]:
        return self.array.shape

    @property
    def dtype(self):
        return self.array.dtype

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.BASIC, self._getitem
        )

    def _getitem(self, key):
        # pull the data from the array attribute if possible, to avoid
        # downloading coordinate data twice
        array = getattr(self.array, "array", self.array)
        result = robust_getitem(array, key, catch=ValueError)
        result = np.asarray(result)
        # in some cases, pydap doesn't squeeze axes automatically like numpy
        axis = tuple(n for n, k in enumerate(key) if isinstance(k, integer_types))
        if result.ndim + len(axis) != array.ndim and axis:
            result = np.squeeze(result, axis)

        return result


def _fix_attributes(attributes):
    attributes = dict(attributes)
    for k in list(attributes):
        if k.lower() == "global" or k.lower().endswith("_global"):
            # move global attributes to the top level, like the netcdf-C
            # DAP client
            attributes.update(attributes.pop(k))
        elif is_dict_like(attributes[k]):
            # Make Hierarchical attributes to a single level with a
            # dot-separated key
            attributes.update(
                {
                    f"{k}.{k_child}": v_child
                    for k_child, v_child in attributes.pop(k).items()
                }
            )
    return attributes


class PydapDataStore(AbstractDataStore):
    """Store for accessing OpenDAP datasets with pydap.

    This store provides an alternative way to access OpenDAP datasets that may
    be useful if the netCDF4 library is not available.
    """

    def __init__(self, ds):
        """
        Parameters
        ----------
        ds : pydap DatasetType
        """
        self.ds = ds

    @classmethod
    def open(
        cls,
        url,
        application=None,
        session=None,
        output_grid=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ):
        import pydap.client
        import pydap.lib

        if timeout is None:
            from pydap.lib import DEFAULT_TIMEOUT

            timeout = DEFAULT_TIMEOUT

        kwargs = {
            "url": url,
            "application": application,
            "session": session,
            "output_grid": output_grid or True,
            "timeout": timeout,
        }
        if verify is not None:
            kwargs.update({"verify": verify})
        if user_charset is not None:
            kwargs.update({"user_charset": user_charset})
        ds = pydap.client.open_url(**kwargs)
        return cls(ds)

    def open_store_variable(self, var):
        data = indexing.LazilyIndexedArray(PydapArrayWrapper(var))
        return Variable(var.dimensions, data, _fix_attributes(var.attributes))

    def get_variables(self):
        return FrozenDict(
            (k, self.open_store_variable(self.ds[k])) for k in self.ds.keys()
        )

    def get_attrs(self):
        return Frozen(_fix_attributes(self.ds.attributes))

    def get_dimensions(self):
        return Frozen(self.ds.dimensions)


class PydapBackendEntrypoint(BackendEntrypoint):
    """
    Backend for steaming datasets over the internet using
    the Data Access Protocol, also known as DODS or OPeNDAP
    based on the pydap package.

    This backend is selected by default for urls.

    For more information about the underlying library, visit:
    https://www.pydap.org

    See Also
    --------
    backends.PydapDataStore
    """

    description = "Open remote datasets via OPeNDAP using pydap in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.PydapBackendEntrypoint.html"

    def guess_can_open(
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
    ) -> bool:
        return isinstance(filename_or_obj, str) and is_remote_uri(filename_or_obj)

    def open_dataset(  # type: ignore[override]  # allow LSP violation, not supporting **kwargs
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        application=None,
        session=None,
        output_grid=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ) -> Dataset:
        store = PydapDataStore.open(
            url=filename_or_obj,
            application=application,
            session=session,
            output_grid=output_grid,
            timeout=timeout,
            verify=verify,
            user_charset=user_charset,
        )

        store_entrypoint = StoreBackendEntrypoint()
        with close_on_error(store):
            ds = store_entrypoint.open_dataset(
                store,
                mask_and_scale=mask_and_scale,
                decode_times=decode_times,
                concat_characters=concat_characters,
                decode_coords=decode_coords,
                drop_variables=drop_variables,
                use_cftime=use_cftime,
                decode_timedelta=decode_timedelta,
            )
            return ds


BACKEND_ENTRYPOINTS["pydap"] = ("pydap", PydapBackendEntrypoint)
