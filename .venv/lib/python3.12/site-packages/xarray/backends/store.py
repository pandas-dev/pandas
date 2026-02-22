from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

from xarray import conventions
from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    AbstractDataStore,
    BackendEntrypoint,
    T_PathFileOrDataStore,
)
from xarray.core.coordinates import Coordinates
from xarray.core.dataset import Dataset

if TYPE_CHECKING:
    pass


class StoreBackendEntrypoint(BackendEntrypoint):
    description = "Open AbstractDataStore instances in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.StoreBackendEntrypoint.html"

    def guess_can_open(self, filename_or_obj: T_PathFileOrDataStore) -> bool:
        return isinstance(filename_or_obj, AbstractDataStore)

    def open_dataset(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        set_indexes: bool = True,
        use_cftime=None,
        decode_timedelta=None,
    ) -> Dataset:
        assert isinstance(filename_or_obj, AbstractDataStore)

        vars, attrs = filename_or_obj.load()
        encoding = filename_or_obj.get_encoding()

        vars, attrs, coord_names = conventions.decode_cf_variables(
            vars,
            attrs,
            mask_and_scale=mask_and_scale,
            decode_times=decode_times,
            concat_characters=concat_characters,
            decode_coords=decode_coords,
            drop_variables=drop_variables,
            use_cftime=use_cftime,
            decode_timedelta=decode_timedelta,
        )

        # split data and coordinate variables (promote dimension coordinates)
        data_vars = {}
        coord_vars = {}
        for name, var in vars.items():
            if name in coord_names or var.dims == (name,):
                coord_vars[name] = var
            else:
                data_vars[name] = var

        # explicit Coordinates object with no index passed
        coords = Coordinates(coord_vars, indexes={})

        ds = Dataset(data_vars, coords=coords, attrs=attrs)
        ds.set_close(filename_or_obj.close)
        ds.encoding = encoding

        return ds


BACKEND_ENTRYPOINTS["store"] = (None, StoreBackendEntrypoint)
