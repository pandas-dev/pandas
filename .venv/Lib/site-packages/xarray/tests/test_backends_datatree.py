from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import pytest

import xarray as xr
from xarray.backends.api import open_datatree, open_groups
from xarray.core.datatree import DataTree
from xarray.testing import assert_equal, assert_identical
from xarray.tests import (
    requires_h5netcdf,
    requires_netCDF4,
    requires_zarr,
)

if TYPE_CHECKING:
    from xarray.core.datatree_io import T_DataTreeNetcdfEngine

try:
    import netCDF4 as nc4
except ImportError:
    pass


class DatatreeIOBase:
    engine: T_DataTreeNetcdfEngine | None = None

    def test_to_netcdf(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree
        original_dt.to_netcdf(filepath, engine=self.engine)

        roundtrip_dt = open_datatree(filepath, engine=self.engine)
        assert_equal(original_dt, roundtrip_dt)

    def test_to_netcdf_inherited_coords(self, tmpdir):
        filepath = tmpdir / "test.nc"
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"a": (("x",), [1, 2])}, coords={"x": [3, 4]}),
                "/sub": xr.Dataset({"b": (("x",), [5, 6])}),
            }
        )
        original_dt.to_netcdf(filepath, engine=self.engine)

        roundtrip_dt = open_datatree(filepath, engine=self.engine)
        assert_equal(original_dt, roundtrip_dt)
        subtree = cast(DataTree, roundtrip_dt["/sub"])
        assert "x" not in subtree.to_dataset(inherited=False).coords

    def test_netcdf_encoding(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree

        # add compression
        comp = dict(zlib=True, complevel=9)
        enc = {"/set2": {var: comp for var in original_dt["/set2"].dataset.data_vars}}

        original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)
        roundtrip_dt = open_datatree(filepath, engine=self.engine)

        assert roundtrip_dt["/set2/a"].encoding["zlib"] == comp["zlib"]
        assert roundtrip_dt["/set2/a"].encoding["complevel"] == comp["complevel"]

        enc["/not/a/group"] = {"foo": "bar"}  # type: ignore[dict-item]
        with pytest.raises(ValueError, match="unexpected encoding group.*"):
            original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)


@requires_netCDF4
class TestNetCDF4DatatreeIO(DatatreeIOBase):
    engine: T_DataTreeNetcdfEngine | None = "netcdf4"

    def test_open_datatree(self, tmpdir) -> None:
        """Create a test netCDF4 file with this unaligned structure:
        Group: /
        │   Dimensions:        (lat: 1, lon: 2)
        │   Dimensions without coordinates: lat, lon
        │   Data variables:
        │       root_variable  (lat, lon) float64 16B ...
        └── Group: /Group1
            │   Dimensions:      (lat: 1, lon: 2)
            │   Dimensions without coordinates: lat, lon
            │   Data variables:
            │       group_1_var  (lat, lon) float64 16B ...
            └── Group: /Group1/subgroup1
                    Dimensions:        (lat: 2, lon: 2)
                    Dimensions without coordinates: lat, lon
                    Data variables:
                        subgroup1_var  (lat, lon) float64 32B ...
        """
        filepath = tmpdir + "/unaligned_subgroups.nc"
        with nc4.Dataset(filepath, "w", format="NETCDF4") as root_group:
            group_1 = root_group.createGroup("/Group1")
            subgroup_1 = group_1.createGroup("/subgroup1")

            root_group.createDimension("lat", 1)
            root_group.createDimension("lon", 2)
            root_group.createVariable("root_variable", np.float64, ("lat", "lon"))

            group_1_var = group_1.createVariable(
                "group_1_var", np.float64, ("lat", "lon")
            )
            group_1_var[:] = np.array([[0.1, 0.2]])
            group_1_var.units = "K"
            group_1_var.long_name = "air_temperature"

            subgroup_1.createDimension("lat", 2)

            subgroup1_var = subgroup_1.createVariable(
                "subgroup1_var", np.float64, ("lat", "lon")
            )
            subgroup1_var[:] = np.array([[0.1, 0.2]])
        with pytest.raises(ValueError):
            open_datatree(filepath)

    def test_open_groups(self, tmpdir) -> None:
        """Test `open_groups` with netCDF4 file with the same unaligned structure:
        Group: /
        │   Dimensions:        (lat: 1, lon: 2)
        │   Dimensions without coordinates: lat, lon
        │   Data variables:
        │       root_variable  (lat, lon) float64 16B ...
        └── Group: /Group1
            │   Dimensions:      (lat: 1, lon: 2)
            │   Dimensions without coordinates: lat, lon
            │   Data variables:
            │       group_1_var  (lat, lon) float64 16B ...
            └── Group: /Group1/subgroup1
                    Dimensions:        (lat: 2, lon: 2)
                    Dimensions without coordinates: lat, lon
                    Data variables:
                        subgroup1_var  (lat, lon) float64 32B ...
        """
        filepath = tmpdir + "/unaligned_subgroups.nc"
        with nc4.Dataset(filepath, "w", format="NETCDF4") as root_group:
            group_1 = root_group.createGroup("/Group1")
            subgroup_1 = group_1.createGroup("/subgroup1")

            root_group.createDimension("lat", 1)
            root_group.createDimension("lon", 2)
            root_group.createVariable("root_variable", np.float64, ("lat", "lon"))

            group_1_var = group_1.createVariable(
                "group_1_var", np.float64, ("lat", "lon")
            )
            group_1_var[:] = np.array([[0.1, 0.2]])
            group_1_var.units = "K"
            group_1_var.long_name = "air_temperature"

            subgroup_1.createDimension("lat", 2)

            subgroup1_var = subgroup_1.createVariable(
                "subgroup1_var", np.float64, ("lat", "lon")
            )
            subgroup1_var[:] = np.array([[0.1, 0.2]])

        unaligned_dict_of_datasets = open_groups(filepath)

        # Check that group names are keys in the dictionary of `xr.Datasets`
        assert "/" in unaligned_dict_of_datasets.keys()
        assert "/Group1" in unaligned_dict_of_datasets.keys()
        assert "/Group1/subgroup1" in unaligned_dict_of_datasets.keys()
        # Check that group name returns the correct datasets
        assert_identical(
            unaligned_dict_of_datasets["/"], xr.open_dataset(filepath, group="/")
        )
        assert_identical(
            unaligned_dict_of_datasets["/Group1"],
            xr.open_dataset(filepath, group="Group1"),
        )
        assert_identical(
            unaligned_dict_of_datasets["/Group1/subgroup1"],
            xr.open_dataset(filepath, group="/Group1/subgroup1"),
        )

    def test_open_groups_to_dict(self, tmpdir) -> None:
        """Create a an aligned netCDF4 with the following structure to test `open_groups`
        and `DataTree.from_dict`.
        Group: /
        │   Dimensions:        (lat: 1, lon: 2)
        │   Dimensions without coordinates: lat, lon
        │   Data variables:
        │       root_variable  (lat, lon) float64 16B ...
        └── Group: /Group1
            │   Dimensions:      (lat: 1, lon: 2)
            │   Dimensions without coordinates: lat, lon
            │   Data variables:
            │       group_1_var  (lat, lon) float64 16B ...
            └── Group: /Group1/subgroup1
                    Dimensions:        (lat: 1, lon: 2)
                    Dimensions without coordinates: lat, lon
                    Data variables:
                        subgroup1_var  (lat, lon) float64 16B ...
        """
        filepath = tmpdir + "/all_aligned_child_nodes.nc"
        with nc4.Dataset(filepath, "w", format="NETCDF4") as root_group:
            group_1 = root_group.createGroup("/Group1")
            subgroup_1 = group_1.createGroup("/subgroup1")

            root_group.createDimension("lat", 1)
            root_group.createDimension("lon", 2)
            root_group.createVariable("root_variable", np.float64, ("lat", "lon"))

            group_1_var = group_1.createVariable(
                "group_1_var", np.float64, ("lat", "lon")
            )
            group_1_var[:] = np.array([[0.1, 0.2]])
            group_1_var.units = "K"
            group_1_var.long_name = "air_temperature"

            subgroup1_var = subgroup_1.createVariable(
                "subgroup1_var", np.float64, ("lat", "lon")
            )
            subgroup1_var[:] = np.array([[0.1, 0.2]])

        aligned_dict_of_datasets = open_groups(filepath)
        aligned_dt = DataTree.from_dict(aligned_dict_of_datasets)

        assert open_datatree(filepath).identical(aligned_dt)


@requires_h5netcdf
class TestH5NetCDFDatatreeIO(DatatreeIOBase):
    engine: T_DataTreeNetcdfEngine | None = "h5netcdf"


@requires_zarr
class TestZarrDatatreeIO:
    engine = "zarr"

    def test_to_zarr(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.zarr"
        original_dt = simple_datatree
        original_dt.to_zarr(filepath)

        roundtrip_dt = open_datatree(filepath, engine="zarr")
        assert_equal(original_dt, roundtrip_dt)

    def test_zarr_encoding(self, tmpdir, simple_datatree):
        import zarr

        filepath = tmpdir / "test.zarr"
        original_dt = simple_datatree

        comp = {"compressor": zarr.Blosc(cname="zstd", clevel=3, shuffle=2)}
        enc = {"/set2": {var: comp for var in original_dt["/set2"].dataset.data_vars}}
        original_dt.to_zarr(filepath, encoding=enc)
        roundtrip_dt = open_datatree(filepath, engine="zarr")

        print(roundtrip_dt["/set2/a"].encoding)
        assert roundtrip_dt["/set2/a"].encoding["compressor"] == comp["compressor"]

        enc["/not/a/group"] = {"foo": "bar"}  # type: ignore[dict-item]
        with pytest.raises(ValueError, match="unexpected encoding group.*"):
            original_dt.to_zarr(filepath, encoding=enc, engine="zarr")

    def test_to_zarr_zip_store(self, tmpdir, simple_datatree):
        from zarr.storage import ZipStore

        filepath = tmpdir / "test.zarr.zip"
        original_dt = simple_datatree
        store = ZipStore(filepath)
        original_dt.to_zarr(store)

        roundtrip_dt = open_datatree(store, engine="zarr")
        assert_equal(original_dt, roundtrip_dt)

    def test_to_zarr_not_consolidated(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.zarr"
        zmetadata = filepath / ".zmetadata"
        s1zmetadata = filepath / "set1" / ".zmetadata"
        filepath = str(filepath)  # casting to str avoids a pathlib bug in xarray
        original_dt = simple_datatree
        original_dt.to_zarr(filepath, consolidated=False)
        assert not zmetadata.exists()
        assert not s1zmetadata.exists()

        with pytest.warns(RuntimeWarning, match="consolidated"):
            roundtrip_dt = open_datatree(filepath, engine="zarr")
        assert_equal(original_dt, roundtrip_dt)

    def test_to_zarr_default_write_mode(self, tmpdir, simple_datatree):
        import zarr

        simple_datatree.to_zarr(tmpdir)

        # with default settings, to_zarr should not overwrite an existing dir
        with pytest.raises(zarr.errors.ContainsGroupError):
            simple_datatree.to_zarr(tmpdir)

    def test_to_zarr_inherited_coords(self, tmpdir):
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"a": (("x",), [1, 2])}, coords={"x": [3, 4]}),
                "/sub": xr.Dataset({"b": (("x",), [5, 6])}),
            }
        )
        filepath = tmpdir / "test.zarr"
        original_dt.to_zarr(filepath)

        roundtrip_dt = open_datatree(filepath, engine="zarr")
        assert_equal(original_dt, roundtrip_dt)
        subtree = cast(DataTree, roundtrip_dt["/sub"])
        assert "x" not in subtree.to_dataset(inherited=False).coords
