from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from xarray.backends.api import open_datatree
from xarray.testing import assert_equal
from xarray.tests import (
    requires_h5netcdf,
    requires_netCDF4,
    requires_zarr,
)

if TYPE_CHECKING:
    from xarray.backends.api import T_NetcdfEngine


class DatatreeIOBase:
    engine: T_NetcdfEngine | None = None

    def test_to_netcdf(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree
        original_dt.to_netcdf(filepath, engine=self.engine)

        roundtrip_dt = open_datatree(filepath, engine=self.engine)
        assert_equal(original_dt, roundtrip_dt)

    def test_netcdf_encoding(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree

        # add compression
        comp = dict(zlib=True, complevel=9)
        enc = {"/set2": {var: comp for var in original_dt["/set2"].ds.data_vars}}

        original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)
        roundtrip_dt = open_datatree(filepath, engine=self.engine)

        assert roundtrip_dt["/set2/a"].encoding["zlib"] == comp["zlib"]
        assert roundtrip_dt["/set2/a"].encoding["complevel"] == comp["complevel"]

        enc["/not/a/group"] = {"foo": "bar"}  # type: ignore
        with pytest.raises(ValueError, match="unexpected encoding group.*"):
            original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)


@requires_netCDF4
class TestNetCDF4DatatreeIO(DatatreeIOBase):
    engine: T_NetcdfEngine | None = "netcdf4"


@requires_h5netcdf
class TestH5NetCDFDatatreeIO(DatatreeIOBase):
    engine: T_NetcdfEngine | None = "h5netcdf"


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
        enc = {"/set2": {var: comp for var in original_dt["/set2"].ds.data_vars}}
        original_dt.to_zarr(filepath, encoding=enc)
        roundtrip_dt = open_datatree(filepath, engine="zarr")

        print(roundtrip_dt["/set2/a"].encoding)
        assert roundtrip_dt["/set2/a"].encoding["compressor"] == comp["compressor"]

        enc["/not/a/group"] = {"foo": "bar"}  # type: ignore
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
