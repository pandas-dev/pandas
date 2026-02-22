from __future__ import annotations

import io
import re
import sys
from numbers import Number

import numpy as np
import pytest

import xarray as xr
from xarray.backends.writers import get_default_netcdf_write_engine
from xarray.tests import (
    assert_identical,
    assert_no_warnings,
    requires_dask,
    requires_h5netcdf,
    requires_netCDF4,
    requires_scipy,
)


@requires_netCDF4
@requires_scipy
@requires_h5netcdf
def test_get_default_netcdf_write_engine() -> None:
    assert xr.get_options()["netcdf_engine_order"] == ("netcdf4", "h5netcdf", "scipy")

    engine = get_default_netcdf_write_engine("", format=None)
    assert engine == "netcdf4"

    engine = get_default_netcdf_write_engine("", format="NETCDF4")
    assert engine == "netcdf4"

    engine = get_default_netcdf_write_engine("", format="NETCDF4_CLASSIC")
    assert engine == "netcdf4"

    engine = get_default_netcdf_write_engine("", format="NETCDF3_CLASSIC")
    assert engine == "netcdf4"

    engine = get_default_netcdf_write_engine(io.BytesIO(), format=None)
    assert engine == "h5netcdf"

    engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF4")
    assert engine == "h5netcdf"

    engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF3_CLASSIC")
    assert engine == "scipy"

    engine = get_default_netcdf_write_engine("path.zarr#mode=nczarr", format=None)
    assert engine == "netcdf4"

    with xr.set_options(netcdf_engine_order=["netcdf4", "scipy", "h5netcdf"]):
        engine = get_default_netcdf_write_engine(io.BytesIO(), format=None)
        assert engine == "scipy"

        engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF4")
        assert engine == "h5netcdf"

        engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF3_CLASSIC")
        assert engine == "scipy"

    with xr.set_options(netcdf_engine_order=["h5netcdf", "scipy", "netcdf4"]):
        engine = get_default_netcdf_write_engine("", format=None)
        assert engine == "h5netcdf"

        engine = get_default_netcdf_write_engine("", format="NETCDF4")
        assert engine == "h5netcdf"

        engine = get_default_netcdf_write_engine("", format="NETCDF4_CLASSIC")
        assert engine == "netcdf4"

        engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF4")
        assert engine == "h5netcdf"

        engine = get_default_netcdf_write_engine("", format="NETCDF3_CLASSIC")
        assert engine == "scipy"

        engine = get_default_netcdf_write_engine(io.BytesIO(), format="NETCDF3_CLASSIC")
        assert engine == "scipy"


@requires_h5netcdf
def test_default_engine_h5netcdf(monkeypatch):
    """Test the default netcdf engine when h5netcdf is the only importable module."""

    monkeypatch.delitem(sys.modules, "netCDF4", raising=False)
    monkeypatch.delitem(sys.modules, "scipy", raising=False)
    monkeypatch.setattr(sys, "meta_path", [])

    engine = get_default_netcdf_write_engine("", format=None)
    assert engine == "h5netcdf"

    with pytest.raises(
        ValueError,
        match=re.escape(
            "cannot write NetCDF files with format='NETCDF3_CLASSIC' because "
            "none of the suitable backend libraries (SUITABLE_BACKENDS) are installed"
        ).replace("SUITABLE_BACKENDS", r"(scipy, netCDF4)|(netCDF4, scipy)"),
    ):
        get_default_netcdf_write_engine("", format="NETCDF3_CLASSIC")


def test_default_engine_nczarr_no_netcdf4_python(monkeypatch):
    monkeypatch.delitem(sys.modules, "netCDF4", raising=False)
    monkeypatch.setattr(sys, "meta_path", [])

    with pytest.raises(
        ValueError,
        match=re.escape(
            "cannot write NetCDF files in NCZarr format because "
            "none of the suitable backend libraries (netCDF4) are installed"
        ),
    ):
        get_default_netcdf_write_engine("#mode=nczarr", format=None)


def test_custom_engine() -> None:
    expected = xr.Dataset(
        dict(a=2 * np.arange(5)), coords=dict(x=("x", np.arange(5), dict(units="s")))
    )

    class CustomBackend(xr.backends.BackendEntrypoint):
        def open_dataset(
            self,
            filename_or_obj,
            drop_variables=None,
            **kwargs,
        ) -> xr.Dataset:
            return expected.copy(deep=True)

    actual = xr.open_dataset("fake_filename", engine=CustomBackend)
    assert_identical(expected, actual)


def test_multiindex() -> None:
    # GH7139
    # Check that we properly handle backends that change index variables
    dataset = xr.Dataset(coords={"coord1": ["A", "B"], "coord2": [1, 2]})
    dataset = dataset.stack(z=["coord1", "coord2"])

    class MultiindexBackend(xr.backends.BackendEntrypoint):
        def open_dataset(
            self,
            filename_or_obj,
            drop_variables=None,
            **kwargs,
        ) -> xr.Dataset:
            return dataset.copy(deep=True)

    loaded = xr.open_dataset("fake_filename", engine=MultiindexBackend)
    assert_identical(dataset, loaded)


class PassThroughBackendEntrypoint(xr.backends.BackendEntrypoint):
    """Access an object passed to the `open_dataset` method."""

    def open_dataset(self, dataset, *, drop_variables=None):
        """Return the first argument."""
        return dataset


def explicit_chunks(chunks, shape):
    """Return explicit chunks, expanding any integer member to a tuple of integers."""
    # Emulate `dask.array.core.normalize_chunks` but for simpler inputs.
    return tuple(
        (
            (
                (size // chunk) * (chunk,)
                + ((size % chunk,) if size % chunk or size == 0 else ())
            )
            if isinstance(chunk, Number)
            else chunk
        )
        for chunk, size in zip(chunks, shape, strict=True)
    )


@requires_dask
class TestPreferredChunks:
    """Test behaviors related to the backend's preferred chunks."""

    var_name = "data"

    def create_dataset(self, shape, pref_chunks):
        """Return a dataset with a variable with the given shape and preferred chunks."""
        dims = tuple(f"dim_{idx}" for idx in range(len(shape)))
        return xr.Dataset(
            {
                self.var_name: xr.Variable(
                    dims,
                    np.empty(shape, dtype=np.dtype("V1")),
                    encoding={
                        "preferred_chunks": dict(zip(dims, pref_chunks, strict=True))
                    },
                )
            }
        )

    def check_dataset(self, initial, final, expected_chunks):
        assert_identical(initial, final)
        assert final[self.var_name].chunks == expected_chunks

    @pytest.mark.parametrize(
        "shape,pref_chunks",
        [
            # Represent preferred chunking with int.
            ((5,), (2,)),
            # Represent preferred chunking with tuple.
            ((5,), ((2, 2, 1),)),
            # Represent preferred chunking with int in two dims.
            ((5, 6), (4, 2)),
            # Represent preferred chunking with tuple in second dim.
            ((5, 6), (4, (2, 2, 2))),
        ],
    )
    @pytest.mark.parametrize("request_with_empty_map", [False, True])
    def test_honor_chunks(self, shape, pref_chunks, request_with_empty_map):
        """Honor the backend's preferred chunks when opening a dataset."""
        initial = self.create_dataset(shape, pref_chunks)
        # To keep the backend's preferred chunks, the `chunks` argument must be an
        # empty mapping or map dimensions to `None`.
        chunks = (
            {}
            if request_with_empty_map
            else dict.fromkeys(initial[self.var_name].dims, None)
        )
        final = xr.open_dataset(
            initial, engine=PassThroughBackendEntrypoint, chunks=chunks
        )
        self.check_dataset(initial, final, explicit_chunks(pref_chunks, shape))

    @pytest.mark.parametrize(
        "shape,pref_chunks,req_chunks",
        [
            # Preferred chunking is int; requested chunking is int.
            ((5,), (2,), (3,)),
            # Preferred chunking is int; requested chunking is tuple.
            ((5,), (2,), ((2, 1, 1, 1),)),
            # Preferred chunking is tuple; requested chunking is int.
            ((5,), ((2, 2, 1),), (3,)),
            # Preferred chunking is tuple; requested chunking is tuple.
            ((5,), ((2, 2, 1),), ((2, 1, 1, 1),)),
            # Split chunks along a dimension other than the first.
            ((1, 5), (1, 2), (1, 3)),
        ],
    )
    def test_split_chunks(self, shape, pref_chunks, req_chunks):
        """Warn when the requested chunks separate the backend's preferred chunks."""
        initial = self.create_dataset(shape, pref_chunks)
        with pytest.warns(UserWarning):
            final = xr.open_dataset(
                initial,
                engine=PassThroughBackendEntrypoint,
                chunks=dict(zip(initial[self.var_name].dims, req_chunks, strict=True)),
            )
        self.check_dataset(initial, final, explicit_chunks(req_chunks, shape))

    @pytest.mark.parametrize(
        "shape,pref_chunks,req_chunks",
        [
            # Keep preferred chunks using int representation.
            ((5,), (2,), (2,)),
            # Keep preferred chunks using tuple representation.
            ((5,), (2,), ((2, 2, 1),)),
            # Join chunks, leaving a final short chunk.
            ((5,), (2,), (4,)),
            # Join all chunks with an int larger than the dimension size.
            ((5,), (2,), (6,)),
            # Join one chunk using tuple representation.
            ((5,), (1,), ((1, 1, 2, 1),)),
            # Join one chunk using int representation.
            ((5,), ((1, 1, 2, 1),), (2,)),
            # Join multiple chunks using tuple representation.
            ((5,), ((1, 1, 2, 1),), ((2, 3),)),
            # Join chunks in multiple dimensions.
            ((5, 5), (2, (1, 1, 2, 1)), (4, (2, 3))),
        ],
    )
    def test_join_chunks(self, shape, pref_chunks, req_chunks):
        """Don't warn when the requested chunks join or keep the preferred chunks."""
        initial = self.create_dataset(shape, pref_chunks)
        with assert_no_warnings():
            final = xr.open_dataset(
                initial,
                engine=PassThroughBackendEntrypoint,
                chunks=dict(zip(initial[self.var_name].dims, req_chunks, strict=True)),
            )
        self.check_dataset(initial, final, explicit_chunks(req_chunks, shape))

    @pytest.mark.parametrize("create_default_indexes", [True, False])
    def test_default_indexes(self, create_default_indexes):
        """Create default indexes if the backend does not create them."""
        coords = xr.Coordinates({"x": ("x", [0, 1]), "y": list("abc")}, indexes={})
        initial = xr.Dataset({"a": ("x", [1, 2])}, coords=coords)

        with assert_no_warnings():
            final = xr.open_dataset(
                initial,
                engine=PassThroughBackendEntrypoint,
                create_default_indexes=create_default_indexes,
            )

        if create_default_indexes:
            assert all(name in final.xindexes for name in ["x", "y"])
        else:
            assert len(final.xindexes) == 0

    @pytest.mark.parametrize("create_default_indexes", [True, False])
    def test_default_indexes_passthrough(self, create_default_indexes):
        """Allow creating indexes in the backend."""

        initial = xr.Dataset(
            {"a": (["x", "y"], [[1, 2, 3], [4, 5, 6]])},
            coords={"x": ("x", [0, 1]), "y": ("y", list("abc"))},
        ).stack(z=["x", "y"])

        with assert_no_warnings():
            final = xr.open_dataset(
                initial,
                engine=PassThroughBackendEntrypoint,
                create_default_indexes=create_default_indexes,
            )

        assert initial.coords.equals(final.coords)
