from __future__ import annotations

import contextlib
import re
import sys
from collections.abc import Callable, Generator, Hashable
from pathlib import Path
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
import pytest
from packaging.version import Version

import xarray as xr
from xarray import DataTree, load_datatree, open_datatree, open_groups
from xarray.testing import assert_equal, assert_identical
from xarray.tests import (
    has_zarr_v3,
    network,
    parametrize_zarr_format,
    requires_dask,
    requires_h5netcdf,
    requires_h5netcdf_or_netCDF4,
    requires_netCDF4,
    requires_pydap,
    requires_zarr,
)
from xarray.tests.test_backends import TestNetCDF4Data as _TestNetCDF4Data

if TYPE_CHECKING:
    from xarray.backends.writers import T_DataTreeNetcdfEngine

with contextlib.suppress(ImportError):
    import netCDF4 as nc4


ON_WINDOWS = sys.platform == "win32"


class TestNetCDF4DataTree(_TestNetCDF4Data):
    @contextlib.contextmanager
    def open(self, path, **kwargs):
        with open_datatree(path, engine=self.engine, **kwargs) as ds:
            yield ds.to_dataset()

    def test_child_group_with_inconsistent_dimensions(self) -> None:
        with pytest.raises(
            ValueError, match=r"group '/child' is not aligned with its parents"
        ):
            super().test_child_group_with_inconsistent_dimensions()


def diff_chunks(
    comparison: dict[tuple[str, Hashable], bool], tree1: DataTree, tree2: DataTree
) -> str:
    mismatching_variables = [loc for loc, equals in comparison.items() if not equals]

    variable_messages = [
        "\n".join(
            [
                f"L  {path}:{name}: {tree1[path].variables[name].chunksizes}",
                f"R  {path}:{name}: {tree2[path].variables[name].chunksizes}",
            ]
        )
        for path, name in mismatching_variables
    ]
    return "\n".join(["Differing chunk sizes:"] + variable_messages)


def assert_chunks_equal(
    actual: DataTree, expected: DataTree, enforce_dask: bool = False
) -> None:
    __tracebackhide__ = True

    from xarray.namedarray.pycompat import array_type

    dask_array_type = array_type("dask")

    comparison = {
        (path, name): (
            (
                not enforce_dask
                or isinstance(node1.variables[name].data, dask_array_type)
            )
            and node1.variables[name].chunksizes == node2.variables[name].chunksizes
        )
        for path, (node1, node2) in xr.group_subtrees(actual, expected)
        for name in node1.variables.keys()
    }

    assert all(comparison.values()), diff_chunks(comparison, actual, expected)


@pytest.fixture(scope="module")
def unaligned_datatree_nc(tmp_path_factory):
    """Creates a test netCDF4 file with the following unaligned structure, writes it to a /tmp directory
    and returns the file path of the netCDF4 file.

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
    filepath = tmp_path_factory.mktemp("data") / "unaligned_subgroups.nc"
    with nc4.Dataset(filepath, "w", format="NETCDF4") as root_group:
        group_1 = root_group.createGroup("/Group1")
        subgroup_1 = group_1.createGroup("/subgroup1")

        root_group.createDimension("lat", 1)
        root_group.createDimension("lon", 2)
        root_group.createVariable("root_variable", np.float64, ("lat", "lon"))

        group_1_var = group_1.createVariable("group_1_var", np.float64, ("lat", "lon"))
        group_1_var[:] = np.array([[0.1, 0.2]])
        group_1_var.units = "K"
        group_1_var.long_name = "air_temperature"

        subgroup_1.createDimension("lat", 2)

        subgroup1_var = subgroup_1.createVariable(
            "subgroup1_var", np.float64, ("lat", "lon")
        )
        subgroup1_var[:] = np.array([[0.1, 0.2]])

    yield filepath


@pytest.fixture(scope="module")
def unaligned_datatree_zarr_factory(
    tmp_path_factory,
) -> Generator[
    Callable[[Literal[2, 3]], Path],
    None,
    None,
]:
    """Creates a zarr store with the following unaligned group hierarchy:
    Group: /
    │   Dimensions:  (y: 3, x: 2)
    │   Dimensions without coordinates: y, x
    │   Data variables:
    │       a        (y) int64 24B ...
    │       set0     (x) int64 16B ...
    └── Group: /Group1
    │   │   Dimensions:  ()
    │   │   Data variables:
    │   │       a        int64 8B ...
    │   │       b        int64 8B ...
    │   └── /Group1/subgroup1
    │           Dimensions:  ()
    │           Data variables:
    │               a        int64 8B ...
    │               b        int64 8B ...
    └── Group: /Group2
            Dimensions:  (y: 2, x: 2)
            Dimensions without coordinates: y, x
            Data variables:
                a        (y) int64 16B ...
                b        (x) float64 16B ...
    """

    def _unaligned_datatree_zarr(zarr_format: Literal[2, 3]) -> Path:
        filepath = tmp_path_factory.mktemp("data") / "unaligned_simple_datatree.zarr"
        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": 0, "b": 1})
        set2_data = xr.Dataset({"a": ("y", [2, 3]), "b": ("x", [0.1, 0.2])})

        root_data.to_zarr(
            filepath,
            mode="w",
            zarr_format=zarr_format,
        )
        set1_data.to_zarr(
            filepath,
            group="/Group1",
            mode="a",
            zarr_format=zarr_format,
        )
        set2_data.to_zarr(
            filepath,
            group="/Group2",
            mode="a",
            zarr_format=zarr_format,
        )
        set1_data.to_zarr(
            filepath,
            group="/Group1/subgroup1",
            mode="a",
            zarr_format=zarr_format,
        )

        return filepath

    yield _unaligned_datatree_zarr


class NetCDFIOBase:
    engine: T_DataTreeNetcdfEngine | None

    def test_to_netcdf(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree
        original_dt.to_netcdf(filepath, engine=self.engine)

        with open_datatree(filepath, engine=self.engine) as roundtrip_dt:
            assert roundtrip_dt._close is not None
            assert_equal(original_dt, roundtrip_dt)

    def test_decode_cf(self, tmpdir):
        filepath = tmpdir / "test-cf-convention.nc"
        original_dt = xr.DataTree(
            xr.Dataset(
                {
                    "test": xr.DataArray(
                        data=np.array([0, 1, 2], dtype=np.uint16),
                        attrs={"_FillValue": 99},
                    ),
                }
            )
        )
        original_dt.to_netcdf(filepath, engine=self.engine)
        with open_datatree(
            filepath, engine=self.engine, decode_cf=False
        ) as roundtrip_dt:
            assert original_dt["test"].dtype == roundtrip_dt["test"].dtype

    def test_to_netcdf_inherited_coords(self, tmpdir) -> None:
        filepath = tmpdir / "test.nc"
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"a": (("x",), [1, 2])}, coords={"x": [3, 4]}),
                "/sub": xr.Dataset({"b": (("x",), [5, 6])}),
            }
        )
        original_dt.to_netcdf(filepath, engine=self.engine)

        with open_datatree(filepath, engine=self.engine) as roundtrip_dt:
            assert_equal(original_dt, roundtrip_dt)
            subtree = cast(DataTree, roundtrip_dt["/sub"])
            assert "x" not in subtree.to_dataset(inherit=False).coords

    def test_netcdf_encoding(self, tmpdir, simple_datatree) -> None:
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree

        # add compression
        comp = dict(zlib=True, complevel=9)
        enc = {"/set2": dict.fromkeys(original_dt["/set2"].dataset.data_vars, comp)}

        original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)
        with open_datatree(filepath, engine=self.engine) as roundtrip_dt:
            assert roundtrip_dt["/set2/a"].encoding["zlib"] == comp["zlib"]
            assert roundtrip_dt["/set2/a"].encoding["complevel"] == comp["complevel"]

            enc["/not/a/group"] = {"foo": "bar"}  # type: ignore[dict-item]
            with pytest.raises(ValueError, match=r"unexpected encoding group.*"):
                original_dt.to_netcdf(filepath, encoding=enc, engine=self.engine)

    def test_write_subgroup(self, tmpdir) -> None:
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        ).children["child"]

        expected_dt = original_dt.copy()
        expected_dt.name = None

        filepath = tmpdir / "test.zarr"
        original_dt.to_netcdf(filepath, engine=self.engine)

        with open_datatree(filepath, engine=self.engine) as roundtrip_dt:
            assert_equal(original_dt, roundtrip_dt)
            assert_identical(expected_dt, roundtrip_dt)

    @requires_netCDF4
    def test_no_redundant_dimensions(self, tmpdir) -> None:
        # regression test for https://github.com/pydata/xarray/issues/10241
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        )
        filepath = tmpdir / "test.zarr"
        original_dt.to_netcdf(filepath, engine=self.engine)

        root = nc4.Dataset(str(filepath))
        child = root.groups["child"]
        assert list(root.dimensions) == ["x"]
        assert list(child.dimensions) == []

    @requires_dask
    def test_compute_false(self, tmpdir, simple_datatree):
        filepath = tmpdir / "test.nc"
        original_dt = simple_datatree.chunk()
        result = original_dt.to_netcdf(filepath, engine=self.engine, compute=False)

        if not ON_WINDOWS:
            # File at filepath is not closed until .compute() is called. On
            # Windows, this means we can't open it yet.
            with open_datatree(filepath, engine=self.engine) as in_progress_dt:
                assert in_progress_dt.isomorphic(original_dt)
                assert not in_progress_dt.equals(original_dt)

        result.compute()
        with open_datatree(filepath, engine=self.engine) as written_dt:
            assert_identical(written_dt, original_dt)

    def test_default_write_engine(self, tmpdir, simple_datatree, monkeypatch):
        # Ensure the other netCDF library are not installed
        exclude = "netCDF4" if self.engine == "h5netcdf" else "h5netcdf"
        monkeypatch.delitem(sys.modules, exclude, raising=False)
        monkeypatch.setattr(sys, "meta_path", [])

        filepath = tmpdir + "/phony_dims.nc"
        original_dt = simple_datatree
        original_dt.to_netcdf(filepath)  # should not raise

    @requires_dask
    def test_open_datatree_chunks(self, tmpdir) -> None:
        filepath = tmpdir / "test.nc"

        chunks = {"x": 2, "y": 1}

        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": ("y", [-1, 0, 1]), "b": ("x", [-10, 6])})
        set2_data = xr.Dataset({"a": ("y", [1, 2, 3]), "b": ("x", [0.1, 0.2])})
        original_tree = DataTree.from_dict(
            {
                "/": root_data.chunk(chunks),
                "/group1": set1_data.chunk(chunks),
                "/group2": set2_data.chunk(chunks),
            }
        )
        original_tree.to_netcdf(filepath, engine=self.engine)

        with open_datatree(filepath, engine=self.engine, chunks=chunks) as tree:
            xr.testing.assert_identical(tree, original_tree)

            assert_chunks_equal(tree, original_tree, enforce_dask=True)

    def test_roundtrip_via_memoryview(self, simple_datatree) -> None:
        original_dt = simple_datatree
        memview = original_dt.to_netcdf(engine=self.engine)
        roundtrip_dt = load_datatree(memview, engine=self.engine)
        assert_equal(original_dt, roundtrip_dt)

    def test_to_memoryview_compute_false(self, simple_datatree) -> None:
        original_dt = simple_datatree
        with pytest.raises(
            NotImplementedError,
            match=re.escape("to_netcdf() with compute=False is not yet implemented"),
        ):
            original_dt.to_netcdf(engine=self.engine, compute=False)

    def test_open_datatree_specific_group(self, tmpdir, simple_datatree) -> None:
        """Test opening a specific group within a NetCDF file using `open_datatree`."""
        filepath = tmpdir / "test.nc"
        group = "/set1"
        original_dt = simple_datatree
        original_dt.to_netcdf(filepath, engine=self.engine)
        expected_subtree = original_dt[group].copy()
        expected_subtree.orphan()
        with open_datatree(filepath, group=group, engine=self.engine) as subgroup_tree:
            assert subgroup_tree.root.parent is None
            assert_equal(subgroup_tree, expected_subtree)


@requires_h5netcdf_or_netCDF4
class TestGenericNetCDFIO(NetCDFIOBase):
    engine: T_DataTreeNetcdfEngine | None = None

    @requires_netCDF4
    def test_open_netcdf3(self, tmpdir) -> None:
        filepath = tmpdir / "test.nc"
        ds = xr.Dataset({"foo": 1})
        ds.to_netcdf(filepath, format="NETCDF3_CLASSIC")

        expected_dt = DataTree(ds)
        roundtrip_dt = load_datatree(filepath)  # must use netCDF4 engine
        assert_equal(expected_dt, roundtrip_dt)

    @requires_h5netcdf
    @requires_netCDF4
    def test_memoryview_write_h5netcdf_read_netcdf4(self, simple_datatree) -> None:
        original_dt = simple_datatree
        memview = original_dt.to_netcdf(engine="h5netcdf")
        roundtrip_dt = load_datatree(memview, engine="netcdf4")
        assert_equal(original_dt, roundtrip_dt)

    @requires_h5netcdf
    @requires_netCDF4
    def test_memoryview_write_netcdf4_read_h5netcdf(self, simple_datatree) -> None:
        original_dt = simple_datatree
        memview = original_dt.to_netcdf(engine="netcdf4")
        roundtrip_dt = load_datatree(memview, engine="h5netcdf")
        assert_equal(original_dt, roundtrip_dt)

    def test_open_datatree_unaligned_hierarchy(self, unaligned_datatree_nc) -> None:
        with pytest.raises(
            ValueError,
            match=(
                re.escape(
                    "group '/Group1/subgroup1' is not aligned with its parents:\nGroup:\n"
                )
                + ".*"
            ),
        ):
            open_datatree(unaligned_datatree_nc)

    def test_open_groups(self, unaligned_datatree_nc) -> None:
        """Test `open_groups` with a netCDF4 file with an unaligned group hierarchy."""
        unaligned_dict_of_datasets = open_groups(unaligned_datatree_nc)

        # Check that group names are keys in the dictionary of `xr.Datasets`
        assert "/" in unaligned_dict_of_datasets.keys()
        assert "/Group1" in unaligned_dict_of_datasets.keys()
        assert "/Group1/subgroup1" in unaligned_dict_of_datasets.keys()
        # Check that group name returns the correct datasets
        with xr.open_dataset(unaligned_datatree_nc, group="/") as expected:
            assert_identical(unaligned_dict_of_datasets["/"], expected)
        with xr.open_dataset(unaligned_datatree_nc, group="Group1") as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1"], expected)
        with xr.open_dataset(
            unaligned_datatree_nc, group="/Group1/subgroup1"
        ) as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1/subgroup1"], expected)

        for ds in unaligned_dict_of_datasets.values():
            ds.close()

    @requires_dask
    def test_open_groups_chunks(self, tmpdir) -> None:
        """Test `open_groups` with chunks on a netcdf4 file."""

        chunks = {"x": 2, "y": 1}
        filepath = tmpdir / "test.nc"

        chunks = {"x": 2, "y": 1}

        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": ("y", [-1, 0, 1]), "b": ("x", [-10, 6])})
        set2_data = xr.Dataset({"a": ("y", [1, 2, 3]), "b": ("x", [0.1, 0.2])})
        original_tree = DataTree.from_dict(
            {
                "/": root_data.chunk(chunks),
                "/group1": set1_data.chunk(chunks),
                "/group2": set2_data.chunk(chunks),
            }
        )
        original_tree.to_netcdf(filepath, mode="w")

        dict_of_datasets = open_groups(filepath, chunks=chunks)

        for path, ds in dict_of_datasets.items():
            assert {k: max(vs) for k, vs in ds.chunksizes.items()} == chunks, (
                f"unexpected chunking for {path}"
            )

        for ds in dict_of_datasets.values():
            ds.close()


@requires_netCDF4
class TestNetCDF4DatatreeIO(NetCDFIOBase):
    engine: T_DataTreeNetcdfEngine | None = "netcdf4"

    def test_open_groups_to_dict(self, tmpdir) -> None:
        """Create an aligned netCDF4 with the following structure to test `open_groups`
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
        with open_datatree(filepath) as opened_tree:
            assert opened_tree.identical(aligned_dt)
        for ds in aligned_dict_of_datasets.values():
            ds.close()


@requires_h5netcdf
class TestH5NetCDFDatatreeIO(NetCDFIOBase):
    engine: T_DataTreeNetcdfEngine | None = "h5netcdf"

    def test_phony_dims_warning(self, tmpdir) -> None:
        filepath = tmpdir + "/phony_dims.nc"
        import h5py

        foo_data = np.arange(125).reshape(5, 5, 5)
        bar_data = np.arange(625).reshape(25, 5, 5)
        var = {"foo1": foo_data, "foo2": bar_data, "foo3": foo_data, "foo4": bar_data}
        with h5py.File(filepath, "w") as f:
            grps = ["bar", "baz"]
            for grp in grps:
                fx = f.create_group(grp)
                for k, v in var.items():
                    fx.create_dataset(k, data=v)

        with pytest.warns(UserWarning, match="The 'phony_dims' kwarg"):
            with open_datatree(filepath, engine=self.engine) as tree:
                assert tree.bar.dims == {
                    "phony_dim_0": 5,
                    "phony_dim_1": 5,
                    "phony_dim_2": 5,
                    "phony_dim_3": 25,
                }

    def test_roundtrip_using_filelike_object(self, tmpdir, simple_datatree) -> None:
        original_dt = simple_datatree
        filepath = tmpdir + "/test.nc"
        # h5py requires both read and write access when writing, it will
        # work with file-like objects provided they support both, and are
        # seekable.
        with open(filepath, "wb+") as file:
            original_dt.to_netcdf(file, engine=self.engine)
        with open(filepath, "rb") as file:
            with open_datatree(file, engine=self.engine) as roundtrip_dt:
                assert_equal(original_dt, roundtrip_dt)


@network
@requires_pydap
class TestPyDAPDatatreeIO:
    """Test PyDAP backend for DataTree."""

    engine: T_DataTreeNetcdfEngine | None = "pydap"
    # you can check these by adding a .dmr to urls, and replacing dap4 with http
    unaligned_datatree_url = (
        "dap4://test.opendap.org/opendap/dap4/unaligned_simple_datatree.nc.h5"
    )
    all_aligned_child_nodes_url = (
        "dap4://test.opendap.org/opendap/dap4/all_aligned_child_nodes.nc.h5"
    )
    simplegroup_datatree_url = "dap4://test.opendap.org/opendap/dap4/SimpleGroup.nc4.h5"

    def test_open_datatree_unaligned_hierarchy(
        self,
        url=unaligned_datatree_url,
    ) -> None:
        with pytest.raises(
            ValueError,
            match=(
                re.escape(
                    "group '/Group1/subgroup1' is not aligned with its parents:\nGroup:\n"
                )
                + ".*"
            ),
        ):
            open_datatree(url, engine=self.engine)

    def test_open_groups(self, url=unaligned_datatree_url) -> None:
        """Test `open_groups` with a netCDF4/HDF5 file with an unaligned group hierarchy."""
        unaligned_dict_of_datasets = open_groups(url, engine=self.engine)

        # Check that group names are keys in the dictionary of `xr.Datasets`
        assert "/" in unaligned_dict_of_datasets.keys()
        assert "/Group1" in unaligned_dict_of_datasets.keys()
        assert "/Group1/subgroup1" in unaligned_dict_of_datasets.keys()
        # Check that group name returns the correct datasets
        with xr.open_dataset(url, engine=self.engine, group="/") as expected:
            assert_identical(unaligned_dict_of_datasets["/"], expected)
        with xr.open_dataset(url, group="Group1", engine=self.engine) as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1"], expected)
        with xr.open_dataset(
            url,
            group="/Group1/subgroup1",
            engine=self.engine,
        ) as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1/subgroup1"], expected)

    def test_inherited_coords(self, tmpdir, url=simplegroup_datatree_url) -> None:
        """Test that `open_datatree` inherits coordinates from root tree.

        This particular h5 file is a test file that inherits the time coordinate from the root
        dataset to the child dataset.

        Group: /
        │   Dimensions:        (time: 1, Z: 1000, nv: 2)
        │   Coordinates:
        |       time: (time)    float32 0.5
        |       Z:    (Z)       float32 -0.0 -1.0 -2.0 ...
        │   Data variables:
        │       Pressure  (Z)   float32 ...
        |       time_bnds (time, nv) float32 ...
        └── Group: /SimpleGroup
            │   Dimensions:      (time: 1, Z: 1000, nv: 2, Y: 40, X: 40)
            │   Coordinates:
            |      Y:   (Y)     int16 1 2 3 4 ...
            |      X:   (X)     int16 1 2 3 4 ...
            |   Inherited coordinates:
            |      time: (time)    float32 0.5
            |      Z:    (Z)       float32 -0.0 -1.0 -2.0 ...
            │   Data variables:
            │       Temperature  (time, Z, Y, X) float32 ...
            |       Salinity     (time, Z, Y, X) float32 ...
        """
        import pydap
        from pydap.net import create_session

        # Create a session with pre-set retry params in pydap backend, to cache urls
        cache_name = tmpdir / "debug"
        session = create_session(
            use_cache=True, cache_kwargs={"cache_name": cache_name}
        )
        session.cache.clear()

        _version_ = Version(pydap.__version__)

        tree = open_datatree(url, engine=self.engine, session=session)
        assert set(tree.dims) == {"time", "Z", "nv"}
        assert tree["/SimpleGroup"].coords["time"].dims == ("time",)
        assert tree["/SimpleGroup"].coords["Z"].dims == ("Z",)
        assert tree["/SimpleGroup"].coords["Y"].dims == ("Y",)
        assert tree["/SimpleGroup"].coords["X"].dims == ("X",)
        with xr.open_dataset(url, engine=self.engine, group="/SimpleGroup") as expected:
            assert set(tree["/SimpleGroup"].dims) == set(
                list(expected.dims) + ["Z", "nv"]
            )

        if _version_ > Version("3.5.5"):
            # Total downloads are: 1 dmr, + 1 dap url for all dimensions for each group
            assert len(session.cache.urls()) == 3
        else:
            # 1 dmr + 1 dap url per dimension (total there are 4 dimension arrays)
            assert len(session.cache.urls()) == 5

    def test_open_groups_to_dict(self, url=all_aligned_child_nodes_url) -> None:
        aligned_dict_of_datasets = open_groups(url, engine=self.engine)
        aligned_dt = DataTree.from_dict(aligned_dict_of_datasets)
        with open_datatree(url, engine=self.engine) as opened_tree:
            assert opened_tree.identical(aligned_dt)


@requires_zarr
@parametrize_zarr_format
class TestZarrDatatreeIO:
    engine = "zarr"

    def test_to_zarr(self, tmpdir, simple_datatree, zarr_format) -> None:
        filepath = str(tmpdir / "test.zarr")
        original_dt = simple_datatree
        original_dt.to_zarr(filepath, zarr_format=zarr_format)

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            assert_equal(original_dt, roundtrip_dt)

    @pytest.mark.filterwarnings(
        "ignore:Numcodecs codecs are not in the Zarr version 3 specification"
    )
    def test_zarr_encoding(self, tmpdir, simple_datatree, zarr_format) -> None:
        filepath = str(tmpdir / "test.zarr")
        original_dt = simple_datatree

        if zarr_format == 2:
            from numcodecs.blosc import Blosc

            codec = Blosc(cname="zstd", clevel=3, shuffle=2)
            comp = {"compressors": (codec,)} if has_zarr_v3 else {"compressor": codec}
        elif zarr_format == 3:
            # specifying codecs in zarr_format=3 requires importing from zarr 3 namespace
            from zarr.registry import get_codec_class

            Blosc = get_codec_class("numcodecs.blosc")
            comp = {"compressors": (Blosc(cname="zstd", clevel=3),)}  # type: ignore[call-arg]

        enc = {"/set2": dict.fromkeys(original_dt["/set2"].dataset.data_vars, comp)}
        original_dt.to_zarr(filepath, encoding=enc, zarr_format=zarr_format)

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            compressor_key = "compressors" if has_zarr_v3 else "compressor"
            assert (
                roundtrip_dt["/set2/a"].encoding[compressor_key] == comp[compressor_key]
            )

            enc["/not/a/group"] = {"foo": "bar"}  # type: ignore[dict-item]
            with pytest.raises(ValueError, match=r"unexpected encoding group.*"):
                original_dt.to_zarr(filepath, encoding=enc, zarr_format=zarr_format)

    @pytest.mark.xfail(reason="upstream zarr read-only changes have broken this test")
    @pytest.mark.filterwarnings("ignore:Duplicate name")
    def test_to_zarr_zip_store(self, tmpdir, simple_datatree, zarr_format) -> None:
        from zarr.storage import ZipStore

        filepath = str(tmpdir / "test.zarr.zip")
        original_dt = simple_datatree
        store = ZipStore(filepath, mode="w")
        original_dt.to_zarr(store, zarr_format=zarr_format)

        with open_datatree(store, engine="zarr") as roundtrip_dt:  # type: ignore[arg-type, unused-ignore]
            assert_equal(original_dt, roundtrip_dt)

    def test_to_zarr_not_consolidated(
        self, tmpdir, simple_datatree, zarr_format
    ) -> None:
        filepath = tmpdir / "test.zarr"
        zmetadata = filepath / ".zmetadata"
        s1zmetadata = filepath / "set1" / ".zmetadata"
        filepath = str(filepath)  # casting to str avoids a pathlib bug in xarray
        original_dt = simple_datatree
        original_dt.to_zarr(filepath, consolidated=False, zarr_format=zarr_format)
        assert not zmetadata.exists()
        assert not s1zmetadata.exists()

        with pytest.warns(RuntimeWarning, match="consolidated"):
            with open_datatree(filepath, engine="zarr") as roundtrip_dt:
                assert_equal(original_dt, roundtrip_dt)

    def test_to_zarr_default_write_mode(
        self, tmpdir, simple_datatree, zarr_format
    ) -> None:
        simple_datatree.to_zarr(str(tmpdir), zarr_format=zarr_format)

        import zarr

        # expected exception type changed in zarr-python v2->v3, see https://github.com/zarr-developers/zarr-python/issues/2821
        expected_exception_type = (
            FileExistsError if has_zarr_v3 else zarr.errors.ContainsGroupError
        )

        # with default settings, to_zarr should not overwrite an existing dir
        with pytest.raises(expected_exception_type):
            simple_datatree.to_zarr(str(tmpdir))

    @requires_dask
    def test_to_zarr_compute_false(
        self, tmp_path: Path, simple_datatree: DataTree, zarr_format: Literal[2, 3]
    ) -> None:
        import dask.array as da

        storepath = tmp_path / "test.zarr"
        original_dt = simple_datatree.chunk()
        result = original_dt.to_zarr(
            str(storepath), compute=False, zarr_format=zarr_format
        )

        def assert_expected_zarr_files_exist(
            arr_dir: Path,
            chunks_expected: bool,
            is_scalar: bool,
            zarr_format: Literal[2, 3],
        ) -> None:
            """For one zarr array, check that all expected metadata and chunk data files exist."""
            # TODO: This function is now so complicated that it's practically checking compliance with the whole zarr spec...
            # TODO: Perhaps it would be better to instead trust that zarr-python is spec-compliant and check `DataTree` against zarr-python?
            # TODO: The way to do that would ideally be to use zarr-pythons ability to determine how many chunks have been initialized.

            if zarr_format == 2:
                zarray_file, zattrs_file = (arr_dir / ".zarray"), (arr_dir / ".zattrs")

                assert zarray_file.exists() and zarray_file.is_file()
                assert zattrs_file.exists() and zattrs_file.is_file()

                chunk_file = arr_dir / "0"
                if chunks_expected:
                    # assumes empty chunks were written
                    # (i.e. they did not contain only fill_value and write_empty_chunks was False)
                    assert chunk_file.exists() and chunk_file.is_file()
                else:
                    # either dask array or array of all fill_values
                    assert not chunk_file.exists()
            elif zarr_format == 3:
                metadata_file = arr_dir / "zarr.json"
                assert metadata_file.exists() and metadata_file.is_file()

                chunks_dir = arr_dir / "c"
                chunk_file = chunks_dir / "0"
                if chunks_expected:
                    # assumes empty chunks were written
                    # (i.e. they did not contain only fill_value and write_empty_chunks was False)
                    if is_scalar:
                        # this is the expected behaviour for storing scalars in zarr 3, see https://github.com/pydata/xarray/issues/10147
                        assert chunks_dir.exists() and chunks_dir.is_file()
                    else:
                        assert chunks_dir.exists() and chunks_dir.is_dir()
                        assert chunk_file.exists() and chunk_file.is_file()
                else:
                    assert not chunks_dir.exists()
                    assert not chunk_file.exists()

        DEFAULT_ZARR_FILL_VALUE = 0
        # The default value of write_empty_chunks changed from True->False in zarr-python v2->v3
        WRITE_EMPTY_CHUNKS_DEFAULT = not has_zarr_v3

        for node in original_dt.subtree:
            # inherited variables aren't meant to be written to zarr
            local_node_variables = node.to_dataset(inherit=False).variables
            for name, var in local_node_variables.items():
                var_dir = storepath / node.path.removeprefix("/") / name  # type: ignore[operator]

                assert_expected_zarr_files_exist(
                    arr_dir=var_dir,
                    # don't expect dask.Arrays to be written to disk, as compute=False
                    # also don't expect numpy arrays containing only zarr's fill_value to be written to disk
                    chunks_expected=(
                        not isinstance(var.data, da.Array)
                        and (
                            var.data != DEFAULT_ZARR_FILL_VALUE
                            or WRITE_EMPTY_CHUNKS_DEFAULT
                        )
                    ),
                    is_scalar=not bool(var.dims),
                    zarr_format=zarr_format,
                )

        in_progress_dt = load_datatree(str(storepath), engine="zarr")
        assert not in_progress_dt.equals(original_dt)

        result.compute()
        written_dt = load_datatree(str(storepath), engine="zarr")
        assert_identical(written_dt, original_dt)

    @requires_dask
    def test_rplus_mode(
        self, tmp_path: Path, simple_datatree: DataTree, zarr_format: Literal[2, 3]
    ) -> None:
        storepath = tmp_path / "test.zarr"
        original_dt = simple_datatree.chunk()
        original_dt.to_zarr(storepath, compute=False, zarr_format=zarr_format)
        original_dt.to_zarr(storepath, mode="r+")
        with open_datatree(str(storepath), engine="zarr") as written_dt:
            assert_identical(written_dt, original_dt)

    @requires_dask
    def test_to_zarr_no_redundant_computation(self, tmpdir, zarr_format) -> None:
        import dask.array as da

        eval_count = 0

        def expensive_func(x):
            nonlocal eval_count
            eval_count += 1
            return x + 1

        base = da.random.random((), chunks=())
        derived1 = da.map_blocks(expensive_func, base, meta=np.array((), np.float64))
        derived2 = derived1 + 1  # depends on derived1
        tree = DataTree.from_dict(
            {
                "group1": xr.Dataset({"derived": derived1}),
                "group2": xr.Dataset({"derived": derived2}),
            }
        )

        filepath = str(tmpdir / "test.zarr")
        tree.to_zarr(filepath, zarr_format=zarr_format)
        assert eval_count == 1  # not 2

    def test_to_zarr_inherited_coords(self, tmpdir, zarr_format):
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"a": (("x",), [1, 2])}, coords={"x": [3, 4]}),
                "/sub": xr.Dataset({"b": (("x",), [5, 6])}),
            }
        )
        filepath = str(tmpdir / "test.zarr")
        original_dt.to_zarr(filepath, zarr_format=zarr_format)

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            assert_equal(original_dt, roundtrip_dt)
            subtree = cast(DataTree, roundtrip_dt["/sub"])
            assert "x" not in subtree.to_dataset(inherit=False).coords

    def test_open_groups_round_trip(self, tmpdir, simple_datatree, zarr_format) -> None:
        """Test `open_groups` opens a zarr store with the `simple_datatree` structure."""
        filepath = str(tmpdir / "test.zarr")
        original_dt = simple_datatree
        original_dt.to_zarr(filepath, zarr_format=zarr_format)

        roundtrip_dict = open_groups(filepath, engine="zarr")
        roundtrip_dt = DataTree.from_dict(roundtrip_dict)

        with open_datatree(filepath, engine="zarr") as opened_tree:
            assert opened_tree.identical(roundtrip_dt)

        for ds in roundtrip_dict.values():
            ds.close()

    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    def test_open_datatree_unaligned_hierarchy(
        self, unaligned_datatree_zarr_factory, zarr_format
    ) -> None:
        storepath = unaligned_datatree_zarr_factory(zarr_format=zarr_format)

        with pytest.raises(
            ValueError,
            match=(
                re.escape("group '/Group2' is not aligned with its parents:") + ".*"
            ),
        ):
            open_datatree(storepath, engine="zarr")

    @requires_dask
    def test_open_datatree_chunks(self, tmpdir, zarr_format) -> None:
        filepath = str(tmpdir / "test.zarr")

        chunks = {"x": 2, "y": 1}

        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": ("y", [-1, 0, 1]), "b": ("x", [-10, 6])})
        set2_data = xr.Dataset({"a": ("y", [1, 2, 3]), "b": ("x", [0.1, 0.2])})
        original_tree = DataTree.from_dict(
            {
                "/": root_data.chunk(chunks),
                "/group1": set1_data.chunk(chunks),
                "/group2": set2_data.chunk(chunks),
            }
        )
        original_tree.to_zarr(filepath, zarr_format=zarr_format)

        with open_datatree(filepath, engine="zarr", chunks=chunks) as tree:
            xr.testing.assert_identical(tree, original_tree)
            assert_chunks_equal(tree, original_tree, enforce_dask=True)
            # https://github.com/pydata/xarray/issues/10098
            # If the open tasks are not give unique tokens per node, and the
            # dask graph is computed in one go, data won't be uniquely loaded
            # from each node.
            xr.testing.assert_identical(tree.compute(), original_tree)

    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    def test_open_groups(self, unaligned_datatree_zarr_factory, zarr_format) -> None:
        """Test `open_groups` with a zarr store of an unaligned group hierarchy."""

        storepath = unaligned_datatree_zarr_factory(zarr_format=zarr_format)
        unaligned_dict_of_datasets = open_groups(storepath, engine="zarr")

        assert "/" in unaligned_dict_of_datasets.keys()
        assert "/Group1" in unaligned_dict_of_datasets.keys()
        assert "/Group1/subgroup1" in unaligned_dict_of_datasets.keys()
        assert "/Group2" in unaligned_dict_of_datasets.keys()
        # Check that group name returns the correct datasets
        with xr.open_dataset(storepath, group="/", engine="zarr") as expected:
            assert_identical(unaligned_dict_of_datasets["/"], expected)
        with xr.open_dataset(storepath, group="Group1", engine="zarr") as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1"], expected)
        with xr.open_dataset(
            storepath, group="/Group1/subgroup1", engine="zarr"
        ) as expected:
            assert_identical(unaligned_dict_of_datasets["/Group1/subgroup1"], expected)
        with xr.open_dataset(storepath, group="/Group2", engine="zarr") as expected:
            assert_identical(unaligned_dict_of_datasets["/Group2"], expected)

        for ds in unaligned_dict_of_datasets.values():
            ds.close()

    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    @pytest.mark.parametrize("write_consolidated_metadata", [True, False, None])
    def test_open_datatree_specific_group(
        self,
        tmpdir,
        simple_datatree,
        write_consolidated_metadata,
        zarr_format,
    ) -> None:
        """Test opening a specific group within a Zarr store using `open_datatree`."""
        filepath = str(tmpdir / "test.zarr")
        group = "/set2"
        original_dt = simple_datatree
        original_dt.to_zarr(
            filepath, consolidated=write_consolidated_metadata, zarr_format=zarr_format
        )
        expected_subtree = original_dt[group].copy()
        expected_subtree.orphan()
        with open_datatree(filepath, group=group, engine=self.engine) as subgroup_tree:
            assert subgroup_tree.root.parent is None
            assert_equal(subgroup_tree, expected_subtree)

    @requires_dask
    def test_open_groups_chunks(self, tmpdir, zarr_format) -> None:
        """Test `open_groups` with chunks on a zarr store."""

        chunks = {"x": 2, "y": 1}
        filepath = str(tmpdir / "test.zarr")
        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": ("y", [-1, 0, 1]), "b": ("x", [-10, 6])})
        set2_data = xr.Dataset({"a": ("y", [1, 2, 3]), "b": ("x", [0.1, 0.2])})
        original_tree = DataTree.from_dict(
            {
                "/": root_data.chunk(chunks),
                "/group1": set1_data.chunk(chunks),
                "/group2": set2_data.chunk(chunks),
            }
        )
        original_tree.to_zarr(filepath, mode="w", zarr_format=zarr_format)

        dict_of_datasets = open_groups(filepath, engine="zarr", chunks=chunks)

        for path, ds in dict_of_datasets.items():
            assert {k: max(vs) for k, vs in ds.chunksizes.items()} == chunks, (
                f"unexpected chunking for {path}"
            )

        for ds in dict_of_datasets.values():
            ds.close()

    def test_write_subgroup(self, tmpdir, zarr_format) -> None:
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        ).children["child"]

        expected_dt = original_dt.copy()
        expected_dt.name = None

        filepath = str(tmpdir / "test.zarr")
        original_dt.to_zarr(filepath, zarr_format=zarr_format)

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            assert_equal(original_dt, roundtrip_dt)
            assert_identical(expected_dt, roundtrip_dt)

    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    def test_write_inherited_coords_false(self, tmpdir, zarr_format) -> None:
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        )

        filepath = str(tmpdir / "test.zarr")
        original_dt.to_zarr(
            filepath, write_inherited_coords=False, zarr_format=zarr_format
        )

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            assert_identical(original_dt, roundtrip_dt)

        expected_child = original_dt.children["child"].copy(inherit=False)
        expected_child.name = None
        with open_datatree(filepath, group="child", engine="zarr") as roundtrip_child:
            assert_identical(expected_child, roundtrip_child)

    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    def test_write_inherited_coords_true(self, tmpdir, zarr_format) -> None:
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        )

        filepath = str(tmpdir / "test.zarr")
        original_dt.to_zarr(
            filepath, write_inherited_coords=True, zarr_format=zarr_format
        )

        with open_datatree(filepath, engine="zarr") as roundtrip_dt:
            assert_identical(original_dt, roundtrip_dt)

        expected_child = original_dt.children["child"].copy(inherit=True)
        expected_child.name = None
        with open_datatree(filepath, group="child", engine="zarr") as roundtrip_child:
            assert_identical(expected_child, roundtrip_child)

    @pytest.mark.xfail(
        ON_WINDOWS,
        reason="Permission errors from Zarr: https://github.com/pydata/xarray/pull/10793",
    )
    @pytest.mark.filterwarnings(
        "ignore:Failed to open Zarr store with consolidated metadata:RuntimeWarning"
    )
    def test_zarr_engine_recognised(self, tmpdir, zarr_format) -> None:
        """Test that xarray can guess the zarr backend when the engine is not specified"""
        original_dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2, 3]}),
                "/child": xr.Dataset({"foo": ("x", [4, 5, 6])}),
            }
        )

        filepath = str(tmpdir / "test.zarr")
        original_dt.to_zarr(
            filepath, write_inherited_coords=True, zarr_format=zarr_format
        )

        with open_datatree(filepath) as roundtrip_dt:
            assert_identical(original_dt, roundtrip_dt)
