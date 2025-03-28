from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from os import PathLike
from typing import Any, Literal, get_args

from xarray.core.datatree import DataTree
from xarray.core.types import NetcdfWriteModes, ZarrWriteModes

T_DataTreeNetcdfEngine = Literal["netcdf4", "h5netcdf"]
T_DataTreeNetcdfTypes = Literal["NETCDF4"]


def _get_nc_dataset_class(engine: T_DataTreeNetcdfEngine | None):
    if engine == "netcdf4":
        from netCDF4 import Dataset
    elif engine == "h5netcdf":
        from h5netcdf.legacyapi import Dataset
    elif engine is None:
        try:
            from netCDF4 import Dataset
        except ImportError:
            from h5netcdf.legacyapi import Dataset
    else:
        raise ValueError(f"unsupported engine: {engine}")
    return Dataset


def _create_empty_netcdf_group(
    filename: str | PathLike,
    group: str,
    mode: NetcdfWriteModes,
    engine: T_DataTreeNetcdfEngine | None,
):
    ncDataset = _get_nc_dataset_class(engine)

    with ncDataset(filename, mode=mode) as rootgrp:
        rootgrp.createGroup(group)


def _datatree_to_netcdf(
    dt: DataTree,
    filepath: str | PathLike,
    mode: NetcdfWriteModes = "w",
    encoding: Mapping[str, Any] | None = None,
    unlimited_dims: Mapping | None = None,
    format: T_DataTreeNetcdfTypes | None = None,
    engine: T_DataTreeNetcdfEngine | None = None,
    group: str | None = None,
    compute: bool = True,
    **kwargs,
):
    """This function creates an appropriate datastore for writing a datatree to
    disk as a netCDF file.

    See `DataTree.to_netcdf` for full API docs.
    """

    if format not in [None, *get_args(T_DataTreeNetcdfTypes)]:
        raise ValueError("to_netcdf only supports the NETCDF4 format")

    if engine not in [None, *get_args(T_DataTreeNetcdfEngine)]:
        raise ValueError("to_netcdf only supports the netcdf4 and h5netcdf engines")

    if group is not None:
        raise NotImplementedError(
            "specifying a root group for the tree has not been implemented"
        )

    if not compute:
        raise NotImplementedError("compute=False has not been implemented yet")

    if encoding is None:
        encoding = {}

    # In the future, we may want to expand this check to insure all the provided encoding
    # options are valid. For now, this simply checks that all provided encoding keys are
    # groups in the datatree.
    if set(encoding) - set(dt.groups):
        raise ValueError(
            f"unexpected encoding group name(s) provided: {set(encoding) - set(dt.groups)}"
        )

    if unlimited_dims is None:
        unlimited_dims = {}

    for node in dt.subtree:
        ds = node.to_dataset(inherited=False)
        group_path = node.path
        if ds is None:
            _create_empty_netcdf_group(filepath, group_path, mode, engine)
        else:
            ds.to_netcdf(
                filepath,
                group=group_path,
                mode=mode,
                encoding=encoding.get(node.path),
                unlimited_dims=unlimited_dims.get(node.path),
                engine=engine,
                format=format,
                compute=compute,
                **kwargs,
            )
        mode = "a"


def _create_empty_zarr_group(
    store: MutableMapping | str | PathLike[str], group: str, mode: ZarrWriteModes
):
    import zarr

    root = zarr.open_group(store, mode=mode)
    root.create_group(group, overwrite=True)


def _datatree_to_zarr(
    dt: DataTree,
    store: MutableMapping | str | PathLike[str],
    mode: ZarrWriteModes = "w-",
    encoding: Mapping[str, Any] | None = None,
    consolidated: bool = True,
    group: str | None = None,
    compute: Literal[True] = True,
    **kwargs,
):
    """This function creates an appropriate datastore for writing a datatree
    to a zarr store.

    See `DataTree.to_zarr` for full API docs.
    """

    from zarr.convenience import consolidate_metadata

    if group is not None:
        raise NotImplementedError(
            "specifying a root group for the tree has not been implemented"
        )

    if not compute:
        raise NotImplementedError("compute=False has not been implemented yet")

    if encoding is None:
        encoding = {}

    # In the future, we may want to expand this check to insure all the provided encoding
    # options are valid. For now, this simply checks that all provided encoding keys are
    # groups in the datatree.
    if set(encoding) - set(dt.groups):
        raise ValueError(
            f"unexpected encoding group name(s) provided: {set(encoding) - set(dt.groups)}"
        )

    for node in dt.subtree:
        ds = node.to_dataset(inherited=False)
        group_path = node.path
        if ds is None:
            _create_empty_zarr_group(store, group_path, mode)
        else:
            ds.to_zarr(
                store,
                group=group_path,
                mode=mode,
                encoding=encoding.get(node.path),
                consolidated=False,
                **kwargs,
            )
        if "w" in mode:
            mode = "a"

    if consolidated:
        consolidate_metadata(store)
