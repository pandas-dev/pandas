from __future__ import annotations

import io
from collections.abc import Hashable, Mapping, MutableMapping
from os import PathLike
from typing import TYPE_CHECKING, Any, Literal, get_args

from xarray.backends.api import (
    _normalize_path,
    delayed_close_after_writes,
    dump_to_store,
    get_default_netcdf_write_engine,
    get_writable_netcdf_store,
    get_writable_zarr_store,
)
from xarray.backends.common import ArrayWriter, BytesIOProxy
from xarray.backends.locks import get_dask_scheduler
from xarray.core.datatree import DataTree
from xarray.core.types import NetcdfWriteModes, ZarrWriteModes

T_DataTreeNetcdfEngine = Literal["netcdf4", "h5netcdf", "pydap"]
T_DataTreeNetcdfTypes = Literal["NETCDF4"]

if TYPE_CHECKING:
    from dask.delayed import Delayed

    from xarray.backends import ZarrStore
    from xarray.core.types import ZarrStoreLike


def _datatree_to_netcdf(
    dt: DataTree,
    filepath: str | PathLike | io.IOBase | None = None,
    mode: NetcdfWriteModes = "w",
    encoding: Mapping[str, Any] | None = None,
    unlimited_dims: Mapping | None = None,
    format: T_DataTreeNetcdfTypes | None = None,
    engine: T_DataTreeNetcdfEngine | None = None,
    group: str | None = None,
    write_inherited_coords: bool = False,
    compute: bool = True,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> None | memoryview | Delayed:
    """Implementation of `DataTree.to_netcdf`."""

    if format not in [None, *get_args(T_DataTreeNetcdfTypes)]:
        raise ValueError("DataTree.to_netcdf only supports the NETCDF4 format")

    if engine not in [None, *get_args(T_DataTreeNetcdfEngine)]:
        raise ValueError(
            "DataTree.to_netcdf only supports the netcdf4 and h5netcdf engines"
        )

    filepath = _normalize_path(filepath)

    if engine is None:
        to_fileobject_or_memoryview = not isinstance(filepath, str)
        engine = get_default_netcdf_write_engine(
            format="NETCDF4",  # required for supporting groups
            to_fileobject_or_memoryview=to_fileobject_or_memoryview,
        )  # type: ignore[assignment]

    if group is not None:
        raise NotImplementedError(
            "specifying a root group for the tree has not been implemented"
        )

    if encoding is None:
        encoding = {}

    # In the future, we may want to expand this check to insure all the provided encoding
    # options are valid. For now, this simply checks that all provided encoding keys are
    # groups in the datatree.
    if set(encoding) - set(dt.groups):
        raise ValueError(
            f"unexpected encoding group name(s) provided: {set(encoding) - set(dt.groups)}"
        )

    if filepath is None:
        if not compute:
            raise NotImplementedError(
                "to_netcdf() with compute=False is not yet implemented when "
                "returning a memoryview"
            )
        target = BytesIOProxy()
    else:
        target = filepath  # type: ignore[assignment]

    if unlimited_dims is None:
        unlimited_dims = {}

    scheduler = get_dask_scheduler()
    have_chunks = any(
        v.chunks is not None for node in dt.subtree for v in node.variables.values()
    )
    autoclose = have_chunks and scheduler in ["distributed", "multiprocessing"]

    root_store = get_writable_netcdf_store(
        target,
        engine,  # type: ignore[arg-type]
        mode=mode,
        format=format,
        autoclose=autoclose,
        invalid_netcdf=invalid_netcdf,
        auto_complex=auto_complex,
    )

    writer = ArrayWriter()

    # TODO: allow this work (setting up the file for writing array data)
    # to be parallelized with dask
    try:
        for node in dt.subtree:
            at_root = node is dt
            dataset = node.to_dataset(inherit=write_inherited_coords or at_root)
            node_store = (
                root_store if at_root else root_store.get_child_store(node.path)
            )
            dump_to_store(
                dataset,
                node_store,
                writer,
                encoding=encoding.get(node.path),
                unlimited_dims=unlimited_dims.get(node.path),
            )

        if autoclose:
            root_store.close()

        writes = writer.sync(compute=compute)

    finally:
        if compute:
            root_store.close()
        else:
            root_store.sync()

    if filepath is None:
        assert isinstance(target, BytesIOProxy)  # created in this function
        return target.getbuffer()

    if not compute:
        return delayed_close_after_writes(writes, root_store)

    return None


def _datatree_to_zarr(
    dt: DataTree,
    store: ZarrStoreLike,
    mode: ZarrWriteModes = "w-",
    encoding: Mapping[str, Any] | None = None,
    synchronizer=None,
    group: str | None = None,
    write_inherited_coords: bool = False,
    *,
    chunk_store: MutableMapping | str | PathLike | None = None,
    compute: bool = True,
    consolidated: bool | None = None,
    append_dim: Hashable | None = None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
    safe_chunks: bool = True,
    align_chunks: bool = False,
    storage_options: dict[str, str] | None = None,
    zarr_version: int | None = None,
    zarr_format: int | None = None,
    write_empty_chunks: bool | None = None,
    chunkmanager_store_kwargs: dict[str, Any] | None = None,
) -> ZarrStore | Delayed:
    """Implementation of `DataTree.to_zarr`."""

    if group is not None:
        raise NotImplementedError(
            "specifying a root group for the tree has not been implemented"
        )

    if append_dim is not None:
        raise NotImplementedError(
            "specifying ``append_dim`` with ``DataTree.to_zarr`` has not been implemented"
        )

    if encoding is None:
        encoding = {}

    # In the future, we may want to expand this check to insure all the provided encoding
    # options are valid. For now, this simply checks that all provided encoding keys are
    # groups in the datatree.
    if set(encoding) - set(dt.groups):
        raise ValueError(
            f"unexpected encoding group name(s) provided: {set(encoding) - set(dt.groups)}"
        )

    root_store = get_writable_zarr_store(
        store,
        chunk_store=chunk_store,
        mode=mode,
        synchronizer=synchronizer,
        group=group,
        consolidated=consolidated,
        append_dim=append_dim,
        region=region,
        safe_chunks=safe_chunks,
        align_chunks=align_chunks,
        storage_options=storage_options,
        zarr_version=zarr_version,
        zarr_format=zarr_format,
        write_empty_chunks=write_empty_chunks,
    )

    writer = ArrayWriter()

    try:
        for rel_path, node in dt.subtree_with_keys:
            at_root = node is dt
            dataset = node.to_dataset(inherit=write_inherited_coords or at_root)
            # Use a relative path for group, because absolute paths are broken
            # with consolidated metadata in zarr 3.1.2 and earlier:
            # https://github.com/zarr-developers/zarr-python/pull/3428
            node_store = root_store if at_root else root_store.get_child_store(rel_path)

            dataset = node_store._validate_and_autodetect_region(dataset)
            node_store._validate_encoding(encoding)

            dump_to_store(
                dataset,
                node_store,
                writer,
                encoding=encoding.get(node.path),
            )
        writes = writer.sync(
            compute=compute, chunkmanager_store_kwargs=chunkmanager_store_kwargs
        )
    finally:
        if compute:
            root_store.close()

    if not compute:
        return delayed_close_after_writes(writes, root_store)

    return root_store
