from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

import numpy as np

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    AbstractDataStore,
    BackendArray,
    BackendEntrypoint,
    T_PathFileOrDataStore,
    _normalize_path,
    datatree_from_dict_with_io_cleanup,
    robust_getitem,
)
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core import indexing
from xarray.core.utils import (
    Frozen,
    FrozenDict,
    close_on_error,
    is_remote_uri,
)
from xarray.core.variable import Variable
from xarray.namedarray.pycompat import integer_types

if TYPE_CHECKING:
    import os

    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree
    from xarray.core.types import ReadBuffer


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
        result = robust_getitem(self.array, key, catch=ValueError)
        # in some cases, pydap doesn't squeeze axes automatically like numpy
        result = np.asarray(result)
        axis = tuple(n for n, k in enumerate(key) if isinstance(k, integer_types))
        if result.ndim + len(axis) != self.array.ndim and axis:
            result = np.squeeze(result, axis)

        return result


def get_group(ds, group):
    if group in {None, "", "/"}:
        # use the root group
        return ds
    else:
        try:
            return ds[group]
        except KeyError as e:
            # wrap error to provide slightly more helpful message
            raise KeyError(f"group not found: {group}", e) from e


class PydapDataStore(AbstractDataStore):
    """Store for accessing OpenDAP datasets with pydap.

    This store provides an alternative way to access OpenDAP datasets that may
    be useful if the netCDF4 library is not available.
    """

    def __init__(self, dataset, group=None):
        """
        Parameters
        ----------
        ds : pydap DatasetType
        group: str or None (default None)
            The group to open. If None, the root group is opened.
        """
        self.dataset = dataset
        self.group = group

    @classmethod
    def open(
        cls,
        url,
        group=None,
        application=None,
        session=None,
        output_grid=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ):
        from pydap.client import open_url
        from pydap.net import DEFAULT_TIMEOUT

        if output_grid is not None:
            # output_grid is no longer passed to pydap.client.open_url
            from xarray.core.utils import emit_user_level_warning

            emit_user_level_warning(
                "`output_grid` is deprecated and will be removed in a future version"
                " of xarray. Will be set to `None`, the new default. ",
                DeprecationWarning,
            )
            output_grid = False  # new default behavior
        kwargs = {
            "url": url,
            "application": application,
            "session": session,
            "output_grid": output_grid or False,
            "timeout": timeout or DEFAULT_TIMEOUT,
            "verify": verify or True,
            "user_charset": user_charset,
        }
        if isinstance(url, str):
            # check uit begins with an acceptable scheme
            dataset = open_url(**kwargs)
        elif hasattr(url, "ds"):
            # pydap dataset
            dataset = url.ds
        args = {"dataset": dataset}
        if group:
            # only then, change the default
            args["group"] = group
        return cls(**args)

    def open_store_variable(self, var):
        data = indexing.LazilyIndexedArray(PydapArrayWrapper(var))
        try:
            dimensions = [
                dim.split("/")[-1] if dim.startswith("/") else dim for dim in var.dims
            ]
        except AttributeError:
            # GridType does not have a dims attribute - instead get `dimensions`
            # see https://github.com/pydap/pydap/issues/485
            dimensions = var.dimensions
        return Variable(dimensions, data, var.attributes)

    def get_variables(self):
        # get first all variables arrays, excluding any container type like,
        # `Groups`, `Sequence` or `Structure` types
        try:
            _vars = list(self.ds.variables())
            _vars += list(self.ds.grids())  # dap2 objects
        except AttributeError:
            from pydap.model import GroupType

            _vars = [
                var
                for var in self.ds.keys()
                # check the key is not a BaseType or GridType
                if not isinstance(self.ds[var], GroupType)
            ]
        return FrozenDict((k, self.open_store_variable(self.ds[k])) for k in _vars)

    def get_attrs(self):
        """Remove any opendap specific attributes"""
        opendap_attrs = (
            "configuration",
            "build_dmrpp",
            "bes",
            "libdap",
            "invocation",
            "dimensions",
        )
        attrs = self.ds.attributes
        list(map(attrs.pop, opendap_attrs, [None] * 6))
        return Frozen(attrs)

    def get_dimensions(self):
        return Frozen(self.ds.dimensions)

    @property
    def ds(self):
        return get_group(self.dataset, self.group)


class PydapBackendEntrypoint(BackendEntrypoint):
    """
    Backend for steaming datasets over the internet using
    the Data Access Protocol, also known as DODS or OPeNDAP
    based on the pydap package.

    This backend is selected by default for urls.

    For more information about the underlying library, visit:
    https://pydap.github.io/pydap/en/intro.html

    See Also
    --------
    backends.PydapDataStore
    """

    description = "Open remote datasets via OPeNDAP using pydap in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.PydapBackendEntrypoint.html"

    def guess_can_open(self, filename_or_obj: T_PathFileOrDataStore) -> bool:
        return isinstance(filename_or_obj, str) and is_remote_uri(filename_or_obj)

    def open_dataset(
        self,
        filename_or_obj: (
            str | os.PathLike[Any] | ReadBuffer | bytes | memoryview | AbstractDataStore
        ),
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group=None,
        application=None,
        session=None,
        output_grid=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ) -> Dataset:
        store = PydapDataStore.open(
            url=filename_or_obj,
            group=group,
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

    def open_datatree(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group: str | None = None,
        application=None,
        session=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ) -> DataTree:
        groups_dict = self.open_groups_as_dict(
            filename_or_obj,
            mask_and_scale=mask_and_scale,
            decode_times=decode_times,
            concat_characters=concat_characters,
            decode_coords=decode_coords,
            drop_variables=drop_variables,
            use_cftime=use_cftime,
            decode_timedelta=decode_timedelta,
            group=group,
            application=None,
            session=None,
            timeout=None,
            verify=None,
            user_charset=None,
        )

        return datatree_from_dict_with_io_cleanup(groups_dict)

    def open_groups_as_dict(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group: str | None = None,
        application=None,
        session=None,
        timeout=None,
        verify=None,
        user_charset=None,
    ) -> dict[str, Dataset]:
        from xarray.core.treenode import NodePath

        filename_or_obj = _normalize_path(filename_or_obj)
        store = PydapDataStore.open(
            url=filename_or_obj,
            application=application,
            session=session,
            timeout=timeout,
            verify=verify,
            user_charset=user_charset,
        )

        # Check for a group and make it a parent if it exists
        if group:
            parent = str(NodePath("/") / NodePath(group))
        else:
            parent = str(NodePath("/"))

        groups_dict = {}
        group_names = [parent]
        # construct fully qualified path to group
        try:
            # this works for pydap >= 3.5.1
            Groups = store.ds[parent].groups()
        except AttributeError:
            # THIS IS ONLY NEEDED FOR `pydap == 3.5.0`
            # `pydap>= 3.5.1` has a new method `groups()`
            # that returns a dict of group names and their paths
            def group_fqn(store, path=None, g_fqn=None) -> dict[str, str]:
                """To be removed for pydap > 3.5.0.
                Derives the fully qualifying name of a Group."""
                from pydap.model import GroupType

                if not path:
                    path = "/"  # parent
                if not g_fqn:
                    g_fqn = {}
                groups = [
                    store[key].id
                    for key in store.keys()
                    if isinstance(store[key], GroupType)
                ]
                for g in groups:
                    g_fqn.update({g: path})
                    subgroups = [
                        var for var in store[g] if isinstance(store[g][var], GroupType)
                    ]
                    if len(subgroups) > 0:
                        npath = path + g
                        g_fqn = group_fqn(store[g], npath, g_fqn)
                return g_fqn

            Groups = group_fqn(store.ds)
        group_names += [
            str(NodePath(path_to_group) / NodePath(group))
            for group, path_to_group in Groups.items()
        ]
        for path_group in group_names:
            # get a group from the store
            store.group = path_group
            store_entrypoint = StoreBackendEntrypoint()
            with close_on_error(store):
                group_ds = store_entrypoint.open_dataset(
                    store,
                    mask_and_scale=mask_and_scale,
                    decode_times=decode_times,
                    concat_characters=concat_characters,
                    decode_coords=decode_coords,
                    drop_variables=drop_variables,
                    use_cftime=use_cftime,
                    decode_timedelta=decode_timedelta,
                )
            if group:
                group_name = str(NodePath(path_group).relative_to(parent))
            else:
                group_name = str(NodePath(path_group))
            groups_dict[group_name] = group_ds

        return groups_dict


BACKEND_ENTRYPOINTS["pydap"] = ("pydap", PydapBackendEntrypoint)
