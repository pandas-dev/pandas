from __future__ import annotations

import logging
import os
import time
import traceback
from collections.abc import Iterable
from glob import glob
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from xarray.conventions import cf_encoder
from xarray.core import indexing
from xarray.core.utils import FrozenDict, NdimSizeLenMixin, is_remote_uri
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array

if TYPE_CHECKING:
    from io import BufferedIOBase

    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree
    from xarray.core.types import NestedSequence

# Create a logger object, but don't add any handlers. Leave that to user code.
logger = logging.getLogger(__name__)


NONE_VAR_NAME = "__values__"


def _normalize_path(path):
    """
    Normalize pathlikes to string.

    Parameters
    ----------
    path :
        Path to file.

    Examples
    --------
    >>> from pathlib import Path

    >>> directory = Path(xr.backends.common.__file__).parent
    >>> paths_path = Path(directory).joinpath("comm*n.py")
    >>> paths_str = xr.backends.common._normalize_path(paths_path)
    >>> print([type(p) for p in (paths_str,)])
    [<class 'str'>]
    """
    if isinstance(path, os.PathLike):
        path = os.fspath(path)

    if isinstance(path, str) and not is_remote_uri(path):
        path = os.path.abspath(os.path.expanduser(path))

    return path


def _find_absolute_paths(
    paths: str | os.PathLike | NestedSequence[str | os.PathLike], **kwargs
) -> list[str]:
    """
    Find absolute paths from the pattern.

    Parameters
    ----------
    paths :
        Path(s) to file(s). Can include wildcards like * .
    **kwargs :
        Extra kwargs. Mainly for fsspec.

    Examples
    --------
    >>> from pathlib import Path

    >>> directory = Path(xr.backends.common.__file__).parent
    >>> paths = str(Path(directory).joinpath("comm*n.py"))  # Find common with wildcard
    >>> paths = xr.backends.common._find_absolute_paths(paths)
    >>> [Path(p).name for p in paths]
    ['common.py']
    """
    if isinstance(paths, str):
        if is_remote_uri(paths) and kwargs.get("engine", None) == "zarr":
            try:
                from fsspec.core import get_fs_token_paths
            except ImportError as e:
                raise ImportError(
                    "The use of remote URLs for opening zarr requires the package fsspec"
                ) from e

            fs, _, _ = get_fs_token_paths(
                paths,
                mode="rb",
                storage_options=kwargs.get("backend_kwargs", {}).get(
                    "storage_options", {}
                ),
                expand=False,
            )
            tmp_paths = fs.glob(fs._strip_protocol(paths))  # finds directories
            paths = [fs.get_mapper(path) for path in tmp_paths]
        elif is_remote_uri(paths):
            raise ValueError(
                "cannot do wild-card matching for paths that are remote URLs "
                f"unless engine='zarr' is specified. Got paths: {paths}. "
                "Instead, supply paths as an explicit list of strings."
            )
        else:
            paths = sorted(glob(_normalize_path(paths)))
    elif isinstance(paths, os.PathLike):
        paths = [os.fspath(paths)]
    else:
        paths = [os.fspath(p) if isinstance(p, os.PathLike) else p for p in paths]

    return paths


def _encode_variable_name(name):
    if name is None:
        name = NONE_VAR_NAME
    return name


def _decode_variable_name(name):
    if name == NONE_VAR_NAME:
        name = None
    return name


def _iter_nc_groups(root, parent="/"):
    from xarray.core.treenode import NodePath

    parent = NodePath(parent)
    for path, group in root.groups.items():
        gpath = parent / path
        yield str(gpath)
        yield from _iter_nc_groups(group, parent=gpath)


def find_root_and_group(ds):
    """Find the root and group name of a netCDF4/h5netcdf dataset."""
    hierarchy = ()
    while ds.parent is not None:
        hierarchy = (ds.name.split("/")[-1],) + hierarchy
        ds = ds.parent
    group = "/" + "/".join(hierarchy)
    return ds, group


def robust_getitem(array, key, catch=Exception, max_retries=6, initial_delay=500):
    """
    Robustly index an array, using retry logic with exponential backoff if any
    of the errors ``catch`` are raised. The initial_delay is measured in ms.

    With the default settings, the maximum delay will be in the range of 32-64
    seconds.
    """
    assert max_retries >= 0
    for n in range(max_retries + 1):
        try:
            return array[key]
        except catch:
            if n == max_retries:
                raise
            base_delay = initial_delay * 2**n
            next_delay = base_delay + np.random.randint(base_delay)
            msg = (
                f"getitem failed, waiting {next_delay} ms before trying again "
                f"({max_retries - n} tries remaining). Full traceback: {traceback.format_exc()}"
            )
            logger.debug(msg)
            time.sleep(1e-3 * next_delay)


class BackendArray(NdimSizeLenMixin, indexing.ExplicitlyIndexed):
    __slots__ = ()

    def get_duck_array(self, dtype: np.typing.DTypeLike = None):
        key = indexing.BasicIndexer((slice(None),) * self.ndim)
        return self[key]  # type: ignore [index]


class AbstractDataStore:
    __slots__ = ()

    def get_dimensions(self):  # pragma: no cover
        raise NotImplementedError()

    def get_attrs(self):  # pragma: no cover
        raise NotImplementedError()

    def get_variables(self):  # pragma: no cover
        raise NotImplementedError()

    def get_encoding(self):
        return {}

    def load(self):
        """
        This loads the variables and attributes simultaneously.
        A centralized loading function makes it easier to create
        data stores that do automatic encoding/decoding.

        For example::

            class SuffixAppendingDataStore(AbstractDataStore):

                def load(self):
                    variables, attributes = AbstractDataStore.load(self)
                    variables = {'%s_suffix' % k: v
                                 for k, v in variables.items()}
                    attributes = {'%s_suffix' % k: v
                                  for k, v in attributes.items()}
                    return variables, attributes

        This function will be called anytime variables or attributes
        are requested, so care should be taken to make sure its fast.
        """
        variables = FrozenDict(
            (_decode_variable_name(k), v) for k, v in self.get_variables().items()
        )
        attributes = FrozenDict(self.get_attrs())
        return variables, attributes

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()


class ArrayWriter:
    __slots__ = ("sources", "targets", "regions", "lock")

    def __init__(self, lock=None):
        self.sources = []
        self.targets = []
        self.regions = []
        self.lock = lock

    def add(self, source, target, region=None):
        if is_chunked_array(source):
            self.sources.append(source)
            self.targets.append(target)
            self.regions.append(region)
        else:
            if region:
                target[region] = source
            else:
                target[...] = source

    def sync(self, compute=True, chunkmanager_store_kwargs=None):
        if self.sources:
            chunkmanager = get_chunked_array_type(*self.sources)

            # TODO: consider wrapping targets with dask.delayed, if this makes
            # for any discernible difference in performance, e.g.,
            # targets = [dask.delayed(t) for t in self.targets]

            if chunkmanager_store_kwargs is None:
                chunkmanager_store_kwargs = {}

            delayed_store = chunkmanager.store(
                self.sources,
                self.targets,
                lock=self.lock,
                compute=compute,
                flush=True,
                regions=self.regions,
                **chunkmanager_store_kwargs,
            )
            self.sources = []
            self.targets = []
            self.regions = []
            return delayed_store


class AbstractWritableDataStore(AbstractDataStore):
    __slots__ = ()

    def encode(self, variables, attributes):
        """
        Encode the variables and attributes in this store

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        attributes : dict-like
            Dictionary of key/value (attribute name / attribute) pairs

        Returns
        -------
        variables : dict-like
        attributes : dict-like

        """
        variables = {k: self.encode_variable(v) for k, v in variables.items()}
        attributes = {k: self.encode_attribute(v) for k, v in attributes.items()}
        return variables, attributes

    def encode_variable(self, v):
        """encode one variable"""
        return v

    def encode_attribute(self, a):
        """encode one attribute"""
        return a

    def set_dimension(self, dim, length):  # pragma: no cover
        raise NotImplementedError()

    def set_attribute(self, k, v):  # pragma: no cover
        raise NotImplementedError()

    def set_variable(self, k, v):  # pragma: no cover
        raise NotImplementedError()

    def store_dataset(self, dataset):
        """
        in stores, variables are all variables AND coordinates
        in xarray.Dataset variables are variables NOT coordinates,
        so here we pass the whole dataset in instead of doing
        dataset.variables
        """
        self.store(dataset, dataset.attrs)

    def store(
        self,
        variables,
        attributes,
        check_encoding_set=frozenset(),
        writer=None,
        unlimited_dims=None,
    ):
        """
        Top level method for putting data on this store, this method:
          - encodes variables/attributes
          - sets dimensions
          - sets variables

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        attributes : dict-like
            Dictionary of key/value (attribute name / attribute) pairs
        check_encoding_set : list-like
            List of variables that should be checked for invalid encoding
            values
        writer : ArrayWriter
        unlimited_dims : list-like
            List of dimension names that should be treated as unlimited
            dimensions.
        """
        if writer is None:
            writer = ArrayWriter()

        variables, attributes = self.encode(variables, attributes)

        self.set_attributes(attributes)
        self.set_dimensions(variables, unlimited_dims=unlimited_dims)
        self.set_variables(
            variables, check_encoding_set, writer, unlimited_dims=unlimited_dims
        )

    def set_attributes(self, attributes):
        """
        This provides a centralized method to set the dataset attributes on the
        data store.

        Parameters
        ----------
        attributes : dict-like
            Dictionary of key/value (attribute name / attribute) pairs
        """
        for k, v in attributes.items():
            self.set_attribute(k, v)

    def set_variables(self, variables, check_encoding_set, writer, unlimited_dims=None):
        """
        This provides a centralized method to set the variables on the data
        store.

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        check_encoding_set : list-like
            List of variables that should be checked for invalid encoding
            values
        writer : ArrayWriter
        unlimited_dims : list-like
            List of dimension names that should be treated as unlimited
            dimensions.
        """

        for vn, v in variables.items():
            name = _encode_variable_name(vn)
            check = vn in check_encoding_set
            target, source = self.prepare_variable(
                name, v, check, unlimited_dims=unlimited_dims
            )

            writer.add(source, target)

    def set_dimensions(self, variables, unlimited_dims=None):
        """
        This provides a centralized method to set the dimensions on the data
        store.

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        unlimited_dims : list-like
            List of dimension names that should be treated as unlimited
            dimensions.
        """
        if unlimited_dims is None:
            unlimited_dims = set()

        existing_dims = self.get_dimensions()

        dims = {}
        for v in unlimited_dims:  # put unlimited_dims first
            dims[v] = None
        for v in variables.values():
            dims.update(dict(zip(v.dims, v.shape)))

        for dim, length in dims.items():
            if dim in existing_dims and length != existing_dims[dim]:
                raise ValueError(
                    "Unable to update size for existing dimension"
                    f"{dim!r} ({length} != {existing_dims[dim]})"
                )
            elif dim not in existing_dims:
                is_unlimited = dim in unlimited_dims
                self.set_dimension(dim, length, is_unlimited)


class WritableCFDataStore(AbstractWritableDataStore):
    __slots__ = ()

    def encode(self, variables, attributes):
        # All NetCDF files get CF encoded by default, without this attempting
        # to write times, for example, would fail.
        variables, attributes = cf_encoder(variables, attributes)
        variables = {k: self.encode_variable(v) for k, v in variables.items()}
        attributes = {k: self.encode_attribute(v) for k, v in attributes.items()}
        return variables, attributes


class BackendEntrypoint:
    """
    ``BackendEntrypoint`` is a class container and it is the main interface
    for the backend plugins, see :ref:`RST backend_entrypoint`.
    It shall implement:

    - ``open_dataset`` method: it shall implement reading from file, variables
      decoding and it returns an instance of :py:class:`~xarray.Dataset`.
      It shall take in input at least ``filename_or_obj`` argument and
      ``drop_variables`` keyword argument.
      For more details see :ref:`RST open_dataset`.
    - ``guess_can_open`` method: it shall return ``True`` if the backend is able to open
      ``filename_or_obj``, ``False`` otherwise. The implementation of this
      method is not mandatory.
    - ``open_datatree`` method: it shall implement reading from file, variables
      decoding and it returns an instance of :py:class:`~datatree.DataTree`.
      It shall take in input at least ``filename_or_obj`` argument. The
      implementation of this method is not mandatory.  For more details see
      <reference to open_datatree documentation>.

    Attributes
    ----------

    open_dataset_parameters : tuple, default: None
        A list of ``open_dataset`` method parameters.
        The setting of this attribute is not mandatory.
    description : str, default: ""
        A short string describing the engine.
        The setting of this attribute is not mandatory.
    url : str, default: ""
        A string with the URL to the backend's documentation.
        The setting of this attribute is not mandatory.
    """

    open_dataset_parameters: ClassVar[tuple | None] = None
    description: ClassVar[str] = ""
    url: ClassVar[str] = ""

    def __repr__(self) -> str:
        txt = f"<{type(self).__name__}>"
        if self.description:
            txt += f"\n  {self.description}"
        if self.url:
            txt += f"\n  Learn more at {self.url}"
        return txt

    def open_dataset(
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
        *,
        drop_variables: str | Iterable[str] | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Backend open_dataset method used by Xarray in :py:func:`~xarray.open_dataset`.
        """

        raise NotImplementedError()

    def guess_can_open(
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
    ) -> bool:
        """
        Backend open_dataset method used by Xarray in :py:func:`~xarray.open_dataset`.
        """

        return False

    def open_datatree(
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
        **kwargs: Any,
    ) -> DataTree:
        """
        Backend open_datatree method used by Xarray in :py:func:`~xarray.open_datatree`.
        """

        raise NotImplementedError()


# mapping of engine name to (module name, BackendEntrypoint Class)
BACKEND_ENTRYPOINTS: dict[str, tuple[str | None, type[BackendEntrypoint]]] = {}
