from __future__ import annotations

import logging
import os
import time
import traceback
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from glob import glob
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Self,
    TypeVar,
    Union,
    overload,
)

import numpy as np
import pandas as pd

from xarray.coding import strings, variables
from xarray.coding.variables import SerializationWarning
from xarray.conventions import cf_encoder
from xarray.core import indexing
from xarray.core.datatree import DataTree, Variable
from xarray.core.types import ReadBuffer
from xarray.core.utils import (
    FrozenDict,
    NdimSizeLenMixin,
    attempt_import,
    emit_user_level_warning,
    is_remote_uri,
)
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array
from xarray.namedarray.utils import is_duck_dask_array

if TYPE_CHECKING:
    from xarray.core.dataset import Dataset
    from xarray.core.types import NestedSequence

    T_Name = Union[Hashable, None]

# Create a logger object, but don't add any handlers. Leave that to user code.
logger = logging.getLogger(__name__)


NONE_VAR_NAME = "__values__"

T = TypeVar("T")


@overload
def _normalize_path(path: os.PathLike) -> str: ...


@overload
def _normalize_path(path: str) -> str: ...


@overload
def _normalize_path(path: T) -> T: ...


def _normalize_path(path: os.PathLike | str | T) -> str | T:
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

    return path  # type: ignore[return-value]


@overload
def _find_absolute_paths(
    paths: str | os.PathLike | Sequence[str | os.PathLike],
    **kwargs,
) -> list[str]: ...


@overload
def _find_absolute_paths(
    paths: ReadBuffer | Sequence[ReadBuffer],
    **kwargs,
) -> list[ReadBuffer]: ...


@overload
def _find_absolute_paths(
    paths: NestedSequence[str | os.PathLike], **kwargs
) -> NestedSequence[str]: ...


@overload
def _find_absolute_paths(
    paths: NestedSequence[ReadBuffer], **kwargs
) -> NestedSequence[ReadBuffer]: ...


@overload
def _find_absolute_paths(
    paths: str
    | os.PathLike
    | ReadBuffer
    | NestedSequence[str | os.PathLike | ReadBuffer],
    **kwargs,
) -> NestedSequence[str | ReadBuffer]: ...


def _find_absolute_paths(
    paths: str
    | os.PathLike
    | ReadBuffer
    | NestedSequence[str | os.PathLike | ReadBuffer],
    **kwargs,
) -> NestedSequence[str | ReadBuffer]:
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
        if is_remote_uri(paths) and kwargs.get("engine") == "zarr":
            if TYPE_CHECKING:
                import fsspec
            else:
                fsspec = attempt_import("fsspec")

            fs, _, _ = fsspec.core.get_fs_token_paths(
                paths,
                mode="rb",
                storage_options=kwargs.get("backend_kwargs", {}).get(
                    "storage_options", {}
                ),
                expand=False,
            )
            tmp_paths = fs.glob(fs._strip_protocol(paths))  # finds directories
            return [fs.get_mapper(path) for path in tmp_paths]
        elif is_remote_uri(paths):
            raise ValueError(
                "cannot do wild-card matching for paths that are remote URLs "
                f"unless engine='zarr' is specified. Got paths: {paths}. "
                "Instead, supply paths as an explicit list of strings."
            )
        else:
            return sorted(glob(_normalize_path(paths)))
    elif isinstance(paths, os.PathLike):
        return [_normalize_path(paths)]
    elif isinstance(paths, ReadBuffer):
        return [paths]

    def _normalize_path_list(
        lpaths: NestedSequence[str | os.PathLike | ReadBuffer],
    ) -> NestedSequence[str | ReadBuffer]:
        paths = []
        for p in lpaths:
            if isinstance(p, str | os.PathLike):
                paths.append(_normalize_path(p))
            elif isinstance(p, list):
                paths.append(_normalize_path_list(p))  # type: ignore[arg-type]
            else:
                paths.append(p)  # type: ignore[arg-type]
        return paths

    return _normalize_path_list(paths)


@dataclass
class BytesIOProxy:
    """Proxy object for a write that a memoryview."""

    getvalue: Callable[[], memoryview] | None = None

    def getbuffer(self) -> memoryview:
        """Get the value of this write as bytes or memory."""
        if self.getvalue is None:
            raise ValueError("must set getvalue before fetching value")
        return self.getvalue()


def _open_remote_file(file, mode, storage_options=None):
    import fsspec

    fs, _, paths = fsspec.get_fs_token_paths(
        file, mode=mode, storage_options=storage_options
    )
    return fs.open(paths[0], mode=mode)


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
    yield str(parent)
    for path, group in root.groups.items():
        gpath = parent / path
        yield from _iter_nc_groups(group, parent=gpath)


def find_root_and_group(ds):
    """Find the root and group name of a netCDF4/h5netcdf dataset."""
    hierarchy = ()
    while ds.parent is not None:
        hierarchy = (ds.name.split("/")[-1],) + hierarchy
        ds = ds.parent
    group = "/" + "/".join(hierarchy)
    return ds, group


def collect_ancestor_dimensions(group) -> dict[str, int]:
    """Returns dimensions defined in parent groups.

    If dimensions are defined in multiple ancestors, use the size of the closest
    ancestor.
    """
    dims = {}
    while (group := group.parent) is not None:
        for k, v in group.dimensions.items():
            if k not in dims:
                dims[k] = len(v)
    return dims


def datatree_from_dict_with_io_cleanup(groups_dict: Mapping[str, Dataset]) -> DataTree:
    """DataTree.from_dict with file clean-up."""
    try:
        tree = DataTree.from_dict(groups_dict)
    except Exception:
        for ds in groups_dict.values():
            ds.close()
        raise
    for path, ds in groups_dict.items():
        tree[path].set_close(ds._close)
    return tree


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

    async def async_getitem(self, key: indexing.ExplicitIndexer) -> np.typing.ArrayLike:
        raise NotImplementedError("Backend does not support asynchronous loading")

    def get_duck_array(self, dtype: np.typing.DTypeLike | None = None):
        key = indexing.BasicIndexer((slice(None),) * self.ndim)
        return self[key]  # type: ignore[index]

    async def async_get_duck_array(self, dtype: np.typing.DTypeLike | None = None):
        key = indexing.BasicIndexer((slice(None),) * self.ndim)
        return await self.async_getitem(key)


class AbstractDataStore:
    __slots__ = ()

    def get_child_store(self, group: str) -> Self:  # pragma: no cover
        """Get a store corresponding to the indicated child group."""
        raise NotImplementedError()

    def get_dimensions(self):  # pragma: no cover
        raise NotImplementedError()

    def get_parent_dimensions(self):  # pragma: no cover
        return {}

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
                    variables = {"%s_suffix" % k: v for k, v in variables.items()}
                    attributes = {"%s_suffix" % k: v for k, v in attributes.items()}
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


T_PathFileOrDataStore = (
    str | os.PathLike[Any] | ReadBuffer | bytes | memoryview | AbstractDataStore
)


class ArrayWriter:
    __slots__ = ("lock", "regions", "sources", "targets")

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
        elif region:
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
        encoded_variables = {}
        for k, v in variables.items():
            try:
                encoded_variables[k] = self.encode_variable(v)
            except Exception as e:
                e.add_note(f"Raised while encoding variable {k!r} with value {v!r}")
                raise

        encoded_attributes = {}
        for k, v in attributes.items():
            try:
                encoded_attributes[k] = self.encode_attribute(v)
            except Exception as e:
                e.add_note(f"Raised while encoding attribute {k!r} with value {v!r}")
                raise

        return encoded_variables, encoded_attributes

    def encode_variable(self, v, name=None):
        """encode one variable"""
        return v

    def encode_attribute(self, a):
        """encode one attribute"""
        return a

    def prepare_variable(self, name, variable, check_encoding, unlimited_dims):
        raise NotImplementedError()

    def set_dimension(self, dim, length, is_unlimited):  # pragma: no cover
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

        parent_dims = self.get_parent_dimensions()
        existing_dims = self.get_dimensions()

        dims = {}
        for v in unlimited_dims:  # put unlimited_dims first
            dims[v] = None
        for v in variables.values():
            dims |= v.sizes

        for dim, length in dims.items():
            if dim in existing_dims and length != existing_dims[dim]:
                raise ValueError(
                    "Unable to update size for existing dimension"
                    f"{dim!r} ({length} != {existing_dims[dim]})"
                )
            elif dim not in existing_dims and length != parent_dims.get(dim):
                is_unlimited = dim in unlimited_dims
                self.set_dimension(dim, length, is_unlimited)

    def sync(self):
        """Write all buffered data to disk."""
        raise NotImplementedError()


def _infer_dtype(array, name=None):
    """Given an object array with no missing values, infer its dtype from all elements."""
    if array.dtype.kind != "O":
        raise TypeError("infer_type must be called on a dtype=object array")

    if array.size == 0:
        return np.dtype(float)

    native_dtypes = set(np.vectorize(type, otypes=[object])(array.ravel()))
    if len(native_dtypes) > 1 and native_dtypes != {bytes, str}:
        native_dtype_names = ", ".join(x.__name__ for x in native_dtypes)
        raise ValueError(
            f"unable to infer dtype on variable {name!r}; object array "
            f"contains mixed native types: {native_dtype_names}"
        )

    element = array[(0,) * array.ndim]
    # We use the base types to avoid subclasses of bytes and str (which might
    # not play nice with e.g. hdf5 datatypes), such as those from numpy
    if isinstance(element, bytes):
        return strings.create_vlen_dtype(bytes)
    elif isinstance(element, str):
        return strings.create_vlen_dtype(str)

    dtype = np.array(element).dtype
    if dtype.kind != "O":
        return dtype

    raise ValueError(
        f"unable to infer dtype on variable {name!r}; xarray "
        "cannot serialize arbitrary Python objects"
    )


def _copy_with_dtype(data, dtype: np.typing.DTypeLike | None):
    """Create a copy of an array with the given dtype.

    We use this instead of np.array() to ensure that custom object dtypes end
    up on the resulting array.
    """
    result = np.empty(data.shape, dtype)
    result[...] = data
    return result


def ensure_dtype_not_object(var: Variable, name: T_Name = None) -> Variable:
    if var.dtype.kind == "O":
        dims, data, attrs, encoding = variables.unpack_for_encoding(var)

        # leave vlen dtypes unchanged
        if strings.check_vlen_dtype(data.dtype) is not None:
            return var

        if is_duck_dask_array(data):
            emit_user_level_warning(
                f"variable {name} has data in the form of a dask array with "
                "dtype=object, which means it is being loaded into memory "
                "to determine a data type that can be safely stored on disk. "
                "To avoid this, coerce this variable to a fixed-size dtype "
                "with astype() before saving it.",
                category=SerializationWarning,
            )
            data = data.compute()

        missing = pd.isnull(data)
        if missing.any():
            # nb. this will fail for dask.array data
            non_missing_values = data[~missing]
            inferred_dtype = _infer_dtype(non_missing_values, name)

            # There is no safe bit-pattern for NA in typical binary string
            # formats, we so can't set a fill_value. Unfortunately, this means
            # we can't distinguish between missing values and empty strings.
            fill_value: bytes | str
            if strings.is_bytes_dtype(inferred_dtype):
                fill_value = b""
            elif strings.is_unicode_dtype(inferred_dtype):
                fill_value = ""
            else:
                # insist on using float for numeric values
                if not np.issubdtype(inferred_dtype, np.floating):
                    inferred_dtype = np.dtype(float)
                fill_value = inferred_dtype.type(np.nan)

            data = _copy_with_dtype(data, dtype=inferred_dtype)
            data[missing] = fill_value
        else:
            data = _copy_with_dtype(data, dtype=_infer_dtype(data, name))

        assert data.dtype.kind != "O" or data.dtype.metadata
        var = Variable(dims, data, attrs, encoding, fastpath=True)
    return var


class WritableCFDataStore(AbstractWritableDataStore):
    __slots__ = ()

    def encode(self, variables, attributes):
        # All NetCDF files get CF encoded by default, without this attempting
        # to write times, for example, would fail.
        variables, attributes = cf_encoder(variables, attributes)
        variables = {
            k: ensure_dtype_not_object(v, name=k) for k, v in variables.items()
        }
        return super().encode(variables, attributes)


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
    supports_groups : bool, default: False
        Whether the backend supports opening groups (via open_datatree and
        open_groups_as_dict) or not.
    """

    open_dataset_parameters: ClassVar[tuple | None] = None
    description: ClassVar[str] = ""
    url: ClassVar[str] = ""
    supports_groups: ClassVar[bool] = False

    def __repr__(self) -> str:
        txt = f"<{type(self).__name__}>"
        if self.description:
            txt += f"\n  {self.description}"
        if self.url:
            txt += f"\n  Learn more at {self.url}"
        return txt

    def open_dataset(
        self,
        filename_or_obj: str
        | os.PathLike[Any]
        | ReadBuffer
        | bytes
        | memoryview
        | AbstractDataStore,
        *,
        drop_variables: str | Iterable[str] | None = None,
    ) -> Dataset:
        """
        Backend open_dataset method used by Xarray in :py:func:`~xarray.open_dataset`.
        """

        raise NotImplementedError()

    def guess_can_open(
        self,
        filename_or_obj: str
        | os.PathLike[Any]
        | ReadBuffer
        | bytes
        | memoryview
        | AbstractDataStore,
    ) -> bool:
        """
        Backend open_dataset method used by Xarray in :py:func:`~xarray.open_dataset`.
        """

        return False

    def open_datatree(
        self,
        filename_or_obj: str
        | os.PathLike[Any]
        | ReadBuffer
        | bytes
        | memoryview
        | AbstractDataStore,
        *,
        drop_variables: str | Iterable[str] | None = None,
    ) -> DataTree:
        """
        Backend open_datatree method used by Xarray in :py:func:`~xarray.open_datatree`.

        If implemented, set the class variable supports_groups to True.
        """

        raise NotImplementedError()

    def open_groups_as_dict(
        self,
        filename_or_obj: str
        | os.PathLike[Any]
        | ReadBuffer
        | bytes
        | memoryview
        | AbstractDataStore,
        *,
        drop_variables: str | Iterable[str] | None = None,
    ) -> dict[str, Dataset]:
        """
        Opens a dictionary mapping from group names to Datasets.

        Called by :py:func:`~xarray.open_groups`.
        This function exists to provide a universal way to open all groups in a file,
        before applying any additional consistency checks or requirements necessary
        to create a `DataTree` object (typically done using :py:meth:`~xarray.DataTree.from_dict`).

        If implemented, set the class variable supports_groups to True.
        """

        raise NotImplementedError()


# mapping of engine name to (module name, BackendEntrypoint Class)
BACKEND_ENTRYPOINTS: dict[str, tuple[str | None, type[BackendEntrypoint]]] = {}


def _is_likely_dap_url(url: str) -> bool:
    """
    Determines if a URL is likely an OPeNDAP (DAP) endpoint based on
    known protocols, server software path patterns, and file extensions.

    Parameters
    ----------
    url : str

    Returns
    -------
        True if the URL matches common DAP patterns, False otherwise.
    """
    if not url:
        return False

    url_lower = url.lower()

    # For remote URIs, check for DAP server software path patterns
    if is_remote_uri(url_lower):
        dap_path_patterns = (
            "/dodsc/",  # THREDDS Data Server (TDS) DAP endpoint (case-insensitive)
            "/dods/",  # GrADS Data Server (GDS) DAP endpoint
            "/opendap/",  # Generic OPeNDAP/Hyrax server
            "/erddap/",  # ERDDAP data server
            "/dap2/",  # Explicit DAP2 version in path
            "/dap4/",  # Explicit DAP4 version in path
            "/dap/",
        )
        return any(pattern in url_lower for pattern in dap_path_patterns)

    return False
