from __future__ import annotations

import functools
import io
import itertools
import textwrap
from collections import ChainMap
from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Iterator,
    Mapping,
)
from html import escape
from os import PathLike
from typing import (
    TYPE_CHECKING,
    Any,
    Concatenate,
    Literal,
    NoReturn,
    ParamSpec,
    TypeVar,
    Union,
    overload,
)

from xarray.core import utils
from xarray.core._aggregations import DataTreeAggregations
from xarray.core._typed_ops import DataTreeOpsMixin
from xarray.core.common import TreeAttrAccessMixin, get_chunksizes
from xarray.core.coordinates import Coordinates, DataTreeCoordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.dataset_variables import DataVariables
from xarray.core.datatree_mapping import (
    _handle_errors_with_path_context,
    map_over_datasets,
)
from xarray.core.formatting import (
    datatree_repr,
    diff_treestructure,
    dims_and_coords_repr,
)
from xarray.core.formatting_html import (
    datatree_repr as datatree_repr_html,
)
from xarray.core.indexes import Index, Indexes
from xarray.core.options import OPTIONS as XR_OPTS
from xarray.core.options import _get_keep_attrs
from xarray.core.treenode import NamedNode, NodePath, zip_subtrees
from xarray.core.types import Self
from xarray.core.utils import (
    Default,
    FilteredMapping,
    Frozen,
    _default,
    drop_dims_from_indexers,
    either_dict_or_kwargs,
    maybe_wrap_array,
    parse_dims_as_set,
)
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array
from xarray.structure.alignment import align
from xarray.structure.merge import dataset_update_method

try:
    from xarray.core.variable import calculate_dimensions
except ImportError:
    # for xarray versions 2022.03.0 and earlier
    from xarray.core.dataset import calculate_dimensions

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd
    from dask.delayed import Delayed

    from xarray.backends import ZarrStore
    from xarray.core.datatree_io import T_DataTreeNetcdfEngine, T_DataTreeNetcdfTypes
    from xarray.core.types import (
        Dims,
        DtCompatible,
        ErrorOptions,
        ErrorOptionsWithWarn,
        NetcdfWriteModes,
        T_ChunkDimFreq,
        T_ChunksFreq,
        ZarrWriteModes,
    )
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint
    from xarray.structure.merge import CoercibleMapping, CoercibleValue

# """
# DEVELOPERS' NOTE
# ----------------
# The idea of this module is to create a `DataTree` class which inherits the tree
# structure from TreeNode, and also copies the entire API of `xarray.Dataset`, but with
# certain methods decorated to instead map the dataset function over every node in the
# tree. As this API is copied without directly subclassing `xarray.Dataset` we instead
# create various Mixin classes (in ops.py) which each define part of `xarray.Dataset`'s
# extensive API.
#
# Some of these methods must be wrapped to map over all nodes in the subtree. Others are
# fine to inherit unaltered (normally because they (a) only call dataset properties and
# (b) don't return a dataset that should be nested into a new tree) and some will get
# overridden by the class definition of DataTree.
# """


T_Path = Union[str, NodePath]
T = TypeVar("T")
P = ParamSpec("P")


def _collect_data_and_coord_variables(
    data: Dataset,
) -> tuple[dict[Hashable, Variable], dict[Hashable, Variable]]:
    data_variables = {}
    coord_variables = {}
    for k, v in data.variables.items():
        if k in data._coord_names:
            coord_variables[k] = v
        else:
            data_variables[k] = v
    return data_variables, coord_variables


def _to_new_dataset(data: Dataset | Coordinates | None) -> Dataset:
    if isinstance(data, Dataset):
        ds = data.copy(deep=False)
    elif isinstance(data, Coordinates):
        ds = data.to_dataset()
    elif data is None:
        ds = Dataset()
    else:
        raise TypeError(f"data object is not an xarray.Dataset, dict, or None: {data}")
    return ds


def _inherited_dataset(ds: Dataset, parent: Dataset) -> Dataset:
    return Dataset._construct_direct(
        variables=parent._variables | ds._variables,
        coord_names=parent._coord_names | ds._coord_names,
        dims=parent._dims | ds._dims,
        attrs=ds._attrs,
        indexes=parent._indexes | ds._indexes,
        encoding=ds._encoding,
        close=ds._close,
    )


def _without_header(text: str) -> str:
    return "\n".join(text.split("\n")[1:])


def _indented(text: str) -> str:
    return textwrap.indent(text, prefix="    ")


def check_alignment(
    path: str,
    node_ds: Dataset,
    parent_ds: Dataset | None,
    children: Mapping[str, DataTree],
) -> None:
    if parent_ds is not None:
        try:
            align(node_ds, parent_ds, join="exact", copy=False)
        except ValueError as e:
            node_repr = _indented(_without_header(repr(node_ds)))
            parent_repr = _indented(dims_and_coords_repr(parent_ds))
            raise ValueError(
                f"group {path!r} is not aligned with its parents:\n"
                f"Group:\n{node_repr}\nFrom parents:\n{parent_repr}"
            ) from e

    if children:
        if parent_ds is not None:
            base_ds = _inherited_dataset(node_ds, parent_ds)
        else:
            base_ds = node_ds

        for child_name, child in children.items():
            child_path = str(NodePath(path) / child_name)
            child_ds = child.to_dataset(inherit=False)
            check_alignment(child_path, child_ds, base_ds, child.children)


def _deduplicate_inherited_coordinates(child: DataTree, parent: DataTree) -> None:
    # This method removes repeated indexes (and corresponding coordinates)
    # that are repeated between a DataTree and its parents.
    removed_something = False
    for name in parent._indexes:
        if name in child._node_indexes:
            # Indexes on a Dataset always have a corresponding coordinate.
            # We already verified that these coordinates match in the
            # check_alignment() call from _pre_attach().
            del child._node_indexes[name]
            del child._node_coord_variables[name]
            removed_something = True

    if removed_something:
        child._node_dims = calculate_dimensions(
            child._data_variables | child._node_coord_variables
        )

    for grandchild in child._children.values():
        _deduplicate_inherited_coordinates(grandchild, child)


def _check_for_slashes_in_names(variables: Iterable[Hashable]) -> None:
    offending_variable_names = [
        name for name in variables if isinstance(name, str) and "/" in name
    ]
    if len(offending_variable_names) > 0:
        raise ValueError(
            "Given variables have names containing the '/' character: "
            f"{offending_variable_names}. "
            "Variables stored in DataTree objects cannot have names containing '/' characters, as this would make path-like access to variables ambiguous."
        )


class DatasetView(Dataset):
    """
    An immutable Dataset-like view onto the data in a single DataTree node.

    In-place operations modifying this object should raise an AttributeError.
    This requires overriding all inherited constructors.

    Operations returning a new result will return a new xarray.Dataset object.
    This includes all API on Dataset, which will be inherited.
    """

    # TODO what happens if user alters (in-place) a DataArray they extracted from this object?

    __slots__ = (
        "_attrs",
        "_cache",  # used by _CachedAccessor
        "_close",
        "_coord_names",
        "_dims",
        "_encoding",
        "_indexes",
        "_variables",
    )

    def __init__(
        self,
        data_vars: Mapping[Any, Any] | None = None,
        coords: Mapping[Any, Any] | None = None,
        attrs: Mapping[Any, Any] | None = None,
    ):
        raise AttributeError("DatasetView objects are not to be initialized directly")

    @classmethod
    def _constructor(
        cls,
        variables: dict[Any, Variable],
        coord_names: set[Hashable],
        dims: dict[Any, int],
        attrs: dict | None,
        indexes: dict[Any, Index],
        encoding: dict | None,
        close: Callable[[], None] | None,
    ) -> DatasetView:
        """Private constructor, from Dataset attributes."""
        # We override Dataset._construct_direct below, so we need a new
        # constructor for creating DatasetView objects.
        obj: DatasetView = object.__new__(cls)
        obj._variables = variables
        obj._coord_names = coord_names
        obj._dims = dims
        obj._indexes = indexes
        obj._attrs = attrs
        obj._close = close
        obj._encoding = encoding
        return obj

    def __setitem__(self, key, val) -> None:
        raise AttributeError(
            "Mutation of the DatasetView is not allowed, please use `.__setitem__` on the wrapping DataTree node, "
            "or use `dt.to_dataset()` if you want a mutable dataset. If calling this from within `map_over_datasets`,"
            "use `.copy()` first to get a mutable version of the input dataset."
        )

    def update(self, other) -> NoReturn:
        raise AttributeError(
            "Mutation of the DatasetView is not allowed, please use `.update` on the wrapping DataTree node, "
            "or use `dt.to_dataset()` if you want a mutable dataset. If calling this from within `map_over_datasets`,"
            "use `.copy()` first to get a mutable version of the input dataset."
        )

    def set_close(self, close: Callable[[], None] | None) -> None:
        raise AttributeError("cannot modify a DatasetView()")

    def close(self) -> None:
        raise AttributeError(
            "cannot close a DatasetView(). Close the associated DataTree node instead"
        )

    # FIXME https://github.com/python/mypy/issues/7328
    @overload  # type: ignore[override]
    def __getitem__(self, key: Mapping) -> Dataset:  # type: ignore[overload-overlap]
        ...

    @overload
    def __getitem__(self, key: Hashable) -> DataArray: ...

    # See: https://github.com/pydata/xarray/issues/8855
    @overload
    def __getitem__(self, key: Any) -> Dataset: ...

    def __getitem__(self, key) -> DataArray | Dataset:
        # TODO call the `_get_item` method of DataTree to allow path-like access to contents of other nodes
        # For now just call Dataset.__getitem__
        return Dataset.__getitem__(self, key)

    @classmethod
    def _construct_direct(  # type: ignore[override]
        cls,
        variables: dict[Any, Variable],
        coord_names: set[Hashable],
        dims: dict[Any, int] | None = None,
        attrs: dict | None = None,
        indexes: dict[Any, Index] | None = None,
        encoding: dict | None = None,
        close: Callable[[], None] | None = None,
    ) -> Dataset:
        """
        Overriding this method (along with ._replace) and modifying it to return a Dataset object
        should hopefully ensure that the return type of any method on this object is a Dataset.
        """
        if dims is None:
            dims = calculate_dimensions(variables)
        if indexes is None:
            indexes = {}
        obj = object.__new__(Dataset)
        obj._variables = variables
        obj._coord_names = coord_names
        obj._dims = dims
        obj._indexes = indexes
        obj._attrs = attrs
        obj._close = close
        obj._encoding = encoding
        return obj

    def _replace(  # type: ignore[override]
        self,
        variables: dict[Hashable, Variable] | None = None,
        coord_names: set[Hashable] | None = None,
        dims: dict[Any, int] | None = None,
        attrs: dict[Hashable, Any] | Default | None = _default,
        indexes: dict[Hashable, Index] | None = None,
        encoding: dict | Default | None = _default,
        inplace: bool = False,
    ) -> Dataset:
        """
        Overriding this method (along with ._construct_direct) and modifying it to return a Dataset object
        should hopefully ensure that the return type of any method on this object is a Dataset.
        """

        if inplace:
            raise AttributeError("In-place mutation of the DatasetView is not allowed")

        return Dataset._replace(
            self,
            variables=variables,
            coord_names=coord_names,
            dims=dims,
            attrs=attrs,
            indexes=indexes,
            encoding=encoding,
            inplace=inplace,
        )

    def map(  # type: ignore[override]
        self,
        func: Callable,
        keep_attrs: bool | None = None,
        args: Iterable[Any] = (),
        **kwargs: Any,
    ) -> Dataset:
        """Apply a function to each data variable in this dataset

        Parameters
        ----------
        func : callable
            Function which can be called in the form `func(x, *args, **kwargs)`
            to transform each DataArray `x` in this dataset into another
            DataArray.
        keep_attrs : bool | None, optional
            If True, both the dataset's and variables' attributes (`attrs`) will be
            copied from the original objects to the new ones. If False, the new dataset
            and variables will be returned without copying the attributes.
        args : iterable, optional
            Positional arguments passed on to `func`.
        **kwargs : Any
            Keyword arguments passed on to `func`.

        Returns
        -------
        applied : Dataset
            Resulting dataset from applying ``func`` to each data variable.

        Examples
        --------
        >>> da = xr.DataArray(np.random.randn(2, 3))
        >>> ds = xr.Dataset({"foo": da, "bar": ("x", [-1, 2])})
        >>> ds
        <xarray.Dataset> Size: 64B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Dimensions without coordinates: dim_0, dim_1, x
        Data variables:
            foo      (dim_0, dim_1) float64 48B 1.764 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2
        >>> ds.map(np.fabs)
        <xarray.Dataset> Size: 64B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Dimensions without coordinates: dim_0, dim_1, x
        Data variables:
            foo      (dim_0, dim_1) float64 48B 1.764 0.4002 0.9787 2.241 1.868 0.9773
            bar      (x) float64 16B 1.0 2.0
        """

        # Copied from xarray.Dataset so as not to call type(self), which causes problems (see https://github.com/xarray-contrib/datatree/issues/188).
        # TODO Refactor xarray upstream to avoid needing to overwrite this.
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)
        variables = {
            k: maybe_wrap_array(v, func(v, *args, **kwargs))
            for k, v in self.data_vars.items()
        }
        if keep_attrs:
            for k, v in variables.items():
                v._copy_attrs_from(self.data_vars[k])
        attrs = self.attrs if keep_attrs else None
        # return type(self)(variables, attrs=attrs)
        return Dataset(variables, attrs=attrs)


class DataTree(
    NamedNode,
    DataTreeAggregations,
    DataTreeOpsMixin,
    TreeAttrAccessMixin,
    Mapping[str, "DataArray | DataTree"],
):
    """
    A tree-like hierarchical collection of xarray objects.

    Attempts to present an API like that of xarray.Dataset, but methods are wrapped to also update all the tree's child nodes.
    """

    # TODO Some way of sorting children by depth

    # TODO do we need a watch out for if methods intended only for root nodes are called on non-root nodes?

    # TODO dataset methods which should not or cannot act over the whole tree, such as .to_array

    # TODO .loc method

    # TODO a lot of properties like .variables could be defined in a DataMapping class which both Dataset and DataTree inherit from

    # TODO all groupby classes

    # TODO a lot of properties like .variables could be defined in a DataMapping class which both Dataset and DataTree inherit from

    # TODO all groupby classes

    _name: str | None
    _parent: DataTree | None
    _children: dict[str, DataTree]
    _cache: dict[str, Any]  # used by _CachedAccessor
    _data_variables: dict[Hashable, Variable]
    _node_coord_variables: dict[Hashable, Variable]
    _node_dims: dict[Hashable, int]
    _node_indexes: dict[Hashable, Index]
    _attrs: dict[Hashable, Any] | None
    _encoding: dict[Hashable, Any] | None
    _close: Callable[[], None] | None

    __slots__ = (
        "_attrs",
        "_cache",  # used by _CachedAccessor
        "_children",
        "_close",
        "_data_variables",
        "_encoding",
        "_name",
        "_node_coord_variables",
        "_node_dims",
        "_node_indexes",
        "_parent",
    )

    def __init__(
        self,
        dataset: Dataset | Coordinates | None = None,
        children: Mapping[str, DataTree] | None = None,
        name: str | None = None,
    ):
        """
        Create a single node of a DataTree.

        The node may optionally contain data in the form of data and coordinate
        variables, stored in the same way as data is stored in an
        xarray.Dataset.

        Parameters
        ----------
        dataset : Dataset, optional
            Data to store directly at this node.
        children : Mapping[str, DataTree], optional
            Any child nodes of this node.
        name : str, optional
            Name for this node of the tree.

        Returns
        -------
        DataTree

        See Also
        --------
        DataTree.from_dict
        """
        self._set_node_data(_to_new_dataset(dataset))

        # comes after setting node data as this will check for clashes between child names and existing variable names
        super().__init__(name=name, children=children)

    def _set_node_data(self, dataset: Dataset):
        _check_for_slashes_in_names(dataset.variables)
        data_vars, coord_vars = _collect_data_and_coord_variables(dataset)
        self._data_variables = data_vars
        self._node_coord_variables = coord_vars
        self._node_dims = dataset._dims
        self._node_indexes = dataset._indexes
        self._encoding = dataset._encoding
        self._attrs = dataset._attrs
        self._close = dataset._close

    def _pre_attach(self: DataTree, parent: DataTree, name: str) -> None:
        super()._pre_attach(parent, name)
        if name in parent.dataset.variables:
            raise KeyError(
                f"parent {parent.name} already contains a variable named {name}"
            )
        path = str(NodePath(parent.path) / name)
        node_ds = self.to_dataset(inherit=False)
        parent_ds = parent._to_dataset_view(rebuild_dims=False, inherit=True)
        check_alignment(path, node_ds, parent_ds, self.children)
        _deduplicate_inherited_coordinates(self, parent)

    @property
    def _node_coord_variables_with_index(self) -> Mapping[Hashable, Variable]:
        return FilteredMapping(
            keys=self._node_indexes, mapping=self._node_coord_variables
        )

    @property
    def _coord_variables(self) -> ChainMap[Hashable, Variable]:
        # ChainMap is incorrected typed in typeshed (only the first argument
        # needs to be mutable)
        # https://github.com/python/typeshed/issues/8430
        return ChainMap(
            self._node_coord_variables,
            *(p._node_coord_variables_with_index for p in self.parents),  # type: ignore[arg-type]
        )

    @property
    def _dims(self) -> ChainMap[Hashable, int]:
        return ChainMap(self._node_dims, *(p._node_dims for p in self.parents))

    @property
    def _indexes(self) -> ChainMap[Hashable, Index]:
        return ChainMap(self._node_indexes, *(p._node_indexes for p in self.parents))

    def _to_dataset_view(self, rebuild_dims: bool, inherit: bool) -> DatasetView:
        coord_vars = self._coord_variables if inherit else self._node_coord_variables
        variables = dict(self._data_variables)
        variables |= coord_vars
        if rebuild_dims:
            dims = calculate_dimensions(variables)
        elif inherit:
            # Note: rebuild_dims=False with inherit=True can create
            # technically invalid Dataset objects because it still includes
            # dimensions that are only defined on parent data variables
            # (i.e. not present on any parent coordinate variables).
            #
            # For example:
            #     >>> tree = DataTree.from_dict(
            #     ...     {
            #     ...         "/": xr.Dataset({"foo": ("x", [1, 2])}),  # x has size 2
            #     ...         "/b": xr.Dataset(),
            #     ...     }
            #     ... )
            #     >>> ds = tree["b"]._to_dataset_view(rebuild_dims=False, inherit=True)
            #     >>> ds
            #     <xarray.DatasetView> Size: 0B
            #     Dimensions:  (x: 2)
            #     Dimensions without coordinates: x
            #     Data variables:
            #         *empty*
            #
            # Notice the "x" dimension is still defined, even though there are no variables
            # or coordinates.
            #
            # Normally this is not supposed to be possible in xarray's data model,
            # but here it is useful internally for use cases where we
            # want to inherit everything from parents nodes, e.g., for align() and repr().
            #
            # The user should never be able to see this dimension via public API.
            dims = dict(self._dims)
        else:
            dims = dict(self._node_dims)
        return DatasetView._constructor(
            variables=variables,
            coord_names=set(self._coord_variables),
            dims=dims,
            attrs=self._attrs,
            indexes=dict(self._indexes if inherit else self._node_indexes),
            encoding=self._encoding,
            close=None,
        )

    @property
    def dataset(self) -> DatasetView:
        """
        An immutable Dataset-like view onto the data in this node.

        Includes inherited coordinates and indexes from parent nodes.

        For a mutable Dataset containing the same data as in this node, use
        `.to_dataset()` instead.

        See Also
        --------
        DataTree.to_dataset
        """
        return self._to_dataset_view(rebuild_dims=True, inherit=True)

    @dataset.setter
    def dataset(self, data: Dataset | None = None) -> None:
        ds = _to_new_dataset(data)
        self._replace_node(ds)

    # soft-deprecated alias, to facilitate the transition from
    # xarray-contrib/datatree
    ds = dataset

    def to_dataset(self, inherit: bool = True) -> Dataset:
        """
        Return the data in this node as a new xarray.Dataset object.

        Parameters
        ----------
        inherit : bool, optional
            If False, only include coordinates and indexes defined at the level
            of this DataTree node, excluding any inherited coordinates and indexes.

        See Also
        --------
        DataTree.dataset
        """
        coord_vars = self._coord_variables if inherit else self._node_coord_variables
        variables = dict(self._data_variables)
        variables |= coord_vars
        dims = calculate_dimensions(variables) if inherit else dict(self._node_dims)
        return Dataset._construct_direct(
            variables,
            set(coord_vars),
            dims,
            None if self._attrs is None else dict(self._attrs),
            dict(self._indexes if inherit else self._node_indexes),
            None if self._encoding is None else dict(self._encoding),
            None,
        )

    @property
    def has_data(self) -> bool:
        """Whether or not there are any variables in this node."""
        return bool(self._data_variables or self._node_coord_variables)

    @property
    def has_attrs(self) -> bool:
        """Whether or not there are any metadata attributes in this node."""
        return len(self.attrs.keys()) > 0

    @property
    def is_empty(self) -> bool:
        """False if node contains any data or attrs. Does not look at children."""
        return not (self.has_data or self.has_attrs)

    @property
    def is_hollow(self) -> bool:
        """True if only leaf nodes contain data."""
        return not any(node.has_data for node in self.subtree if not node.is_leaf)

    @property
    def variables(self) -> Mapping[Hashable, Variable]:
        """Low level interface to node contents as dict of Variable objects.

        This dictionary is frozen to prevent mutation that could violate
        Dataset invariants. It contains all variable objects constituting this
        DataTree node, including both data variables and coordinates.
        """
        return Frozen(self._data_variables | self._coord_variables)

    @property
    def attrs(self) -> dict[Hashable, Any]:
        """Dictionary of global attributes on this node object."""
        if self._attrs is None:
            self._attrs = {}
        return self._attrs

    @attrs.setter
    def attrs(self, value: Mapping[Any, Any]) -> None:
        self._attrs = dict(value)

    @property
    def encoding(self) -> dict:
        """Dictionary of global encoding attributes on this node object."""
        if self._encoding is None:
            self._encoding = {}
        return self._encoding

    @encoding.setter
    def encoding(self, value: Mapping) -> None:
        self._encoding = dict(value)

    @property
    def dims(self) -> Mapping[Hashable, int]:
        """Mapping from dimension names to lengths.

        Cannot be modified directly, but is updated when adding new variables.

        Note that type of this object differs from `DataArray.dims`.
        See `DataTree.sizes`, `Dataset.sizes`, and `DataArray.sizes` for consistently named
        properties.
        """
        return Frozen(self._dims)

    @property
    def sizes(self) -> Mapping[Hashable, int]:
        """Mapping from dimension names to lengths.

        Cannot be modified directly, but is updated when adding new variables.

        This is an alias for `DataTree.dims` provided for the benefit of
        consistency with `DataArray.sizes`.

        See Also
        --------
        DataArray.sizes
        """
        return self.dims

    @property
    def _attr_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for attribute-style access"""
        yield from self._item_sources
        yield self.attrs

    @property
    def _item_sources(self) -> Iterable[Mapping[Any, Any]]:
        """Places to look-up items for key-completion"""
        yield self.data_vars
        yield FilteredMapping(keys=self._coord_variables, mapping=self.coords)

        # virtual coordinates
        yield FilteredMapping(keys=self.dims, mapping=self)

        # immediate child nodes
        yield self.children

    def _ipython_key_completions_(self) -> list[str]:
        """Provide method for the key-autocompletions in IPython.
        See https://ipython.readthedocs.io/en/stable/config/integrating.html#tab-completion
        For the details.
        """

        # TODO allow auto-completing relative string paths, e.g. `dt['path/to/../ <tab> node'`
        # Would require changes to ipython's autocompleter, see https://github.com/ipython/ipython/issues/12420
        # Instead for now we only list direct paths to all node in subtree explicitly

        items_on_this_node = self._item_sources
        paths_to_all_nodes_in_subtree = {
            path: node
            for path, node in self.subtree_with_keys
            if path != "."  # exclude the root node
        }

        all_item_sources = itertools.chain(
            items_on_this_node, [paths_to_all_nodes_in_subtree]
        )

        items = {
            item
            for source in all_item_sources
            for item in source
            if isinstance(item, str)
        }
        return list(items)

    def __contains__(self, key: object) -> bool:
        """The 'in' operator will return true or false depending on whether
        'key' is either an array stored in the datatree or a child node, or neither.
        """
        return key in self.variables or key in self.children

    def __bool__(self) -> bool:
        return bool(self._data_variables) or bool(self._children)

    def __iter__(self) -> Iterator[str]:
        return itertools.chain(self._data_variables, self._children)  # type: ignore[arg-type]

    def __array__(
        self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        raise TypeError(
            "cannot directly convert a DataTree into a "
            "numpy array. Instead, create an xarray.DataArray "
            "first, either with indexing on the DataTree or by "
            "invoking the `to_array()` method."
        )

    def __repr__(self) -> str:  # type: ignore[override]
        return datatree_repr(self)

    def __str__(self) -> str:
        return datatree_repr(self)

    def _repr_html_(self):
        """Make html representation of datatree object"""
        if XR_OPTS["display_style"] == "text":
            return f"<pre>{escape(repr(self))}</pre>"
        return datatree_repr_html(self)

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    # DatasetView does not support close() or set_close(), so we reimplement
    # these methods on DataTree.

    def _close_node(self) -> None:
        if self._close is not None:
            self._close()
        self._close = None

    def close(self) -> None:
        """Close any files associated with this tree."""
        for node in self.subtree:
            node._close_node()

    def set_close(self, close: Callable[[], None] | None) -> None:
        """Set the closer for this node."""
        self._close = close

    def _replace_node(
        self: DataTree,
        data: Dataset | Default = _default,
        children: dict[str, DataTree] | Default = _default,
    ) -> None:
        ds = self.to_dataset(inherit=False) if data is _default else data

        if children is _default:
            children = self._children

        for child_name in children:
            if child_name in ds.variables:
                raise ValueError(f"node already contains a variable named {child_name}")

        parent_ds = (
            self.parent._to_dataset_view(rebuild_dims=False, inherit=True)
            if self.parent is not None
            else None
        )
        check_alignment(self.path, ds, parent_ds, children)

        if data is not _default:
            self._set_node_data(ds)

        if self.parent is not None:
            _deduplicate_inherited_coordinates(self, self.parent)

        self.children = children

    def _copy_node(
        self, inherit: bool, deep: bool = False, memo: dict[int, Any] | None = None
    ) -> Self:
        """Copy just one node of a tree."""
        new_node = super()._copy_node(inherit=inherit, deep=deep, memo=memo)
        data = self._to_dataset_view(rebuild_dims=False, inherit=inherit)._copy(
            deep=deep, memo=memo
        )
        new_node._set_node_data(data)
        return new_node

    def get(  # type: ignore[override]
        self: DataTree, key: str, default: DataTree | DataArray | None = None
    ) -> DataTree | DataArray | None:
        """
        Access child nodes, variables, or coordinates stored in this node.

        Returned object will be either a DataTree or DataArray object depending on whether the key given points to a
        child or variable.

        Parameters
        ----------
        key : str
            Name of variable / child within this node. Must lie in this immediate node (not elsewhere in the tree).
        default : DataTree | DataArray | None, optional
            A value to return if the specified key does not exist. Default return value is None.
        """
        if key in self.children:
            return self.children[key]
        elif key in self.dataset:
            return self.dataset[key]
        else:
            return default

    def __getitem__(self: DataTree, key: str) -> DataTree | DataArray:
        """
        Access child nodes, variables, or coordinates stored anywhere in this tree.

        Returned object will be either a DataTree or DataArray object depending on whether the key given points to a
        child or variable.

        Parameters
        ----------
        key : str
            Name of variable / child within this node, or unix-like path to variable / child within another node.

        Returns
        -------
        DataTree | DataArray
        """

        # Either:
        if utils.is_dict_like(key):
            # dict-like indexing
            raise NotImplementedError("Should this index over whole tree?")
        elif isinstance(key, str):
            # TODO should possibly deal with hashables in general?
            # path-like: a name of a node/variable, or path to a node/variable
            path = NodePath(key)
            return self._get_item(path)
        elif utils.is_list_like(key):
            # iterable of variable names
            raise NotImplementedError(
                "Selecting via tags is deprecated, and selecting multiple items should be "
                "implemented via .subset"
            )
        else:
            raise ValueError(f"Invalid format for key: {key}")

    def _set(self, key: str, val: DataTree | CoercibleValue) -> None:
        """
        Set the child node or variable with the specified key to value.

        Counterpart to the public .get method, and also only works on the immediate node, not other nodes in the tree.
        """
        if isinstance(val, DataTree):
            # create and assign a shallow copy here so as not to alter original name of node in grafted tree
            new_node = val.copy(deep=False)
            new_node.name = key
            new_node._set_parent(new_parent=self, child_name=key)
        else:
            if not isinstance(val, DataArray | Variable):
                # accommodate other types that can be coerced into Variables
                val = DataArray(val)

            self.update({key: val})

    def __setitem__(
        self,
        key: str,
        value: Any,
    ) -> None:
        """
        Add either a child node or an array to the tree, at any position.

        Data can be added anywhere, and new nodes will be created to cross the path to the new location if necessary.

        If there is already a node at the given location, then if value is a Node class or Dataset it will overwrite the
        data already present at that node, and if value is a single array, it will be merged with it.
        """
        # TODO xarray.Dataset accepts other possibilities, how do we exactly replicate all the behaviour?
        if utils.is_dict_like(key):
            raise NotImplementedError
        elif isinstance(key, str):
            # TODO should possibly deal with hashables in general?
            # path-like: a name of a node/variable, or path to a node/variable
            path = NodePath(key)
            if isinstance(value, Dataset):
                value = DataTree(dataset=value)
            return self._set_item(path, value, new_nodes_along_path=True)
        else:
            raise ValueError("Invalid format for key")

    def __delitem__(self, key: str) -> None:
        """Remove a variable or child node from this datatree node."""
        if key in self.children:
            super().__delitem__(key)

        elif key in self._node_coord_variables:
            if key in self._node_indexes:
                del self._node_indexes[key]
            del self._node_coord_variables[key]
            self._node_dims = calculate_dimensions(self.variables)

        elif key in self._data_variables:
            del self._data_variables[key]
            self._node_dims = calculate_dimensions(self.variables)

        else:
            raise KeyError(key)

    @overload
    def update(self, other: Dataset) -> None: ...

    @overload
    def update(self, other: Mapping[Hashable, DataArray | Variable]) -> None: ...

    @overload
    def update(self, other: Mapping[str, DataTree | DataArray | Variable]) -> None: ...

    def update(
        self,
        other: (
            Dataset
            | Mapping[Hashable, DataArray | Variable]
            | Mapping[str, DataTree | DataArray | Variable]
        ),
    ) -> None:
        """
        Update this node's children and / or variables.

        Just like `dict.update` this is an in-place operation.
        """
        new_children: dict[str, DataTree] = {}
        new_variables: CoercibleMapping

        if isinstance(other, Dataset):
            new_variables = other
        else:
            new_variables = {}
            for k, v in other.items():
                if isinstance(v, DataTree):
                    # avoid named node being stored under inconsistent key
                    new_child: DataTree = v.copy()
                    # Datatree's name is always a string until we fix that (#8836)
                    new_child.name = str(k)
                    new_children[str(k)] = new_child
                elif isinstance(v, DataArray | Variable):
                    # TODO this should also accommodate other types that can be coerced into Variables
                    new_variables[k] = v
                else:
                    raise TypeError(f"Type {type(v)} cannot be assigned to a DataTree")

        vars_merge_result = dataset_update_method(
            self.to_dataset(inherit=False), new_variables
        )
        data = Dataset._construct_direct(**vars_merge_result._asdict())

        # TODO are there any subtleties with preserving order of children like this?
        merged_children = {**self.children, **new_children}

        self._replace_node(data, children=merged_children)

    def assign(
        self, items: Mapping[Any, Any] | None = None, **items_kwargs: Any
    ) -> DataTree:
        """
        Assign new data variables or child nodes to a DataTree, returning a new object
        with all the original items in addition to the new ones.

        Parameters
        ----------
        items : mapping of hashable to Any
            Mapping from variable or child node names to the new values. If the new values
            are callable, they are computed on the Dataset and assigned to new
            data variables. If the values are not callable, (e.g. a DataTree, DataArray,
            scalar, or array), they are simply assigned.
        **items_kwargs
            The keyword arguments form of ``variables``.
            One of variables or variables_kwargs must be provided.

        Returns
        -------
        dt : DataTree
            A new DataTree with the new variables or children in addition to all the
            existing items.

        Notes
        -----
        Since ``kwargs`` is a dictionary, the order of your arguments may not
        be preserved, and so the order of the new variables is not well-defined.
        Assigning multiple items within the same ``assign`` is
        possible, but you cannot reference other variables created within the
        same ``assign`` call.

        See Also
        --------
        xarray.Dataset.assign
        pandas.DataFrame.assign
        """
        items = either_dict_or_kwargs(items, items_kwargs, "assign")
        dt = self.copy()
        dt.update(items)
        return dt

    def drop_nodes(
        self: DataTree, names: str | Iterable[str], *, errors: ErrorOptions = "raise"
    ) -> DataTree:
        """
        Drop child nodes from this node.

        Parameters
        ----------
        names : str or iterable of str
            Name(s) of nodes to drop.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a KeyError if any of the node names
            passed are not present as children of this node. If 'ignore',
            any given names that are present are dropped and no error is raised.

        Returns
        -------
        dropped : DataTree
            A copy of the node with the specified children dropped.
        """
        # the Iterable check is required for mypy
        if isinstance(names, str) or not isinstance(names, Iterable):
            names = {names}
        else:
            names = set(names)

        if errors == "raise":
            extra = names - set(self.children)
            if extra:
                raise KeyError(f"Cannot drop all nodes - nodes {extra} not present")

        result = self.copy()
        children_to_keep = {
            name: child for name, child in result.children.items() if name not in names
        }
        result._replace_node(children=children_to_keep)
        return result

    @classmethod
    def from_dict(
        cls,
        d: Mapping[str, Dataset | DataTree | None],
        /,
        name: str | None = None,
    ) -> Self:
        """
        Create a datatree from a dictionary of data objects, organised by paths into the tree.

        Parameters
        ----------
        d : dict-like
            A mapping from path names to xarray.Dataset or DataTree objects.

            Path names are to be given as unix-like path. If path names
            containing more than one part are given, new tree nodes will be
            constructed as necessary.

            To assign data to the root node of the tree use "", ".", "/" or "./"
            as the path.
        name : Hashable | None, optional
            Name for the root node of the tree. Default is None.

        Returns
        -------
        DataTree

        Notes
        -----
        If your dictionary is nested you will need to flatten it before using this method.
        """
        # Find any values corresponding to the root
        d_cast = dict(d)
        root_data = None
        for key in ("", ".", "/", "./"):
            if key in d_cast:
                if root_data is not None:
                    raise ValueError(
                        "multiple entries found corresponding to the root node"
                    )
                root_data = d_cast.pop(key)

        # Create the root node
        if isinstance(root_data, DataTree):
            obj = root_data.copy()
            obj.name = name
        elif root_data is None or isinstance(root_data, Dataset):
            obj = cls(name=name, dataset=root_data, children=None)
        else:
            raise TypeError(
                f'root node data (at "", ".", "/" or "./") must be a Dataset '
                f"or DataTree, got {type(root_data)}"
            )

        def depth(item) -> int:
            pathstr, _ = item
            return len(NodePath(pathstr).parts)

        if d_cast:
            # Populate tree with children determined from data_objects mapping
            # Sort keys by depth so as to insert nodes from root first (see GH issue #9276)
            for path, data in sorted(d_cast.items(), key=depth):
                # Create and set new node
                if isinstance(data, DataTree):
                    new_node = data.copy()
                elif isinstance(data, Dataset) or data is None:
                    new_node = cls(dataset=data)
                else:
                    raise TypeError(f"invalid values: {data}")
                obj._set_item(
                    path,
                    new_node,
                    allow_overwrite=False,
                    new_nodes_along_path=True,
                )

        # TODO: figure out why mypy is raising an error here, likely something
        # to do with the return type of Dataset.copy()
        return obj  # type: ignore[return-value]

    def to_dict(self, relative: bool = False) -> dict[str, Dataset]:
        """
        Create a dictionary mapping of paths to the data contained in those nodes.

        Parameters
        ----------
        relative : bool
            If True, return relative instead of absolute paths.

        Returns
        -------
        dict[str, Dataset]

        See also
        --------
        DataTree.subtree_with_keys
        """
        return {
            node.relative_to(self) if relative else node.path: node.to_dataset()
            for node in self.subtree
        }

    @property
    def nbytes(self) -> int:
        return sum(node.to_dataset().nbytes for node in self.subtree)

    def __len__(self) -> int:
        return len(self.children) + len(self.data_vars)

    @property
    def indexes(self) -> Indexes[pd.Index]:
        """Mapping of pandas.Index objects used for label based indexing.

        Raises an error if this DataTree node has indexes that cannot be coerced
        to pandas.Index objects.

        See Also
        --------
        DataTree.xindexes
        """
        return self.xindexes.to_pandas_indexes()

    @property
    def xindexes(self) -> Indexes[Index]:
        """Mapping of xarray Index objects used for label based indexing."""
        return Indexes(
            self._indexes, {k: self._coord_variables[k] for k in self._indexes}
        )

    @property
    def coords(self) -> DataTreeCoordinates:
        """Dictionary of xarray.DataArray objects corresponding to coordinate
        variables
        """
        return DataTreeCoordinates(self)

    @property
    def data_vars(self) -> DataVariables:
        """Dictionary of DataArray objects corresponding to data variables"""
        return DataVariables(self.to_dataset())

    def isomorphic(self, other: DataTree) -> bool:
        """
        Two DataTrees are considered isomorphic if the set of paths to their
        descendent nodes are the same.

        Nothing about the data in each node is checked.

        Isomorphism is a necessary condition for two trees to be used in a nodewise binary operation,
        such as ``tree1 + tree2``.

        Parameters
        ----------
        other : DataTree
            The other tree object to compare to.

        See Also
        --------
        DataTree.equals
        DataTree.identical
        """
        return diff_treestructure(self, other) is None

    def equals(self, other: DataTree) -> bool:
        """
        Two DataTrees are equal if they have isomorphic node structures, with
        matching node names, and if they have matching variables and
        coordinates, all of which are equal.

        Parameters
        ----------
        other : DataTree
            The other tree object to compare to.

        See Also
        --------
        Dataset.equals
        DataTree.isomorphic
        DataTree.identical
        """
        if not self.isomorphic(other):
            return False

        # Note: by using .dataset, this intentionally does not check that
        # coordinates are defined at the same levels.
        return all(
            node.dataset.equals(other_node.dataset)
            for node, other_node in zip_subtrees(self, other)
        )

    def _inherited_coords_set(self) -> set[str]:
        return set(self.parent.coords if self.parent else [])  # type: ignore[arg-type]

    def identical(self, other: DataTree) -> bool:
        """
        Like equals, but also checks attributes on all datasets, variables and
        coordinates, and requires that any inherited coordinates at the tree
        root are also inherited on the other tree.

        Parameters
        ----------
        other : DataTree
            The other tree object to compare to.

        See Also
        --------
        Dataset.identical
        DataTree.isomorphic
        DataTree.equals
        """
        if not self.isomorphic(other):
            return False

        if self.name != other.name:
            return False

        if self._inherited_coords_set() != other._inherited_coords_set():
            return False

        return all(
            node.dataset.identical(other_node.dataset)
            for node, other_node in zip_subtrees(self, other)
        )

    def filter(self: DataTree, filterfunc: Callable[[DataTree], bool]) -> DataTree:
        """
        Filter nodes according to a specified condition.

        Returns a new tree containing only the nodes in the original tree for which `fitlerfunc(node)` is True.
        Will also contain empty nodes at intermediate positions if required to support leaves.

        Parameters
        ----------
        filterfunc: function
            A function which accepts only one DataTree - the node on which filterfunc will be called.

        Returns
        -------
        DataTree

        See Also
        --------
        match
        pipe
        map_over_datasets
        """
        filtered_nodes = {
            path: node.dataset
            for path, node in self.subtree_with_keys
            if filterfunc(node)
        }
        return DataTree.from_dict(filtered_nodes, name=self.name)

    def filter_like(self, other: DataTree) -> DataTree:
        """
        Filter a datatree like another datatree.

        Returns a new tree containing only the nodes in the original tree which are also present in the other tree.

        Parameters
        ----------
        other : DataTree
            The tree to filter this tree by.

        Returns
        -------
        DataTree

        See Also
        --------
        filter
        isomorphic

        Examples
        --------

        >>> dt = DataTree.from_dict(
        ...     {
        ...         "/a/A": None,
        ...         "/a/B": None,
        ...         "/b/A": None,
        ...         "/b/B": None,
        ...     }
        ... )
        >>> other = DataTree.from_dict(
        ...     {
        ...         "/a/A": None,
        ...         "/b/A": None,
        ...     }
        ... )
        >>> dt.filter_like(other)
        <xarray.DataTree>
        Group: /
        ├── Group: /a
        │   └── Group: /a/A
        └── Group: /b
            └── Group: /b/A
        """
        other_keys = {key for key, _ in other.subtree_with_keys}
        return self.filter(lambda node: node.relative_to(self) in other_keys)

    def prune(self, drop_size_zero_vars: bool = False) -> DataTree:
        """
        Remove empty nodes from the tree.

        Returns a new tree containing only nodes that contain data variables with actual data.
        Intermediate nodes are kept if they are required to support non-empty children.

        Parameters
        ----------
        drop_size_zero_vars : bool, default False
            If True, also considers variables with zero size as empty.
            If False, keeps nodes with data variables even if they have zero size.

        Returns
        -------
        DataTree
            A new tree with empty nodes removed.

        See Also
        --------
        filter

        Examples
        --------
        >>> dt = xr.DataTree.from_dict(
        ...     {
        ...         "/a": xr.Dataset({"foo": ("x", [1, 2])}),
        ...         "/b": xr.Dataset({"bar": ("x", [])}),
        ...         "/c": xr.Dataset(),
        ...     }
        ... )
        >>> dt.prune()  # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
        <xarray.DataTree>
        Group: /
        ├── Group: /a
        │       Dimensions:  (x: 2)
        │       Dimensions without coordinates: x
        │       Data variables:
        │           foo      (x) int64 16B 1 2
        └── Group: /b
                Dimensions:  (x: 0)
                Dimensions without coordinates: x
                Data variables:
                    bar      (x) float64 0B...

        The ``drop_size_zero_vars`` parameter controls whether variables
        with zero size are considered empty:

        >>> dt.prune(drop_size_zero_vars=True)
        <xarray.DataTree>
        Group: /
        └── Group: /a
                Dimensions:  (x: 2)
                Dimensions without coordinates: x
                Data variables:
                    foo      (x) int64 16B 1 2
        """
        non_empty_cond: Callable[[DataTree], bool]
        if drop_size_zero_vars:
            non_empty_cond = lambda node: len(node.data_vars) > 0 and any(
                var.size > 0 for var in node.data_vars.values()
            )
        else:
            non_empty_cond = lambda node: len(node.data_vars) > 0

        return self.filter(non_empty_cond)

    def match(self, pattern: str) -> DataTree:
        """
        Return nodes with paths matching pattern.

        Uses unix glob-like syntax for pattern-matching.

        Parameters
        ----------
        pattern: str
            A pattern to match each node path against.

        Returns
        -------
        DataTree

        See Also
        --------
        filter
        pipe
        map_over_datasets

        Examples
        --------
        >>> dt = DataTree.from_dict(
        ...     {
        ...         "/a/A": None,
        ...         "/a/B": None,
        ...         "/b/A": None,
        ...         "/b/B": None,
        ...     }
        ... )
        >>> dt.match("*/B")
        <xarray.DataTree>
        Group: /
        ├── Group: /a
        │   └── Group: /a/B
        └── Group: /b
            └── Group: /b/B
        """
        matching_nodes = {
            path: node.dataset
            for path, node in self.subtree_with_keys
            if NodePath(node.path).match(pattern)
        }
        return DataTree.from_dict(matching_nodes, name=self.name)

    @overload
    def map_over_datasets(
        self,
        func: Callable[..., Dataset | None],
        *args: Any,
        kwargs: Mapping[str, Any] | None = None,
    ) -> DataTree: ...

    @overload
    def map_over_datasets(
        self,
        func: Callable[..., tuple[Dataset | None, Dataset | None]],
        *args: Any,
        kwargs: Mapping[str, Any] | None = None,
    ) -> tuple[DataTree, DataTree]: ...

    @overload
    def map_over_datasets(
        self,
        func: Callable[..., tuple[Dataset | None, ...]],
        *args: Any,
        kwargs: Mapping[str, Any] | None = None,
    ) -> tuple[DataTree, ...]: ...

    def map_over_datasets(
        self,
        func: Callable[..., Dataset | None | tuple[Dataset | None, ...]],
        *args: Any,
        kwargs: Mapping[str, Any] | None = None,
    ) -> DataTree | tuple[DataTree, ...]:
        """
        Apply a function to every dataset in this subtree, returning a new tree which stores the results.

        The function will be applied to any dataset stored in this node, as well as any dataset stored in any of the
        descendant nodes. The returned tree will have the same structure as the original subtree.

        func needs to return a Dataset in order to rebuild the subtree.

        Parameters
        ----------
        func : callable
            Function to apply to datasets with signature:
            `func(node.dataset, *args, **kwargs) -> Dataset`.

            Function will not be applied to any nodes without datasets.
        *args : tuple, optional
            Positional arguments passed on to `func`. Any DataTree arguments will be
            converted to Dataset objects via `.dataset`.
        kwargs : dict, optional
            Optional keyword arguments passed directly to ``func``.

        Returns
        -------
        subtrees : DataTree, tuple of DataTrees
            One or more subtrees containing results from applying ``func`` to the data at each node.

        See also
        --------
        map_over_datasets
        """
        # TODO this signature means that func has no way to know which node it is being called upon - change?
        return map_over_datasets(func, self, *args, kwargs=kwargs)  # type: ignore[arg-type]

    @overload
    def pipe(
        self,
        func: Callable[Concatenate[Self, P], T],
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> T: ...

    @overload
    def pipe(
        self,
        func: tuple[Callable[..., T], str],
        *args: Any,
        **kwargs: Any,
    ) -> T: ...

    def pipe(
        self,
        func: Callable[Concatenate[Self, P], T] | tuple[Callable[..., T], str],
        *args: Any,
        **kwargs: Any,
    ) -> T:
        """Apply ``func(self, *args, **kwargs)``

        This method replicates the pandas method of the same name.

        Parameters
        ----------
        func : callable
            function to apply to this xarray object (Dataset/DataArray).
            ``args``, and ``kwargs`` are passed into ``func``.
            Alternatively a ``(callable, data_keyword)`` tuple where
            ``data_keyword`` is a string indicating the keyword of
            ``callable`` that expects the xarray object.
        *args
            positional arguments passed into ``func``.
        **kwargs
            a dictionary of keyword arguments passed into ``func``.

        Returns
        -------
        object : T
            the return type of ``func``.

        Notes
        -----
        Use ``.pipe`` when chaining together functions that expect
        xarray or pandas objects, e.g., instead of writing

        .. code:: python

            f(g(h(dt), arg1=a), arg2=b, arg3=c)

        You can write

        .. code:: python

            (dt.pipe(h).pipe(g, arg1=a).pipe(f, arg2=b, arg3=c))

        If you have a function that takes the data as (say) the second
        argument, pass a tuple indicating which keyword expects the
        data. For example, suppose ``f`` takes its data as ``arg2``:

        .. code:: python

            (dt.pipe(h).pipe(g, arg1=a).pipe((f, "arg2"), arg1=a, arg3=c))

        """
        if isinstance(func, tuple):
            # Use different var when unpacking function from tuple because the type
            # signature of the unpacked function differs from the expected type
            # signature in the case where only a function is given, rather than a tuple.
            # This makes type checkers happy at both call sites below.
            f, target = func
            if target in kwargs:
                raise ValueError(
                    f"{target} is both the pipe target and a keyword argument"
                )
            kwargs[target] = self
            return f(*args, **kwargs)

        return func(self, *args, **kwargs)

    # TODO some kind of .collapse() or .flatten() method to merge a subtree

    @property
    def groups(self):
        """Return all groups in the tree, given as a tuple of path-like strings."""
        return tuple(node.path for node in self.subtree)

    def _unary_op(self, f, *args, **kwargs) -> DataTree:
        # TODO do we need to any additional work to avoid duplication etc.? (Similar to aggregations)
        return self.map_over_datasets(functools.partial(f, **kwargs), *args)

    def _binary_op(self, other, f, reflexive=False, join=None) -> DataTree:
        from xarray.core.groupby import GroupBy

        if isinstance(other, GroupBy):
            return NotImplemented

        ds_binop = functools.partial(
            Dataset._binary_op,
            f=f,
            reflexive=reflexive,
            join=join,
        )
        return map_over_datasets(ds_binop, self, other)

    def _inplace_binary_op(self, other, f) -> Self:
        from xarray.core.groupby import GroupBy

        if isinstance(other, GroupBy):
            raise TypeError(
                "in-place operations between a DataTree and "
                "a grouped object are not permitted"
            )

        # TODO see GH issue #9629 for required implementation
        raise NotImplementedError()

    # TODO: dirty workaround for mypy 1.5 error with inherited DatasetOpsMixin vs. Mapping
    # related to https://github.com/python/mypy/issues/9319?
    def __eq__(self, other: DtCompatible) -> Self:  # type: ignore[override]
        return super().__eq__(other)

    # filepath=None writes to a memoryview
    @overload
    def to_netcdf(
        self,
        filepath: None = None,
        mode: NetcdfWriteModes = "w",
        encoding=None,
        unlimited_dims=None,
        format: T_DataTreeNetcdfTypes | None = None,
        engine: T_DataTreeNetcdfEngine | None = None,
        group: str | None = None,
        write_inherited_coords: bool = False,
        compute: bool = True,
        **kwargs,
    ) -> memoryview: ...

    # compute=False returns dask.Delayed
    @overload
    def to_netcdf(
        self,
        filepath: str | PathLike | io.IOBase,
        mode: NetcdfWriteModes = "w",
        encoding=None,
        unlimited_dims=None,
        format: T_DataTreeNetcdfTypes | None = None,
        engine: T_DataTreeNetcdfEngine | None = None,
        group: str | None = None,
        write_inherited_coords: bool = False,
        *,
        compute: Literal[False],
        **kwargs,
    ) -> Delayed: ...

    # default return None
    @overload
    def to_netcdf(
        self,
        filepath: str | PathLike | io.IOBase,
        mode: NetcdfWriteModes = "w",
        encoding=None,
        unlimited_dims=None,
        format: T_DataTreeNetcdfTypes | None = None,
        engine: T_DataTreeNetcdfEngine | None = None,
        group: str | None = None,
        write_inherited_coords: bool = False,
        compute: Literal[True] = True,
        **kwargs,
    ) -> None: ...

    def to_netcdf(
        self,
        filepath: str | PathLike | io.IOBase | None = None,
        mode: NetcdfWriteModes = "w",
        encoding=None,
        unlimited_dims=None,
        format: T_DataTreeNetcdfTypes | None = None,
        engine: T_DataTreeNetcdfEngine | None = None,
        group: str | None = None,
        write_inherited_coords: bool = False,
        compute: bool = True,
        **kwargs,
    ) -> None | memoryview | Delayed:
        """
        Write datatree contents to a netCDF file.

        Parameters
        ----------
        filepath : str or PathLike or file-like object or None
            Path to which to save this datatree, or a file-like object to write
            it to (which must support read and write and be seekable) or None
            to return in-memory bytes as a memoryview.
        mode : {"w", "a"}, default: "w"
            Write ('w') or append ('a') mode. If mode='w', any existing file at
            this location will be overwritten. If mode='a', existing variables
            will be overwritten. Only applies to the root group.
        encoding : dict, optional
            Nested dictionary with variable names as keys and dictionaries of
            variable specific encodings as values, e.g.,
            ``{"root/set1": {"my_variable": {"dtype": "int16", "scale_factor": 0.1,
            "zlib": True}, ...}, ...}``. See ``xarray.Dataset.to_netcdf`` for available
            options.
        unlimited_dims : dict, optional
            Mapping of unlimited dimensions per group that that should be serialized as unlimited dimensions.
            By default, no dimensions are treated as unlimited dimensions.
            Note that unlimited_dims may also be set via
            ``dataset.encoding["unlimited_dims"]``.
        format : {"NETCDF4", }, optional
            File format for the resulting netCDF file:

            * NETCDF4: Data is stored in an HDF5 file, using netCDF4 API features.
        engine : {"netcdf4", "h5netcdf"}, optional
            Engine to use when writing netCDF files. If not provided, the
            default engine is chosen based on available dependencies, with a
            preference for "netcdf4" if writing to a file on disk.
        group : str, optional
            Path to the netCDF4 group in the given file to open as the root group
            of the ``DataTree``. Currently, specifying a group is not supported.
        write_inherited_coords : bool, default: False
            If true, replicate inherited coordinates on all descendant nodes.
            Otherwise, only write coordinates at the level at which they are
            originally defined. This saves disk space, but requires opening the
            full tree to load inherited coordinates.
        compute : bool, default: True
            If true compute immediately, otherwise return a
            ``dask.delayed.Delayed`` object that can be computed later.
        kwargs :
            Additional keyword arguments to be passed to ``xarray.Dataset.to_netcdf``

        Returns
        -------
            * ``memoryview`` if path is None
            * ``dask.delayed.Delayed`` if compute is False
            * ``None`` otherwise

        Note
        ----
            Due to file format specifications the on-disk root group name
            is always ``"/"`` overriding any given ``DataTree`` root node name.
        """
        from xarray.core.datatree_io import _datatree_to_netcdf

        return _datatree_to_netcdf(
            self,
            filepath,
            mode=mode,
            encoding=encoding,
            unlimited_dims=unlimited_dims,
            format=format,
            engine=engine,
            group=group,
            write_inherited_coords=write_inherited_coords,
            compute=compute,
            **kwargs,
        )

    # compute=False returns dask.Delayed
    @overload
    def to_zarr(
        self,
        store,
        mode: ZarrWriteModes = "w-",
        encoding=None,
        consolidated: bool = True,
        group: str | None = None,
        write_inherited_coords: bool = False,
        *,
        compute: Literal[False],
        **kwargs,
    ) -> Delayed: ...

    # default returns ZarrStore
    @overload
    def to_zarr(
        self,
        store,
        mode: ZarrWriteModes = "w-",
        encoding=None,
        consolidated: bool = True,
        group: str | None = None,
        write_inherited_coords: bool = False,
        compute: Literal[True] = True,
        **kwargs,
    ) -> ZarrStore: ...

    def to_zarr(
        self,
        store,
        mode: ZarrWriteModes = "w-",
        encoding=None,
        consolidated: bool = True,
        group: str | None = None,
        write_inherited_coords: bool = False,
        compute: bool = True,
        **kwargs,
    ) -> ZarrStore | Delayed:
        """
        Write datatree contents to a Zarr store.

        Parameters
        ----------
        store : MutableMapping, str or Path, optional
            Store or path to directory in file system
        mode : {{"w", "w-", "a", "r+", None}, default: "w-"
            Persistence mode: “w” means create (overwrite if exists); “w-” means create (fail if exists);
            “a” means override existing variables (create if does not exist); “r+” means modify existing
            array values only (raise an error if any metadata or shapes would change). The default mode
            is “w-”.
        encoding : dict, optional
            Nested dictionary with variable names as keys and dictionaries of
            variable specific encodings as values, e.g.,
            ``{"root/set1": {"my_variable": {"dtype": "int16", "scale_factor": 0.1}, ...}, ...}``.
            See ``xarray.Dataset.to_zarr`` for available options.
        consolidated : bool
            If True, apply zarr's `consolidate_metadata` function to the store
            after writing metadata for all groups.
        group : str, optional
            Group path. (a.k.a. `path` in zarr terminology.)
        write_inherited_coords : bool, default: False
            If true, replicate inherited coordinates on all descendant nodes.
            Otherwise, only write coordinates at the level at which they are
            originally defined. This saves disk space, but requires opening the
            full tree to load inherited coordinates.
        compute : bool, default: True
            If true compute immediately, otherwise return a
            ``dask.delayed.Delayed`` object that can be computed later. Metadata
            is always updated eagerly.
        kwargs :
            Additional keyword arguments to be passed to ``xarray.Dataset.to_zarr``

        Note
        ----
            Due to file format specifications the on-disk root group name
            is always ``"/"`` overriding any given ``DataTree`` root node name.
        """
        from xarray.core.datatree_io import _datatree_to_zarr

        return _datatree_to_zarr(
            self,
            store,
            mode=mode,
            encoding=encoding,
            consolidated=consolidated,
            group=group,
            write_inherited_coords=write_inherited_coords,
            compute=compute,
            **kwargs,
        )

    def _get_all_dims(self) -> set:
        all_dims: set[Any] = set()
        for node in self.subtree:
            all_dims.update(node._node_dims)
        return all_dims

    def reduce(
        self,
        func: Callable,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        numeric_only: bool = False,
        **kwargs: Any,
    ) -> Self:
        """Reduce this tree by applying `func` along some dimension(s)."""
        dims = parse_dims_as_set(dim, self._get_all_dims())
        result = {}
        for path, node in self.subtree_with_keys:
            reduce_dims = [d for d in node._node_dims if d in dims]
            node_result = node.dataset.reduce(
                func,
                reduce_dims,
                keep_attrs=keep_attrs,
                keepdims=keepdims,
                numeric_only=numeric_only,
                **kwargs,
            )
            result[path] = node_result
        return type(self).from_dict(result, name=self.name)

    def _selective_indexing(
        self,
        func: Callable[[Dataset, Mapping[Any, Any]], Dataset],
        indexers: Mapping[Any, Any],
        missing_dims: ErrorOptionsWithWarn = "raise",
    ) -> Self:
        """Apply an indexing operation over the subtree, handling missing
        dimensions and inherited coordinates gracefully by only applying
        indexing at each node selectively.
        """
        all_dims = self._get_all_dims()
        indexers = drop_dims_from_indexers(indexers, all_dims, missing_dims)

        result = {}
        for path, node in self.subtree_with_keys:
            node_indexers = {k: v for k, v in indexers.items() if k in node.dims}
            func_with_error_context = _handle_errors_with_path_context(path)(func)
            node_result = func_with_error_context(node.dataset, node_indexers)
            # Indexing datasets corresponding to each node results in redundant
            # coordinates when indexes from a parent node are inherited.
            # Ideally, we would avoid creating such coordinates in the first
            # place, but that would require implementing indexing operations at
            # the Variable instead of the Dataset level.
            if node is not self:
                for k in node_indexers:
                    if k not in node._node_coord_variables and k in node_result.coords:
                        # We remove all inherited coordinates. Coordinates
                        # corresponding to an index would be de-duplicated by
                        # _deduplicate_inherited_coordinates(), but indexing (e.g.,
                        # with a scalar) can also create scalar coordinates, which
                        # need to be explicitly removed.
                        del node_result.coords[k]
            result[path] = node_result
        return type(self).from_dict(result, name=self.name)

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        drop: bool = False,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new data tree with each array indexed along the specified
        dimension(s).

        This method selects values from each array using its `__getitem__`
        method, except this method does not require knowing the order of
        each array's dimensions.

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by integers, slice objects or arrays.
            indexer can be a integer, slice, array-like or DataArray.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        drop : bool, default: False
            If ``drop=True``, drop coordinates variables indexed by integers
            instead of making them scalar.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Dataset:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : DataTree
            A new DataTree with the same contents as this data tree, except each
            array and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        DataTree.sel
        Dataset.isel
        """

        def apply_indexers(dataset, node_indexers):
            return dataset.isel(node_indexers, drop=drop)

        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")
        return self._selective_indexing(
            apply_indexers, indexers, missing_dims=missing_dims
        )

    def sel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance: int | float | Iterable[int | float] | None = None,
        drop: bool = False,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new data tree with each array indexed by tick labels
        along the specified dimension(s).

        In contrast to `DataTree.isel`, indexers for this method should use
        labels instead of integers.

        Under the hood, this method is powered by using pandas's powerful Index
        objects. This makes label based indexing essentially just as fast as
        using integer indexing.

        It also means this method uses pandas's (well documented) logic for
        indexing. This means you can use string shortcuts for datetime indexes
        (e.g., '2000-01' to select all values in January 2000). It also means
        that slices are treated as inclusive of both the start and stop values,
        unlike normal Python indexing.

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by scalars, slices or arrays of tick labels. For dimensions with
            multi-index, the indexer may also be a dict-like object with keys
            matching index level names.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method to use for inexact matches:

            * None (default): only exact matches
            * pad / ffill: propagate last valid index value forward
            * backfill / bfill: propagate next valid index value backward
            * nearest: use nearest valid index value
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
        drop : bool, optional
            If ``drop=True``, drop coordinates variables in `indexers` instead
            of making them scalar.
        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : DataTree
            A new DataTree with the same contents as this data tree, except each
            variable and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        DataTree.isel
        Dataset.sel
        """

        def apply_indexers(dataset, node_indexers):
            # TODO: reimplement in terms of map_index_queries(), to avoid
            # redundant look-ups of integer positions from labels (via indexes)
            # on child nodes.
            return dataset.sel(
                node_indexers, method=method, tolerance=tolerance, drop=drop
            )

        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "sel")
        return self._selective_indexing(apply_indexers, indexers)

    def load(self, **kwargs) -> Self:
        """Manually trigger loading and/or computation of this datatree's data
        from disk or a remote source into memory and return this datatree.
        Unlike compute, the original datatree is modified and returned.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically. However, this method can be necessary when
        working with many file objects on disk.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        See Also
        --------
        Dataset.load
        dask.compute
        """
        # access .data to coerce everything to numpy or dask arrays
        lazy_data = {
            path: {
                k: v._data
                for k, v in node.variables.items()
                if is_chunked_array(v._data)
            }
            for path, node in self.subtree_with_keys
        }
        flat_lazy_data = {
            (path, var_name): array
            for path, node in lazy_data.items()
            for var_name, array in node.items()
        }
        if flat_lazy_data:
            chunkmanager = get_chunked_array_type(*flat_lazy_data.values())

            # evaluate all the chunked arrays simultaneously
            evaluated_data: tuple[np.ndarray[Any, Any], ...] = chunkmanager.compute(
                *flat_lazy_data.values(), **kwargs
            )

            for (path, var_name), data in zip(
                flat_lazy_data, evaluated_data, strict=False
            ):
                self[path].variables[var_name].data = data

        # load everything else sequentially
        for node in self.subtree:
            for k, v in node.variables.items():
                if k not in lazy_data:
                    v.load()

        return self

    def compute(self, **kwargs) -> Self:
        """Manually trigger loading and/or computation of this datatree's data
        from disk or a remote source into memory and return a new datatree.
        Unlike load, the original datatree is left unaltered.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically. However, this method can be necessary when
        working with many file objects on disk.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        Returns
        -------
        object : DataTree
            New object with lazy data variables and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        """
        new = self.copy(deep=False)
        return new.load(**kwargs)

    def _persist_inplace(self, **kwargs) -> Self:
        """Persist all chunked arrays in memory"""
        # access .data to coerce everything to numpy or dask arrays
        lazy_data = {
            path: {
                k: v._data
                for k, v in node.variables.items()
                if is_chunked_array(v._data)
            }
            for path, node in self.subtree_with_keys
        }
        flat_lazy_data = {
            (path, var_name): array
            for path, node in lazy_data.items()
            for var_name, array in node.items()
        }
        if flat_lazy_data:
            chunkmanager = get_chunked_array_type(*flat_lazy_data.values())

            # evaluate all the dask arrays simultaneously
            evaluated_data = chunkmanager.persist(*flat_lazy_data.values(), **kwargs)

            for (path, var_name), data in zip(
                flat_lazy_data, evaluated_data, strict=False
            ):
                self[path].variables[var_name].data = data

        return self

    def persist(self, **kwargs) -> Self:
        """Trigger computation, keeping data as chunked arrays.

        This operation can be used to trigger computation on underlying dask
        arrays, similar to ``.compute()`` or ``.load()``.  However this
        operation keeps the data as dask arrays. This is particularly useful
        when using the dask.distributed scheduler and you want to load a large
        amount of data into distributed memory.
        Like compute (but unlike load), the original dataset is left unaltered.


        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.persist``.

        Returns
        -------
        object : DataTree
            New object with all dask-backed coordinates and data variables as persisted dask arrays.

        See Also
        --------
        dask.persist
        """
        new = self.copy(deep=False)
        return new._persist_inplace(**kwargs)

    @property
    def chunksizes(self) -> Mapping[str, Mapping[Hashable, tuple[int, ...]]]:
        """
        Mapping from group paths to a mapping of chunksizes.

        If there's no chunked data in a group, the corresponding mapping of chunksizes will be empty.

        Cannot be modified directly, but can be modified by calling .chunk().

        See Also
        --------
        DataTree.chunk
        Dataset.chunksizes
        """
        return Frozen(
            {
                node.path: get_chunksizes(node.variables.values())
                for node in self.subtree
            }
        )

    def chunk(
        self,
        chunks: T_ChunksFreq = {},  # noqa: B006  # {} even though it's technically unsafe, is being used intentionally here (#4667)
        name_prefix: str = "xarray-",
        token: str | None = None,
        lock: bool = False,
        inline_array: bool = False,
        chunked_array_type: str | ChunkManagerEntrypoint | None = None,
        from_array_kwargs=None,
        **chunks_kwargs: T_ChunkDimFreq,
    ) -> Self:
        """Coerce all arrays in all groups in this tree into dask arrays with the given
        chunks.

        Non-dask arrays in this tree will be converted to dask arrays. Dask
        arrays will be rechunked to the given chunk sizes.

        If neither chunks is not provided for one or more dimensions, chunk
        sizes along that dimension will not be updated; non-dask arrays will be
        converted into dask arrays with a single block.

        Along datetime-like dimensions, a :py:class:`groupers.TimeResampler` object is also accepted.

        Parameters
        ----------
        chunks : int, tuple of int, "auto" or mapping of hashable to int or a TimeResampler, optional
            Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, or
            ``{"x": 5, "y": 5}`` or ``{"x": 5, "time": TimeResampler(freq="YE")}``.
        name_prefix : str, default: "xarray-"
            Prefix for the name of any new dask arrays.
        token : str, optional
            Token uniquely identifying this datatree.
        lock : bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        inline_array: bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        chunked_array_type: str, optional
            Which chunked array type to coerce this datatree's arrays to.
            Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntryPoint` system.
            Experimental API that should not be relied upon.
        from_array_kwargs: dict, optional
            Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
            chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
            For example, with dask as the default chunked array type, this method would pass additional kwargs
            to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
        **chunks_kwargs : {dim: chunks, ...}, optional
            The keyword arguments form of ``chunks``.
            One of chunks or chunks_kwargs must be provided

        Returns
        -------
        chunked : xarray.DataTree

        See Also
        --------
        Dataset.chunk
        Dataset.chunksizes
        xarray.unify_chunks
        dask.array.from_array
        """
        # don't support deprecated ways of passing chunks
        if not isinstance(chunks, Mapping):
            raise TypeError(
                f"invalid type for chunks: {type(chunks)}. Only mappings are supported."
            )
        combined_chunks = either_dict_or_kwargs(chunks, chunks_kwargs, "chunk")

        all_dims = self._get_all_dims()

        bad_dims = combined_chunks.keys() - all_dims
        if bad_dims:
            raise ValueError(
                f"chunks keys {tuple(bad_dims)} not found in data dimensions {tuple(all_dims)}"
            )

        rechunked_groups = {
            path: node.dataset.chunk(
                {
                    dim: size
                    for dim, size in combined_chunks.items()
                    if dim in node._node_dims
                },
                name_prefix=name_prefix,
                token=token,
                lock=lock,
                inline_array=inline_array,
                chunked_array_type=chunked_array_type,
                from_array_kwargs=from_array_kwargs,
            )
            for path, node in self.subtree_with_keys
        }

        return self.from_dict(rechunked_groups, name=self.name)
