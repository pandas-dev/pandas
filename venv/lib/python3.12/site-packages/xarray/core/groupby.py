from __future__ import annotations

import copy
import functools
import itertools
import warnings
from collections.abc import Callable, Hashable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Generic, Literal, Union, cast

import numpy as np
import pandas as pd
from packaging.version import Version

from xarray.computation import ops
from xarray.computation.arithmetic import (
    DataArrayGroupbyArithmetic,
    DatasetGroupbyArithmetic,
)
from xarray.core import dtypes, duck_array_ops, nputils
from xarray.core._aggregations import (
    DataArrayGroupByAggregations,
    DatasetGroupByAggregations,
)
from xarray.core.common import ImplementsArrayReduce, ImplementsDatasetReduce
from xarray.core.coordinates import Coordinates, coordinates_from_variable
from xarray.core.duck_array_ops import where
from xarray.core.formatting import format_array_flat
from xarray.core.indexes import (
    PandasMultiIndex,
    filter_indexes_from_coords,
)
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import (
    Dims,
    QuantileMethods,
    T_DataArray,
    T_DataWithCoords,
    T_Xarray,
)
from xarray.core.utils import (
    FrozenMappingWarningOnValuesAccess,
    contains_only_chunked_or_numpy,
    either_dict_or_kwargs,
    emit_user_level_warning,
    hashable,
    is_scalar,
    maybe_wrap_array,
    module_available,
    peek_at,
)
from xarray.core.variable import IndexVariable, Variable
from xarray.namedarray.pycompat import is_chunked_array
from xarray.structure.alignment import align, broadcast
from xarray.structure.concat import concat
from xarray.structure.merge import merge_coords

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import (
        GroupIndex,
        GroupIndices,
        GroupInput,
        GroupKey,
        T_Chunks,
    )
    from xarray.core.utils import Frozen
    from xarray.groupers import EncodedGroups, Grouper


def check_reduce_dims(reduce_dims, dimensions):
    if reduce_dims is not ...:
        if is_scalar(reduce_dims):
            reduce_dims = [reduce_dims]
        if any(dim not in dimensions for dim in reduce_dims):
            raise ValueError(
                f"cannot reduce over dimensions {reduce_dims!r}. expected either '...' "
                f"to reduce over all dimensions or one or more of {dimensions!r}. "
                f"Alternatively, install the `flox` package. "
            )


def _codes_to_group_indices(codes: np.ndarray, N: int) -> GroupIndices:
    """Converts integer codes for groups to group indices."""
    assert codes.ndim == 1
    groups: GroupIndices = tuple([] for _ in range(N))
    for n, g in enumerate(codes):
        if g >= 0:
            groups[g].append(n)
    return groups


def _dummy_copy(xarray_obj):
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    if isinstance(xarray_obj, Dataset):
        res = Dataset(
            {
                k: dtypes.get_fill_value(v.dtype)
                for k, v in xarray_obj.data_vars.items()
            },
            {
                k: dtypes.get_fill_value(v.dtype)
                for k, v in xarray_obj.coords.items()
                if k not in xarray_obj.dims
            },
            xarray_obj.attrs,
        )
    elif isinstance(xarray_obj, DataArray):
        res = DataArray(
            dtypes.get_fill_value(xarray_obj.dtype),
            {
                k: dtypes.get_fill_value(v.dtype)
                for k, v in xarray_obj.coords.items()
                if k not in xarray_obj.dims
            },
            dims=[],
            name=xarray_obj.name,
            attrs=xarray_obj.attrs,
        )
    else:  # pragma: no cover
        raise AssertionError
    return res


def _is_one_or_none(obj) -> bool:
    return obj == 1 or obj is None


def _consolidate_slices(slices: list[slice]) -> list[slice]:
    """Consolidate adjacent slices in a list of slices."""
    result: list[slice] = []
    last_slice = slice(None)
    for slice_ in slices:
        if not isinstance(slice_, slice):
            raise ValueError(f"list element is not a slice: {slice_!r}")
        if (
            result
            and last_slice.stop == slice_.start
            and _is_one_or_none(last_slice.step)
            and _is_one_or_none(slice_.step)
        ):
            last_slice = slice(last_slice.start, slice_.stop, slice_.step)
            result[-1] = last_slice
        else:
            result.append(slice_)
            last_slice = slice_
    return result


def _inverse_permutation_indices(positions, N: int | None = None) -> np.ndarray | None:
    """Like inverse_permutation, but also handles slices.

    Parameters
    ----------
    positions : list of ndarray or slice
        If slice objects, all are assumed to be slices.

    Returns
    -------
    np.ndarray of indices or None, if no permutation is necessary.
    """
    if not positions:
        return None

    if isinstance(positions[0], slice):
        positions = _consolidate_slices(positions)
        if positions == slice(None):
            return None
        positions = [np.arange(sl.start, sl.stop, sl.step) for sl in positions]

    newpositions = nputils.inverse_permutation(np.concatenate(positions), N)
    return newpositions[newpositions != -1]


class _DummyGroup(Generic[T_Xarray]):
    """Class for keeping track of grouped dimensions without coordinates.

    Should not be user visible.
    """

    __slots__ = ("coords", "dataarray", "name", "size")

    def __init__(self, obj: T_Xarray, name: Hashable, coords) -> None:
        self.name = name
        self.coords = coords
        self.size = obj.sizes[name]

    @property
    def dims(self) -> tuple[Hashable]:
        return (self.name,)

    @property
    def ndim(self) -> Literal[1]:
        return 1

    @property
    def values(self) -> range:
        return range(self.size)

    @property
    def data(self) -> np.ndarray:
        return np.arange(self.size, dtype=int)

    def __array__(
        self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        if copy is False:
            raise NotImplementedError(f"An array copy is necessary, got {copy = }.")
        return np.arange(self.size)

    @property
    def shape(self) -> tuple[int, ...]:
        return (self.size,)

    @property
    def attrs(self) -> dict:
        return {}

    def __getitem__(self, key):
        if isinstance(key, tuple):
            (key,) = key
        return self.values[key]

    def to_index(self) -> pd.Index:
        # could be pd.RangeIndex?
        return pd.Index(np.arange(self.size))

    def copy(self, deep: bool = True, data: Any = None):
        raise NotImplementedError

    def to_dataarray(self) -> DataArray:
        from xarray.core.dataarray import DataArray

        return DataArray(
            data=self.data, dims=(self.name,), coords=self.coords, name=self.name
        )

    def to_array(self) -> DataArray:
        """Deprecated version of to_dataarray."""
        return self.to_dataarray()


T_Group = Union["T_DataArray", _DummyGroup]


def _ensure_1d(
    group: T_Group, obj: T_DataWithCoords
) -> tuple[
    T_Group,
    T_DataWithCoords,
    Hashable | None,
    list[Hashable],
]:
    # 1D cases: do nothing
    if isinstance(group, _DummyGroup) or group.ndim == 1:
        return group, obj, None, []

    from xarray.core.dataarray import DataArray

    if isinstance(group, DataArray):
        for dim in set(group.dims) - set(obj.dims):
            obj = obj.expand_dims(dim)
        # try to stack the dims of the group into a single dim
        orig_dims = group.dims
        stacked_dim = "stacked_" + "_".join(map(str, orig_dims))
        # these dimensions get created by the stack operation
        inserted_dims = [dim for dim in group.dims if dim not in group.coords]
        # `newgroup` construction is optimized so we don't create an index unnecessarily,
        # or stack any non-dim coords unnecessarily
        newgroup = DataArray(group.variable.stack({stacked_dim: orig_dims}))
        newobj = obj.stack({stacked_dim: orig_dims})
        return newgroup, newobj, stacked_dim, inserted_dims

    raise TypeError(f"group must be DataArray or _DummyGroup, got {type(group)!r}.")


@dataclass
class ResolvedGrouper(Generic[T_DataWithCoords]):
    """
    Wrapper around a Grouper object.

    The Grouper object represents an abstract instruction to group an object.
    The ResolvedGrouper object is a concrete version that contains all the common
    logic necessary for a GroupBy problem including the intermediates necessary for
    executing a GroupBy calculation. Specialization to the grouping problem at hand,
    is accomplished by calling the `factorize` method on the encapsulated Grouper
    object.

    This class is private API, while Groupers are public.
    """

    grouper: Grouper
    group: T_Group
    obj: T_DataWithCoords
    eagerly_compute_group: Literal[False] | None = field(repr=False, default=None)

    # returned by factorize:
    encoded: EncodedGroups = field(init=False, repr=False)

    @property
    def full_index(self) -> pd.Index:
        return self.encoded.full_index

    @property
    def codes(self) -> DataArray:
        return self.encoded.codes

    @property
    def unique_coord(self) -> Variable | _DummyGroup:
        return self.encoded.unique_coord

    def __post_init__(self) -> None:
        # This copy allows the BinGrouper.factorize() method
        # to update BinGrouper.bins when provided as int, using the output
        # of pd.cut
        # We do not want to modify the original object, since the same grouper
        # might be used multiple times.
        from xarray.groupers import BinGrouper, UniqueGrouper

        self.grouper = copy.deepcopy(self.grouper)

        self.group = _resolve_group(self.obj, self.group)

        if self.eagerly_compute_group:
            raise ValueError(
                f""""Eagerly computing the DataArray you're grouping by ({self.group.name!r}) "
                has been removed.
                Please load this array's data manually using `.compute` or `.load`.
                To intentionally avoid eager loading, either (1) specify
                `.groupby({self.group.name}=UniqueGrouper(labels=...))`
                or (2) pass explicit bin edges using ``bins`` or
                `.groupby({self.group.name}=BinGrouper(bins=...))`; as appropriate."""
            )
        if self.eagerly_compute_group is not None:
            emit_user_level_warning(
                "Passing `eagerly_compute_group` is now deprecated. It has no effect.",
                DeprecationWarning,
            )

        if not isinstance(self.group, _DummyGroup) and is_chunked_array(
            self.group.variable._data
        ):
            # This requires a pass to discover the groups present
            if isinstance(self.grouper, UniqueGrouper) and self.grouper.labels is None:
                raise ValueError(
                    "Please pass `labels` to UniqueGrouper when grouping by a chunked array."
                )
            # this requires a pass to compute the bin edges
            if isinstance(self.grouper, BinGrouper) and isinstance(
                self.grouper.bins, int
            ):
                raise ValueError(
                    "Please pass explicit bin edges to BinGrouper using the ``bins`` kwarg"
                    "when grouping by a chunked array."
                )

        self.encoded = self.grouper.factorize(self.group)

    @property
    def name(self) -> Hashable:
        """Name for the grouped coordinate after reduction."""
        # the name has to come from unique_coord because we need `_bins` suffix for BinGrouper
        (name,) = self.encoded.unique_coord.dims
        return name

    @property
    def size(self) -> int:
        """Number of groups."""
        return len(self)

    def __len__(self) -> int:
        """Number of groups."""
        return len(self.encoded.full_index)


def _parse_group_and_groupers(
    obj: T_Xarray,
    group: GroupInput,
    groupers: dict[str, Grouper],
    *,
    eagerly_compute_group: Literal[False] | None,
) -> tuple[ResolvedGrouper, ...]:
    from xarray.core.dataarray import DataArray
    from xarray.core.variable import Variable
    from xarray.groupers import Grouper, UniqueGrouper

    if group is not None and groupers:
        raise ValueError(
            "Providing a combination of `group` and **groupers is not supported."
        )

    if group is None and not groupers:
        raise ValueError("Either `group` or `**groupers` must be provided.")

    if isinstance(group, np.ndarray | pd.Index):
        raise TypeError(
            f"`group` must be a DataArray. Received {type(group).__name__!r} instead"
        )

    if isinstance(group, Grouper):
        raise TypeError(
            "Cannot group by a Grouper object. "
            f"Instead use `.groupby(var_name={type(group).__name__}(...))`. "
            "You may need to assign the variable you're grouping by as a coordinate using `assign_coords`."
        )

    if isinstance(group, Mapping):
        grouper_mapping = either_dict_or_kwargs(group, groupers, "groupby")
        group = None

    rgroupers: tuple[ResolvedGrouper, ...]
    if isinstance(group, DataArray | Variable):
        rgroupers = (
            ResolvedGrouper(
                UniqueGrouper(), group, obj, eagerly_compute_group=eagerly_compute_group
            ),
        )
    else:
        if group is not None:
            if TYPE_CHECKING:
                assert isinstance(group, str | Sequence)
            group_iter: Sequence[Hashable] = (
                (group,) if isinstance(group, str) else group
            )
            grouper_mapping = {g: UniqueGrouper() for g in group_iter}
        elif groupers:
            grouper_mapping = cast("Mapping[Hashable, Grouper]", groupers)

        rgroupers = tuple(
            ResolvedGrouper(
                grouper, group, obj, eagerly_compute_group=eagerly_compute_group
            )
            for group, grouper in grouper_mapping.items()
        )
    return rgroupers


def _validate_groupby_squeeze(squeeze: Literal[False]) -> None:
    # While we don't generally check the type of every arg, passing
    # multiple dimensions as multiple arguments is common enough, and the
    # consequences hidden enough (strings evaluate as true) to warrant
    # checking here.
    # A future version could make squeeze kwarg only, but would face
    # backward-compat issues.
    if squeeze is not False:
        raise TypeError(f"`squeeze` must be False, but {squeeze!r} was supplied.")


def _resolve_group(
    obj: T_DataWithCoords, group: T_Group | Hashable | IndexVariable
) -> T_Group:
    from xarray.core.dataarray import DataArray

    error_msg = (
        "the group variable's length does not "
        "match the length of this variable along its "
        "dimensions"
    )

    newgroup: T_Group
    if isinstance(group, DataArray):
        try:
            align(obj, group, join="exact", copy=False)
        except ValueError as err:
            raise ValueError(error_msg) from err

        newgroup = group.copy(deep=False)
        newgroup.name = group.name or "group"

    elif isinstance(group, IndexVariable):
        # This assumption is built in to _ensure_1d.
        if group.ndim != 1:
            raise ValueError(
                "Grouping by multi-dimensional IndexVariables is not allowed."
                "Convert to and pass a DataArray instead."
            )
        (group_dim,) = group.dims
        if len(group) != obj.sizes[group_dim]:
            raise ValueError(error_msg)
        newgroup = DataArray(group)

    else:
        if not hashable(group):
            raise TypeError(
                "`group` must be an xarray.DataArray or the "
                "name of an xarray variable or dimension. "
                f"Received {group!r} instead."
            )
        group_da: DataArray = obj[group]
        if group_da.name not in obj._indexes and group_da.name in obj.dims:
            # DummyGroups should not appear on groupby results
            newgroup = _DummyGroup(obj, group_da.name, group_da.coords)
        else:
            newgroup = group_da

    if newgroup.size == 0:
        raise ValueError(f"{newgroup.name} must not be empty")

    return newgroup


@dataclass
class ComposedGrouper:
    """
    Helper class for multi-variable GroupBy.
    This satisfies the Grouper interface, but is awkward to wrap in ResolvedGrouper.
    For one, it simply re-infers a new EncodedGroups using known information
    in existing ResolvedGroupers. So passing in a `group` (hard to define),
    and `obj` (pointless) is not useful.
    """

    groupers: tuple[ResolvedGrouper, ...]

    def factorize(self) -> EncodedGroups:
        from xarray.groupers import EncodedGroups

        groupers = self.groupers

        # At this point all arrays have been factorized.
        codes = tuple(grouper.codes for grouper in groupers)
        shape = tuple(grouper.size for grouper in groupers)
        masks = tuple((code == -1) for code in codes)
        # We broadcast the codes against each other
        broadcasted_codes = broadcast(*codes)
        # This fully broadcasted DataArray is used as a template later
        first_codes = broadcasted_codes[0]
        # Now we convert to a single variable GroupBy problem
        _flatcodes = np.ravel_multi_index(
            tuple(codes.data for codes in broadcasted_codes), shape, mode="wrap"
        )
        # NaNs; as well as values outside the bins are coded by -1
        # Restore these after the raveling
        broadcasted_masks = broadcast(*masks)
        mask = functools.reduce(np.logical_or, broadcasted_masks)  # type: ignore[arg-type]
        _flatcodes = where(mask.data, -1, _flatcodes)

        full_index = pd.MultiIndex.from_product(
            [list(grouper.full_index.values) for grouper in groupers],
            names=tuple(grouper.name for grouper in groupers),
        )
        if not full_index.is_unique:
            raise ValueError(
                "The output index for the GroupBy is non-unique. "
                "This is a bug in the Grouper provided."
            )
        # This will be unused when grouping by dask arrays, so skip..
        if not is_chunked_array(_flatcodes):
            # Constructing an index from the product is wrong when there are missing groups
            # (e.g. binning, resampling). Account for that now.
            midx = full_index[np.sort(pd.unique(_flatcodes[~mask]))]
            group_indices = _codes_to_group_indices(_flatcodes.ravel(), len(full_index))
        else:
            midx = full_index
            group_indices = None

        dim_name = "stacked_" + "_".join(str(grouper.name) for grouper in groupers)

        coords = Coordinates.from_pandas_multiindex(midx, dim=dim_name)
        for grouper in groupers:
            coords.variables[grouper.name].attrs = grouper.group.attrs
        return EncodedGroups(
            codes=first_codes.copy(data=_flatcodes),
            full_index=full_index,
            group_indices=group_indices,
            unique_coord=Variable(dims=(dim_name,), data=midx.values),
            coords=coords,
        )


class GroupBy(Generic[T_Xarray]):
    """A object that implements the split-apply-combine pattern.

    Modeled after `pandas.GroupBy`. The `GroupBy` object can be iterated over
    (unique_value, grouped_array) pairs, but the main way to interact with a
    groupby object are with the `apply` or `reduce` methods. You can also
    directly call numpy methods like `mean` or `std`.

    You should create a GroupBy object by using the `DataArray.groupby` or
    `Dataset.groupby` methods.

    See Also
    --------
    Dataset.groupby
    DataArray.groupby
    """

    __slots__ = (
        "_by_chunked",
        "_codes",
        "_dims",
        "_group_dim",
        # cached properties
        "_groups",
        "_inserted_dims",
        "_len",
        "_obj",
        # Save unstacked object for flox
        "_original_obj",
        "_restore_coord_dims",
        "_sizes",
        "_stacked_dim",
        "encoded",
        # stack nD vars
        "group1d",
        "groupers",
    )
    _obj: T_Xarray
    groupers: tuple[ResolvedGrouper, ...]
    _restore_coord_dims: bool

    _original_obj: T_Xarray
    _group_indices: GroupIndices
    _codes: tuple[DataArray, ...]
    _group_dim: Hashable
    _by_chunked: bool

    _groups: dict[GroupKey, GroupIndex] | None
    _dims: tuple[Hashable, ...] | Frozen[Hashable, int] | None
    _sizes: Mapping[Hashable, int] | None
    _len: int

    # _ensure_1d:
    group1d: T_Group
    _stacked_dim: Hashable | None
    _inserted_dims: list[Hashable]

    encoded: EncodedGroups

    def __init__(
        self,
        obj: T_Xarray,
        groupers: tuple[ResolvedGrouper, ...],
        restore_coord_dims: bool = True,
    ) -> None:
        """Create a GroupBy object

        Parameters
        ----------
        obj : Dataset or DataArray
            Object to group.
        grouper : Grouper
            Grouper object
        restore_coord_dims : bool, default: True
            If True, also restore the dimension order of multi-dimensional
            coordinates.
        """
        self._original_obj = obj
        self._restore_coord_dims = restore_coord_dims
        self.groupers = groupers

        if len(groupers) == 1:
            (grouper,) = groupers
            self.encoded = grouper.encoded
        else:
            if any(
                isinstance(obj._indexes.get(grouper.name, None), PandasMultiIndex)
                for grouper in groupers
            ):
                raise NotImplementedError(
                    "Grouping by multiple variables, one of which "
                    "wraps a Pandas MultiIndex, is not supported yet."
                )
            self.encoded = ComposedGrouper(groupers).factorize()

        # specification for the groupby operation
        # TODO: handle obj having variables that are not present on any of the groupers
        #       simple broadcasting fails for ExtensionArrays.
        codes = self.encoded.codes
        self._by_chunked = is_chunked_array(codes._variable._data)
        if not self._by_chunked:
            (self.group1d, self._obj, self._stacked_dim, self._inserted_dims) = (
                _ensure_1d(group=codes, obj=obj)
            )
            (self._group_dim,) = self.group1d.dims
        else:
            self.group1d = None
            # This transpose preserves dim order behaviour
            self._obj = obj.transpose(..., *codes.dims)
            self._stacked_dim = None
            self._inserted_dims = []
            self._group_dim = None

        # cached attributes
        self._groups = None
        self._dims = None
        self._sizes = None
        self._len = len(self.encoded.full_index)

    @property
    def sizes(self) -> Mapping[Hashable, int]:
        """Ordered mapping from dimension names to lengths.

        Immutable.

        See Also
        --------
        DataArray.sizes
        Dataset.sizes
        """
        if self._sizes is None:
            index = self.encoded.group_indices[0]
            self._sizes = self._obj.isel({self._group_dim: index}).sizes
        return self._sizes

    def shuffle_to_chunks(self, chunks: T_Chunks = None) -> T_Xarray:
        """
        Sort or "shuffle" the underlying object.

        "Shuffle" means the object is sorted so that all group members occur sequentially,
        in the same chunk. Multiple groups may occur in the same chunk.
        This method is particularly useful for chunked arrays (e.g. dask, cubed).
        particularly when you need to map a function that requires all members of a group
        to be present in a single chunk. For chunked array types, the order of appearance
        is not guaranteed, but will depend on the input chunking.

        Parameters
        ----------
        chunks : int, tuple of int, "auto" or mapping of hashable to int or tuple of int, optional
            How to adjust chunks along dimensions not present in the array being grouped by.

        Returns
        -------
        DataArrayGroupBy or DatasetGroupBy

        Examples
        --------
        >>> import dask.array
        >>> da = xr.DataArray(
        ...     dims="x",
        ...     data=dask.array.arange(10, chunks=3),
        ...     coords={"x": [1, 2, 3, 1, 2, 3, 1, 2, 3, 0]},
        ...     name="a",
        ... )
        >>> shuffled = da.groupby("x").shuffle_to_chunks()
        >>> shuffled
        <xarray.DataArray 'a' (x: 10)> Size: 80B
        dask.array<shuffle, shape=(10,), dtype=int64, chunksize=(3,), chunktype=numpy.ndarray>
        Coordinates:
          * x        (x) int64 80B 0 1 1 1 2 2 2 3 3 3

        >>> shuffled.groupby("x").quantile(q=0.5).compute()
        <xarray.DataArray 'a' (x: 4)> Size: 32B
        array([9., 3., 4., 5.])
        Coordinates:
            quantile  float64 8B 0.5
          * x         (x) int64 32B 0 1 2 3

        See Also
        --------
        dask.dataframe.DataFrame.shuffle
        dask.array.shuffle
        """
        self._raise_if_by_is_chunked()
        return self._shuffle_obj(chunks)

    def _shuffle_obj(self, chunks: T_Chunks) -> T_Xarray:
        from xarray.core.dataarray import DataArray

        was_array = isinstance(self._obj, DataArray)
        as_dataset = self._obj._to_temp_dataset() if was_array else self._obj

        for grouper in self.groupers:
            if grouper.name not in as_dataset._variables:
                as_dataset.coords[grouper.name] = grouper.group

        shuffled = as_dataset._shuffle(
            dim=self._group_dim, indices=self.encoded.group_indices, chunks=chunks
        )
        unstacked: Dataset = self._maybe_unstack(shuffled)
        if was_array:
            return self._obj._from_temp_dataset(unstacked)
        else:
            return unstacked  # type: ignore[return-value]

    def map(
        self,
        func: Callable,
        args: tuple[Any, ...] = (),
        shortcut: bool | None = None,
        **kwargs: Any,
    ) -> T_Xarray:
        raise NotImplementedError()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        shortcut: bool = True,
        **kwargs: Any,
    ) -> T_Xarray:
        raise NotImplementedError()

    def _raise_if_by_is_chunked(self):
        if self._by_chunked:
            raise ValueError(
                "This method is not supported when lazily grouping by a chunked array. "
                "Either load the array in to memory prior to grouping using .load or .compute, "
                " or explore another way of applying your function, "
                "potentially using the `flox` package."
            )

    def _raise_if_not_single_group(self):
        if len(self.groupers) != 1:
            raise NotImplementedError(
                "This method is not supported for grouping by multiple variables yet."
            )

    @property
    def groups(self) -> dict[GroupKey, GroupIndex]:
        """
        Mapping from group labels to indices. The indices can be used to index the underlying object.
        """
        # provided to mimic pandas.groupby
        if self._groups is None:
            self._groups = dict(
                zip(
                    self.encoded.unique_coord.data,
                    self.encoded.group_indices,
                    strict=True,
                )
            )
        return self._groups

    def __getitem__(self, key: GroupKey) -> T_Xarray:
        """
        Get DataArray or Dataset corresponding to a particular group label.
        """
        self._raise_if_by_is_chunked()
        return self._obj.isel({self._group_dim: self.groups[key]})

    def __len__(self) -> int:
        return self._len

    def __iter__(self) -> Iterator[tuple[GroupKey, T_Xarray]]:
        return zip(self.encoded.unique_coord.data, self._iter_grouped(), strict=True)

    def __repr__(self) -> str:
        text = (
            f"<{self.__class__.__name__}, "
            f"grouped over {len(self.groupers)} grouper(s),"
            f" {self._len} groups in total:"
        )
        for grouper in self.groupers:
            coord = grouper.unique_coord
            labels = ", ".join(format_array_flat(coord, 30).split())
            text += (
                f"\n    {grouper.name!r}: {type(grouper.grouper).__name__}({grouper.group.name!r}), "
                f"{coord.size}/{grouper.full_index.size} groups with labels {labels}"
            )
        return text + ">"

    def _iter_grouped(self) -> Iterator[T_Xarray]:
        """Iterate over each element in this group"""
        self._raise_if_by_is_chunked()
        for indices in self.encoded.group_indices:
            if indices:
                yield self._obj.isel({self._group_dim: indices})

    def _infer_concat_args(self, applied_example):
        if self._group_dim in applied_example.dims:
            coord = self.group1d
            positions = self.encoded.group_indices
        else:
            coord = self.encoded.unique_coord
            positions = None
        (dim,) = coord.dims
        return dim, positions

    def _binary_op(self, other, f, reflexive=False):
        from xarray.core.dataarray import DataArray
        from xarray.core.dataset import Dataset

        g = f if not reflexive else lambda x, y: f(y, x)

        self._raise_if_not_single_group()
        (grouper,) = self.groupers
        obj = self._original_obj
        name = grouper.name
        group = grouper.group
        codes = self.encoded.codes
        dims = group.dims

        if isinstance(group, _DummyGroup):
            group = coord = group.to_dataarray()
        else:
            coord = grouper.unique_coord
            if isinstance(coord, Variable):
                assert coord.ndim == 1
                (coord_dim,) = coord.dims
                # TODO: explicitly create Index here
                coord = DataArray(coord, coords={coord_dim: coord.data})

        if not isinstance(other, Dataset | DataArray):
            raise TypeError(
                "GroupBy objects only support binary ops "
                "when the other argument is a Dataset or "
                "DataArray"
            )

        if name not in other.dims:
            raise ValueError(
                "incompatible dimensions for a grouped "
                f"binary operation: the group variable {name!r} "
                "is not a dimension on the other argument "
                f"with dimensions {other.dims!r}"
            )

        # Broadcast out scalars for backwards compatibility
        # TODO: get rid of this when fixing GH2145
        for var in other.coords:
            if other[var].ndim == 0:
                other[var] = (
                    other[var].drop_vars(var).expand_dims({name: other.sizes[name]})
                )

        # need to handle NaNs in group or elements that don't belong to any bins
        mask = codes == -1
        if mask.any():
            obj = obj.where(~mask, drop=True)
            group = group.where(~mask, drop=True)
            codes = codes.where(~mask, drop=True).astype(int)

        # if other is dask-backed, that's a hint that the
        # "expanded" dataset is too big to hold in memory.
        # this can be the case when `other` was read from disk
        # and contains our lazy indexing classes
        # We need to check for dask-backed Datasets
        # so utils.is_duck_dask_array does not work for this check
        if obj.chunks and not other.chunks:
            # TODO: What about datasets with some dask vars, and others not?
            # This handles dims other than `name``
            chunks = {k: v for k, v in obj.chunksizes.items() if k in other.dims}
            # a chunk size of 1 seems reasonable since we expect individual elements of
            # other to be repeated multiple times across the reduced dimension(s)
            chunks[name] = 1
            other = other.chunk(chunks)

        # codes are defined for coord, so we align `other` with `coord`
        # before indexing
        other, _ = align(other, coord, join="right", copy=False)
        expanded = other.isel({name: codes})

        result = g(obj, expanded)

        if group.ndim > 1:
            # backcompat:
            # TODO: get rid of this when fixing GH2145
            for var in set(obj.coords) - set(obj.xindexes):
                if set(obj[var].dims) < set(group.dims):
                    result[var] = obj[var].reset_coords(drop=True).broadcast_like(group)

        if isinstance(result, Dataset) and isinstance(obj, Dataset):
            for var in set(result):
                for d in dims:
                    if d not in obj[var].dims:
                        result[var] = result[var].transpose(d, ...)
        return result

    def _restore_dim_order(self, stacked):
        raise NotImplementedError

    def _maybe_reindex(self, combined):
        """Reindexing is needed in two cases:
        1. Our index contained empty groups (e.g., from a resampling or binning). If we
        reduced on that dimension, we want to restore the full index.

        2. We use a MultiIndex for multi-variable GroupBy.
        The MultiIndex stores each level's labels in sorted order
        which are then assigned on unstacking. So we need to restore
        the correct order here.
        """
        has_missing_groups = (
            self.encoded.unique_coord.size != self.encoded.full_index.size
        )
        indexers = {}
        for grouper in self.groupers:
            index = combined._indexes.get(grouper.name, None)
            if (has_missing_groups and index is not None) or (
                len(self.groupers) > 1
                and not isinstance(grouper.full_index, pd.RangeIndex)
                and not index.index.equals(grouper.full_index)
            ):
                indexers[grouper.name] = grouper.full_index
        if indexers:
            combined = combined.reindex(**indexers)
        return combined

    def _maybe_unstack(self, obj):
        """This gets called if we are applying on an array with a
        multidimensional group."""
        from xarray.groupers import UniqueGrouper

        stacked_dim = self._stacked_dim
        if stacked_dim is not None and stacked_dim in obj.dims:
            inserted_dims = self._inserted_dims
            obj = obj.unstack(stacked_dim)
            for dim in inserted_dims:
                if dim in obj.coords:
                    del obj.coords[dim]
            obj._indexes = filter_indexes_from_coords(obj._indexes, set(obj.coords))
        elif len(self.groupers) > 1:
            # TODO: we could clean this up by setting the appropriate `stacked_dim`
            #       and `inserted_dims`
            # if multiple groupers all share the same single dimension, then
            # we don't stack/unstack. Do that manually now.
            dims_to_unstack = self.encoded.unique_coord.dims
            if all(dim in obj.dims for dim in dims_to_unstack):
                obj = obj.unstack(*dims_to_unstack)
            to_drop = [
                grouper.name
                for grouper in self.groupers
                if isinstance(grouper.group, _DummyGroup)
                and isinstance(grouper.grouper, UniqueGrouper)
            ]
            obj = obj.drop_vars(to_drop)

        return obj

    def _flox_reduce(
        self,
        dim: Dims,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> T_Xarray:
        """Adaptor function that translates our groupby API to that of flox."""
        import flox
        from flox.xarray import xarray_reduce

        from xarray.core.dataset import Dataset

        obj = self._original_obj
        variables = (
            {k: v.variable for k, v in obj.data_vars.items()}
            if isinstance(obj, Dataset)  # type: ignore[redundant-expr]  # seems to be a mypy bug
            else obj._coords
        )

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        if Version(flox.__version__) < Version("0.9") and not self._by_chunked:
            # preserve current strategy (approximately) for dask groupby
            # on older flox versions to prevent surprises.
            # flox >=0.9 will choose this on its own.
            kwargs.setdefault("method", "cohorts")

        midx_grouping_vars: tuple[Hashable, ...] = ()
        for grouper in self.groupers:
            name = grouper.name
            maybe_midx = obj._indexes.get(name, None)
            if isinstance(maybe_midx, PandasMultiIndex):
                midx_grouping_vars += tuple(maybe_midx.index.names) + (name,)

        # For datasets, running a numeric-only reduction on non-numeric
        # variable will just drop it.
        non_numeric: dict[Hashable, Variable]
        if kwargs.pop("numeric_only", None):
            non_numeric = {
                name: var
                for name, var in variables.items()
                if (
                    not (np.issubdtype(var.dtype, np.number) or (var.dtype == np.bool_))
                    # this avoids dropping any levels of a MultiIndex, which raises
                    # a warning
                    and name not in midx_grouping_vars
                    and name not in obj.dims
                )
            }
        else:
            non_numeric = {}

        if "min_count" in kwargs:
            if kwargs["func"] not in ["sum", "prod"]:
                raise TypeError("Received an unexpected keyword argument 'min_count'")
            elif kwargs["min_count"] is None:
                # set explicitly to avoid unnecessarily accumulating count
                kwargs["min_count"] = 0

        parsed_dim: tuple[Hashable, ...]
        if isinstance(dim, str):
            parsed_dim = (dim,)
        elif dim is None:
            parsed_dim_list = list()
            # preserve order
            for dim_ in itertools.chain(
                *(grouper.codes.dims for grouper in self.groupers)
            ):
                if dim_ not in parsed_dim_list:
                    parsed_dim_list.append(dim_)
            parsed_dim = tuple(parsed_dim_list)
        elif dim is ...:
            parsed_dim = tuple(obj.dims)
        else:
            parsed_dim = tuple(dim)

        # Do this so we raise the same error message whether flox is present or not.
        # Better to control it here than in flox.
        for grouper in self.groupers:
            if any(
                d not in grouper.codes.dims and d not in obj.dims for d in parsed_dim
            ):
                raise ValueError(f"cannot reduce over dimensions {dim}.")

        has_missing_groups = (
            self.encoded.unique_coord.size != self.encoded.full_index.size
        )
        if self._by_chunked or has_missing_groups or kwargs.get("min_count", 0) > 0:
            # Xarray *always* returns np.nan when there are no observations in a group,
            # We can fake that here by forcing min_count=1 when it is not set.
            # This handles boolean reductions, and count
            # See GH8090, GH9398
            # Note that `has_missing_groups=False` when `self._by_chunked is True`.
            # We *choose* to always do the masking, so that behaviour is predictable
            # in some way. The real solution is to expose fill_value as a kwarg,
            # and set appropriate defaults :/.
            kwargs.setdefault("fill_value", np.nan)
            kwargs.setdefault("min_count", 1)

        # pass RangeIndex as a hint to flox that `by` is already factorized
        expected_groups = tuple(
            pd.RangeIndex(len(grouper)) for grouper in self.groupers
        )

        codes = tuple(g.codes for g in self.groupers)
        result = xarray_reduce(
            obj.drop_vars(non_numeric.keys()),
            *codes,
            dim=parsed_dim,
            expected_groups=expected_groups,
            isbin=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

        # we did end up reducing over dimension(s) that are
        # in the grouped variable
        group_dims = set(grouper.group.dims)
        new_coords = []
        to_drop = []
        if group_dims & set(parsed_dim):
            for grouper in self.groupers:
                output_index = grouper.full_index
                if isinstance(output_index, pd.RangeIndex):
                    # flox always assigns an index so we must drop it here if we don't need it.
                    to_drop.append(grouper.name)
                    continue
                # TODO: We can't simply use `self.encoded.coords` here because it corresponds to `unique_coord`,
                # NOT `full_index`. We would need to construct a new Coordinates object, that corresponds to `full_index`.
                new_coords.append(
                    # Using IndexVariable here ensures we reconstruct PandasMultiIndex with
                    # all associated levels properly.
                    coordinates_from_variable(
                        IndexVariable(
                            dims=grouper.name,
                            data=output_index,
                            attrs=grouper.codes.attrs,
                        )
                    )
                )
            result = result.assign_coords(
                Coordinates._construct_direct(*merge_coords(new_coords))
            ).drop_vars(to_drop)

        # broadcast any non-dim coord variables that don't
        # share all dimensions with the grouper
        result_variables = (
            result._variables if isinstance(result, Dataset) else result._coords
        )
        to_broadcast: dict[Hashable, Variable] = {}
        for name, var in variables.items():
            dims_set = set(var.dims)
            if (
                dims_set <= set(parsed_dim)
                and (dims_set & set(result.dims))
                and name not in result_variables
            ):
                to_broadcast[name] = var
        for name, var in to_broadcast.items():
            if new_dims := tuple(d for d in parsed_dim if d not in var.dims):
                new_sizes = tuple(
                    result.sizes.get(dim, obj.sizes.get(dim)) for dim in new_dims
                )
                result[name] = var.set_dims(
                    new_dims + var.dims, new_sizes + var.shape
                ).transpose(..., *result.dims)

        if not isinstance(result, Dataset):
            # only restore dimension order for arrays
            result = self._restore_dim_order(result)

        return result

    def fillna(self, value: Any) -> T_Xarray:
        """Fill missing values in this object by group.

        This operation follows the normal broadcasting and alignment rules that
        xarray uses for binary arithmetic, except the result is aligned to this
        object (``join='left'``) instead of aligned to the intersection of
        index coordinates (``join='inner'``).

        Parameters
        ----------
        value
            Used to fill all matching missing values by group. Needs
            to be of a valid type for the wrapped object's fillna
            method.

        Returns
        -------
        same type as the grouped object

        See Also
        --------
        Dataset.fillna
        DataArray.fillna
        """
        return ops.fillna(self, value)

    def quantile(
        self,
        q: ArrayLike,
        dim: Dims = None,
        *,
        method: QuantileMethods = "linear",
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
        interpolation: QuantileMethods | None = None,
    ) -> T_Xarray:
        """Compute the qth quantile over each array in the groups and
        concatenate them together into a new array.

        Parameters
        ----------
        q : float or sequence of float
            Quantile to compute, which must be between 0 and 1
            inclusive.
        dim : str or Iterable of Hashable, optional
            Dimension(s) over which to apply quantile.
            Defaults to the grouped dimension.
        method : str, default: "linear"
            This optional parameter specifies the interpolation method to use when the
            desired quantile lies between two data points. The options sorted by their R
            type as summarized in the H&F paper [1]_ are:

                1. "inverted_cdf"
                2. "averaged_inverted_cdf"
                3. "closest_observation"
                4. "interpolated_inverted_cdf"
                5. "hazen"
                6. "weibull"
                7. "linear"  (default)
                8. "median_unbiased"
                9. "normal_unbiased"

            The first three methods are discontiuous.  The following discontinuous
            variations of the default "linear" (7.) option are also available:

                * "lower"
                * "higher"
                * "midpoint"
                * "nearest"

            See :py:func:`numpy.quantile` or [1]_ for details. The "method" argument
            was previously called "interpolation", renamed in accordance with numpy
            version 1.22.0.
        keep_attrs : bool or None, default: None
            If True, the dataarray's attributes (`attrs`) will be copied from
            the original object to the new one.  If False, the new
            object will be returned without attributes.
        skipna : bool or None, default: None
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        quantiles : Variable
            If `q` is a single quantile, then the result is a
            scalar. If multiple percentiles are given, first axis of
            the result corresponds to the quantile. In either case a
            quantile dimension is added to the return array. The other
            dimensions are the dimensions that remain after the
            reduction of the array.

        See Also
        --------
        numpy.nanquantile, numpy.quantile, pandas.Series.quantile, Dataset.quantile
        DataArray.quantile

        Examples
        --------
        >>> da = xr.DataArray(
        ...     [[1.3, 8.4, 0.7, 6.9], [0.7, 4.2, 9.4, 1.5], [6.5, 7.3, 2.6, 1.9]],
        ...     coords={"x": [0, 0, 1], "y": [1, 1, 2, 2]},
        ...     dims=("x", "y"),
        ... )
        >>> ds = xr.Dataset({"a": da})
        >>> da.groupby("x").quantile(0)
        <xarray.DataArray (x: 2, y: 4)> Size: 64B
        array([[0.7, 4.2, 0.7, 1.5],
               [6.5, 7.3, 2.6, 1.9]])
        Coordinates:
          * y         (y) int64 32B 1 1 2 2
            quantile  float64 8B 0.0
          * x         (x) int64 16B 0 1
        >>> ds.groupby("y").quantile(0, dim=...)
        <xarray.Dataset> Size: 40B
        Dimensions:   (y: 2)
        Coordinates:
            quantile  float64 8B 0.0
          * y         (y) int64 16B 1 2
        Data variables:
            a         (y) float64 16B 0.7 0.7
        >>> da.groupby("x").quantile([0, 0.5, 1])
        <xarray.DataArray (x: 2, y: 4, quantile: 3)> Size: 192B
        array([[[0.7 , 1.  , 1.3 ],
                [4.2 , 6.3 , 8.4 ],
                [0.7 , 5.05, 9.4 ],
                [1.5 , 4.2 , 6.9 ]],
        <BLANKLINE>
               [[6.5 , 6.5 , 6.5 ],
                [7.3 , 7.3 , 7.3 ],
                [2.6 , 2.6 , 2.6 ],
                [1.9 , 1.9 , 1.9 ]]])
        Coordinates:
          * y         (y) int64 32B 1 1 2 2
          * quantile  (quantile) float64 24B 0.0 0.5 1.0
          * x         (x) int64 16B 0 1
        >>> ds.groupby("y").quantile([0, 0.5, 1], dim=...)
        <xarray.Dataset> Size: 88B
        Dimensions:   (y: 2, quantile: 3)
        Coordinates:
          * quantile  (quantile) float64 24B 0.0 0.5 1.0
          * y         (y) int64 16B 1 2
        Data variables:
            a         (y, quantile) float64 48B 0.7 5.35 8.4 0.7 2.25 9.4

        References
        ----------
        .. [1] R. J. Hyndman and Y. Fan,
           "Sample quantiles in statistical packages,"
           The American Statistician, 50(4), pp. 361-365, 1996
        """
        # Dataset.quantile does this, do it for flox to ensure same output.
        q = np.asarray(q, dtype=np.float64)

        if (
            method == "linear"
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
            and module_available("flox", minversion="0.9.4")
        ):
            result = self._flox_reduce(
                func="quantile", q=q, dim=dim, keep_attrs=keep_attrs, skipna=skipna
            )
            return result
        else:
            if dim is None:
                dim = (self._group_dim,)
            return self.map(
                self._obj.__class__.quantile,
                shortcut=False,
                q=q,
                dim=dim or self._group_dim,
                method=method,
                keep_attrs=keep_attrs,
                skipna=skipna,
                interpolation=interpolation,
            )

    def where(self, cond, other=dtypes.NA) -> T_Xarray:
        """Return elements from `self` or `other` depending on `cond`.

        Parameters
        ----------
        cond : DataArray or Dataset
            Locations at which to preserve this objects values. dtypes have to be `bool`
        other : scalar, DataArray or Dataset, optional
            Value to use for locations in this object where ``cond`` is False.
            By default, inserts missing values.

        Returns
        -------
        same type as the grouped object

        See Also
        --------
        Dataset.where
        """
        return ops.where_method(self, cond, other)

    def _first_or_last(
        self,
        op: Literal["first" | "last"],
        skipna: bool | None,
        keep_attrs: bool | None,
    ):
        if all(
            isinstance(maybe_slice, slice)
            and (maybe_slice.stop == maybe_slice.start + 1)
            for maybe_slice in self.encoded.group_indices
        ):
            # NB. this is currently only used for reductions along an existing
            # dimension
            return self._obj
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)
        if (
            module_available("flox", minversion="0.10.0")
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            import flox.xrdtypes

            result = self._flox_reduce(
                dim=None,
                func=op,
                skipna=skipna,
                keep_attrs=keep_attrs,
                fill_value=flox.xrdtypes.NA,
            )
        else:
            result = self.reduce(
                getattr(duck_array_ops, op),
                dim=[self._group_dim],
                skipna=skipna,
                keep_attrs=keep_attrs,
            )
        return result

    def first(
        self, skipna: bool | None = None, keep_attrs: bool | None = None
    ) -> T_Xarray:
        """
        Return the first element of each group along the group dimension

        Parameters
        ----------
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.

        """
        return self._first_or_last("first", skipna, keep_attrs)

    def last(
        self, skipna: bool | None = None, keep_attrs: bool | None = None
    ) -> T_Xarray:
        """
        Return the last element of each group along the group dimension

        Parameters
        ----------
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.


        """
        return self._first_or_last("last", skipna, keep_attrs)

    def assign_coords(self, coords=None, **coords_kwargs):
        """Assign coordinates by group.

        See Also
        --------
        Dataset.assign_coords
        Dataset.swap_dims
        """
        coords_kwargs = either_dict_or_kwargs(coords, coords_kwargs, "assign_coords")
        return self.map(lambda ds: ds.assign_coords(**coords_kwargs))


def _maybe_reorder(xarray_obj, dim, positions, N: int | None):
    order = _inverse_permutation_indices(positions, N)

    if order is None or len(order) != xarray_obj.sizes[dim]:
        return xarray_obj
    else:
        return xarray_obj[{dim: order}]


class DataArrayGroupByBase(GroupBy["DataArray"], DataArrayGroupbyArithmetic):
    """GroupBy object specialized to grouping DataArray objects"""

    __slots__ = ()
    _dims: tuple[Hashable, ...] | None

    @property
    def dims(self) -> tuple[Hashable, ...]:
        self._raise_if_by_is_chunked()
        if self._dims is None:
            index = self.encoded.group_indices[0]
            self._dims = self._obj.isel({self._group_dim: index}).dims
        return self._dims

    def _iter_grouped_shortcut(self):
        """Fast version of `_iter_grouped` that yields Variables without
        metadata
        """
        self._raise_if_by_is_chunked()
        var = self._obj.variable
        for _idx, indices in enumerate(self.encoded.group_indices):
            if indices:
                yield var[{self._group_dim: indices}]

    def _concat_shortcut(self, applied, dim, positions=None):
        # nb. don't worry too much about maintaining this method -- it does
        # speed things up, but it's not very interpretable and there are much
        # faster alternatives (e.g., doing the grouped aggregation in a
        # compiled language)
        # TODO: benbovy - explicit indexes: this fast implementation doesn't
        # create an explicit index for the stacked dim coordinate
        stacked = Variable.concat(applied, dim, shortcut=True)
        reordered = _maybe_reorder(stacked, dim, positions, N=self.group1d.size)
        return self._obj._replace_maybe_drop_dims(reordered)

    def _restore_dim_order(self, stacked: DataArray) -> DataArray:
        def lookup_order(dimension):
            for grouper in self.groupers:
                if dimension == grouper.name and grouper.group.ndim == 1:
                    (dimension,) = grouper.group.dims
            if dimension in self._obj.dims:
                axis = self._obj.get_axis_num(dimension)
            else:
                axis = 1e6  # some arbitrarily high value
            return axis

        new_order = sorted(stacked.dims, key=lookup_order)
        stacked = stacked.transpose(
            *new_order, transpose_coords=self._restore_coord_dims
        )
        return stacked

    def map(
        self,
        func: Callable[..., DataArray],
        args: tuple[Any, ...] = (),
        shortcut: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """Apply a function to each array in the group and concatenate them
        together into a new array.

        `func` is called like `func(ar, *args, **kwargs)` for each array `ar`
        in this group.

        Apply uses heuristics (like `pandas.GroupBy.apply`) to figure out how
        to stack together the array. The rule is:

        1. If the dimension along which the group coordinate is defined is
           still in the first grouped array after applying `func`, then stack
           over this dimension.
        2. Otherwise, stack over the new dimension given by name of this
           grouping (the argument to the `groupby` function).

        Parameters
        ----------
        func : callable
            Callable to apply to each array.
        shortcut : bool, optional
            Whether or not to shortcut evaluation under the assumptions that:

            (1) The action of `func` does not depend on any of the array
                metadata (attributes or coordinates) but only on the data and
                dimensions.
            (2) The action of `func` creates arrays with homogeneous metadata,
                that is, with the same dimensions and attributes.

            If these conditions are satisfied `shortcut` provides significant
            speedup. This should be the case for many common groupby operations
            (e.g., applying numpy ufuncs).
        *args : tuple, optional
            Positional arguments passed to `func`.
        **kwargs
            Used to call `func(ar, **kwargs)` for each array `ar`.

        Returns
        -------
        applied : DataArray
            The result of splitting, applying and combining this array.
        """
        grouped = self._iter_grouped_shortcut() if shortcut else self._iter_grouped()
        applied = (maybe_wrap_array(arr, func(arr, *args, **kwargs)) for arr in grouped)
        return self._combine(applied, shortcut=shortcut)

    def apply(self, func, shortcut=False, args=(), **kwargs):
        """
        Backward compatible implementation of ``map``

        See Also
        --------
        DataArrayGroupBy.map
        """
        warnings.warn(
            "GroupBy.apply may be deprecated in the future. Using GroupBy.map is encouraged",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        return self.map(func, shortcut=shortcut, args=args, **kwargs)

    def _combine(self, applied, shortcut=False):
        """Recombine the applied objects like the original."""
        applied_example, applied = peek_at(applied)
        dim, positions = self._infer_concat_args(applied_example)
        if shortcut:
            combined = self._concat_shortcut(applied, dim, positions)
        else:
            combined = concat(
                applied,
                dim,
                data_vars="all",
                coords="different",
                compat="equals",
                join="outer",
            )
            combined = _maybe_reorder(combined, dim, positions, N=self.group1d.size)

        if isinstance(combined, type(self._obj)):
            # only restore dimension order for arrays
            combined = self._restore_dim_order(combined)
        # assign coord and index when the applied function does not return that coord
        if dim not in applied_example.dims:
            combined = combined.assign_coords(self.encoded.coords)
        combined = self._maybe_unstack(combined)
        combined = self._maybe_reindex(combined)
        return combined

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        shortcut: bool = True,
        **kwargs: Any,
    ) -> DataArray:
        """Reduce the items in this group by applying `func` along some
        dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, axis=axis, **kwargs)` to return the result of collapsing
            an np.ndarray over an integer valued axis.
        dim : "...", str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`. If None, apply over the
            groupby dimension, if "..." apply over all dimensions.
        axis : int or sequence of int, optional
            Axis(es) over which to apply `func`. Only one of the 'dimension'
            and 'axis' arguments can be supplied. If neither are supplied, then
            `func` is calculated over all dimension for each group item.
        keep_attrs : bool, optional
            If True, the datasets's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : Array
            Array with summarized data and the indicated dimension(s)
            removed.
        """
        if self._by_chunked:
            raise ValueError(
                "This method is not supported when lazily grouping by a chunked array. "
                "Try installing the `flox` package if you are using one of the standard "
                "reductions (e.g. `mean`). "
            )
        if dim is None:
            dim = [self._group_dim]

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        def reduce_array(ar: DataArray) -> DataArray:
            return ar.reduce(
                func=func,
                dim=dim,
                axis=axis,
                keep_attrs=keep_attrs,
                keepdims=keepdims,
                **kwargs,
            )

        check_reduce_dims(dim, self.dims)

        return self.map(reduce_array, shortcut=shortcut)


class DataArrayGroupBy(
    DataArrayGroupByBase,
    DataArrayGroupByAggregations,
    ImplementsArrayReduce,
):
    __slots__ = ()


class DatasetGroupByBase(GroupBy["Dataset"], DatasetGroupbyArithmetic):
    __slots__ = ()
    _dims: Frozen[Hashable, int] | None

    @property
    def dims(self) -> Frozen[Hashable, int]:
        self._raise_if_by_is_chunked()
        if self._dims is None:
            index = self.encoded.group_indices[0]
            self._dims = self._obj.isel({self._group_dim: index}).dims

        return FrozenMappingWarningOnValuesAccess(self._dims)

    def map(
        self,
        func: Callable[..., Dataset],
        args: tuple[Any, ...] = (),
        shortcut: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """Apply a function to each Dataset in the group and concatenate them
        together into a new Dataset.

        `func` is called like `func(ds, *args, **kwargs)` for each dataset `ds`
        in this group.

        Apply uses heuristics (like `pandas.GroupBy.apply`) to figure out how
        to stack together the datasets. The rule is:

        1. If the dimension along which the group coordinate is defined is
           still in the first grouped item after applying `func`, then stack
           over this dimension.
        2. Otherwise, stack over the new dimension given by name of this
           grouping (the argument to the `groupby` function).

        Parameters
        ----------
        func : callable
            Callable to apply to each sub-dataset.
        args : tuple, optional
            Positional arguments to pass to `func`.
        **kwargs
            Used to call `func(ds, **kwargs)` for each sub-dataset `ar`.

        Returns
        -------
        applied : Dataset
            The result of splitting, applying and combining this dataset.
        """
        # ignore shortcut if set (for now)
        applied = (func(ds, *args, **kwargs) for ds in self._iter_grouped())
        return self._combine(applied)

    def apply(self, func, args=(), shortcut=None, **kwargs):
        """
        Backward compatible implementation of ``map``

        See Also
        --------
        DatasetGroupBy.map
        """

        warnings.warn(
            "GroupBy.apply may be deprecated in the future. Using GroupBy.map is encouraged",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        return self.map(func, shortcut=shortcut, args=args, **kwargs)

    def _combine(self, applied):
        """Recombine the applied objects like the original."""
        applied_example, applied = peek_at(applied)
        dim, positions = self._infer_concat_args(applied_example)
        combined = concat(
            applied,
            dim,
            data_vars="all",
            coords="different",
            compat="equals",
            join="outer",
        )
        combined = _maybe_reorder(combined, dim, positions, N=self.group1d.size)
        # assign coord when the applied function does not return that coord
        if dim not in applied_example.dims:
            combined = combined.assign_coords(self.encoded.coords)
        combined = self._maybe_unstack(combined)
        combined = self._maybe_reindex(combined)
        return combined

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        shortcut: bool = True,
        **kwargs: Any,
    ) -> Dataset:
        """Reduce the items in this group by applying `func` along some
        dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, axis=axis, **kwargs)` to return the result of collapsing
            an np.ndarray over an integer valued axis.
        dim : ..., str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`. By default apply over the
            groupby dimension, with "..." apply over all dimensions.
        axis : int or sequence of int, optional
            Axis(es) over which to apply `func`. Only one of the 'dimension'
            and 'axis' arguments can be supplied. If neither are supplied, then
            `func` is calculated over all dimension for each group item.
        keep_attrs : bool, optional
            If True, the datasets's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : Dataset
            Array with summarized data and the indicated dimension(s)
            removed.
        """

        if self._by_chunked:
            raise ValueError(
                "This method is not supported when lazily grouping by a chunked array. "
                "Try installing the `flox` package if you are using one of the standard "
                "reductions (e.g. `mean`). "
            )

        if dim is None:
            dim = [self._group_dim]

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        def reduce_dataset(ds: Dataset) -> Dataset:
            return ds.reduce(
                func=func,
                dim=dim,
                axis=axis,
                keep_attrs=keep_attrs,
                keepdims=keepdims,
                **kwargs,
            )

        check_reduce_dims(dim, self.dims)

        return self.map(reduce_dataset)

    def assign(self, **kwargs: Any) -> Dataset:
        """Assign data variables by group.

        See Also
        --------
        Dataset.assign
        """
        return self.map(lambda ds: ds.assign(**kwargs))


class DatasetGroupBy(
    DatasetGroupByBase,
    DatasetGroupByAggregations,
    ImplementsDatasetReduce,
):
    __slots__ = ()
