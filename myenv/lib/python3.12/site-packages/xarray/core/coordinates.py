from __future__ import annotations

from collections.abc import Hashable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    cast,
)

import numpy as np
import pandas as pd

from xarray.core import formatting
from xarray.core.alignment import Aligner
from xarray.core.indexes import (
    Index,
    Indexes,
    PandasIndex,
    PandasMultiIndex,
    assert_no_index_corrupted,
    create_default_index_implicit,
)
from xarray.core.merge import merge_coordinates_without_align, merge_coords
from xarray.core.types import DataVars, Self, T_DataArray, T_Xarray
from xarray.core.utils import (
    Frozen,
    ReprObject,
    either_dict_or_kwargs,
    emit_user_level_warning,
)
from xarray.core.variable import Variable, as_variable, calculate_dimensions

if TYPE_CHECKING:
    from xarray.core.common import DataWithCoords
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

# Used as the key corresponding to a DataArray's variable when converting
# arbitrary DataArray objects to datasets
_THIS_ARRAY = ReprObject("<this-array>")


class AbstractCoordinates(Mapping[Hashable, "T_DataArray"]):
    _data: DataWithCoords
    __slots__ = ("_data",)

    def __getitem__(self, key: Hashable) -> T_DataArray:
        raise NotImplementedError()

    @property
    def _names(self) -> set[Hashable]:
        raise NotImplementedError()

    @property
    def dims(self) -> Frozen[Hashable, int] | tuple[Hashable, ...]:
        raise NotImplementedError()

    @property
    def dtypes(self) -> Frozen[Hashable, np.dtype]:
        raise NotImplementedError()

    @property
    def indexes(self) -> Indexes[pd.Index]:
        """Mapping of pandas.Index objects used for label based indexing.

        Raises an error if this Coordinates object has indexes that cannot
        be coerced to pandas.Index objects.

        See Also
        --------
        Coordinates.xindexes
        """
        return self._data.indexes

    @property
    def xindexes(self) -> Indexes[Index]:
        """Mapping of :py:class:`~xarray.indexes.Index` objects
        used for label based indexing.
        """
        return self._data.xindexes

    @property
    def variables(self):
        raise NotImplementedError()

    def _update_coords(self, coords, indexes):
        raise NotImplementedError()

    def _drop_coords(self, coord_names):
        raise NotImplementedError()

    def __iter__(self) -> Iterator[Hashable]:
        # needs to be in the same order as the dataset variables
        for k in self.variables:
            if k in self._names:
                yield k

    def __len__(self) -> int:
        return len(self._names)

    def __contains__(self, key: Hashable) -> bool:
        return key in self._names

    def __repr__(self) -> str:
        return formatting.coords_repr(self)

    def to_dataset(self) -> Dataset:
        raise NotImplementedError()

    def to_index(self, ordered_dims: Sequence[Hashable] | None = None) -> pd.Index:
        """Convert all index coordinates into a :py:class:`pandas.Index`.

        Parameters
        ----------
        ordered_dims : sequence of hashable, optional
            Possibly reordered version of this object's dimensions indicating
            the order in which dimensions should appear on the result.

        Returns
        -------
        pandas.Index
            Index subclass corresponding to the outer-product of all dimension
            coordinates. This will be a MultiIndex if this object is has more
            than more dimension.
        """
        if ordered_dims is None:
            ordered_dims = list(self.dims)
        elif set(ordered_dims) != set(self.dims):
            raise ValueError(
                "ordered_dims must match dims, but does not: "
                f"{ordered_dims} vs {self.dims}"
            )

        if len(ordered_dims) == 0:
            raise ValueError("no valid index for a 0-dimensional object")
        elif len(ordered_dims) == 1:
            (dim,) = ordered_dims
            return self._data.get_index(dim)
        else:
            indexes = [self._data.get_index(k) for k in ordered_dims]

            # compute the sizes of the repeat and tile for the cartesian product
            # (taken from pandas.core.reshape.util)
            index_lengths = np.fromiter(
                (len(index) for index in indexes), dtype=np.intp
            )
            cumprod_lengths = np.cumprod(index_lengths)

            if cumprod_lengths[-1] == 0:
                # if any factor is empty, the cartesian product is empty
                repeat_counts = np.zeros_like(cumprod_lengths)

            else:
                # sizes of the repeats
                repeat_counts = cumprod_lengths[-1] / cumprod_lengths
            # sizes of the tiles
            tile_counts = np.roll(cumprod_lengths, 1)
            tile_counts[0] = 1

            # loop over the indexes
            # for each MultiIndex or Index compute the cartesian product of the codes

            code_list = []
            level_list = []
            names = []

            for i, index in enumerate(indexes):
                if isinstance(index, pd.MultiIndex):
                    codes, levels = index.codes, index.levels
                else:
                    code, level = pd.factorize(index)
                    codes = [code]
                    levels = [level]

                # compute the cartesian product
                code_list += [
                    np.tile(np.repeat(code, repeat_counts[i]), tile_counts[i])
                    for code in codes
                ]
                level_list += levels
                names += index.names

        return pd.MultiIndex(level_list, code_list, names=names)


class Coordinates(AbstractCoordinates):
    """Dictionary like container for Xarray coordinates (variables + indexes).

    This collection is a mapping of coordinate names to
    :py:class:`~xarray.DataArray` objects.

    It can be passed directly to the :py:class:`~xarray.Dataset` and
    :py:class:`~xarray.DataArray` constructors via their `coords` argument. This
    will add both the coordinates variables and their index.

    Coordinates are either:

    - returned via the :py:attr:`Dataset.coords` and :py:attr:`DataArray.coords`
      properties
    - built from Pandas or other index objects
      (e.g., :py:meth:`Coordinates.from_pandas_multiindex`)
    - built directly from coordinate data and Xarray ``Index`` objects (beware that
      no consistency check is done on those inputs)

    Parameters
    ----------
    coords: dict-like, optional
        Mapping where keys are coordinate names and values are objects that
        can be converted into a :py:class:`~xarray.Variable` object
        (see :py:func:`~xarray.as_variable`). If another
        :py:class:`~xarray.Coordinates` object is passed, its indexes
        will be added to the new created object.
    indexes: dict-like, optional
        Mapping where keys are coordinate names and values are
        :py:class:`~xarray.indexes.Index` objects. If None (default),
        pandas indexes will be created for each dimension coordinate.
        Passing an empty dictionary will skip this default behavior.

    Examples
    --------
    Create a dimension coordinate with a default (pandas) index:

    >>> xr.Coordinates({"x": [1, 2]})
    Coordinates:
      * x        (x) int64 16B 1 2

    Create a dimension coordinate with no index:

    >>> xr.Coordinates(coords={"x": [1, 2]}, indexes={})
    Coordinates:
        x        (x) int64 16B 1 2

    Create a new Coordinates object from existing dataset coordinates
    (indexes are passed):

    >>> ds = xr.Dataset(coords={"x": [1, 2]})
    >>> xr.Coordinates(ds.coords)
    Coordinates:
      * x        (x) int64 16B 1 2

    Create indexed coordinates from a ``pandas.MultiIndex`` object:

    >>> midx = pd.MultiIndex.from_product([["a", "b"], [0, 1]])
    >>> xr.Coordinates.from_pandas_multiindex(midx, "x")
    Coordinates:
      * x          (x) object 32B MultiIndex
      * x_level_0  (x) object 32B 'a' 'a' 'b' 'b'
      * x_level_1  (x) int64 32B 0 1 0 1

    Create a new Dataset object by passing a Coordinates object:

    >>> midx_coords = xr.Coordinates.from_pandas_multiindex(midx, "x")
    >>> xr.Dataset(coords=midx_coords)
    <xarray.Dataset> Size: 96B
    Dimensions:    (x: 4)
    Coordinates:
      * x          (x) object 32B MultiIndex
      * x_level_0  (x) object 32B 'a' 'a' 'b' 'b'
      * x_level_1  (x) int64 32B 0 1 0 1
    Data variables:
        *empty*

    """

    _data: DataWithCoords

    __slots__ = ("_data",)

    def __init__(
        self,
        coords: Mapping[Any, Any] | None = None,
        indexes: Mapping[Any, Index] | None = None,
    ) -> None:
        # When coordinates are constructed directly, an internal Dataset is
        # created so that it is compatible with the DatasetCoordinates and
        # DataArrayCoordinates classes serving as a proxy for the data.
        # TODO: refactor DataArray / Dataset so that Coordinates store the data.
        from xarray.core.dataset import Dataset

        if coords is None:
            coords = {}

        variables: dict[Hashable, Variable]
        default_indexes: dict[Hashable, PandasIndex] = {}
        coords_obj_indexes: dict[Hashable, Index] = {}

        if isinstance(coords, Coordinates):
            if indexes is not None:
                raise ValueError(
                    "passing both a ``Coordinates`` object and a mapping of indexes "
                    "to ``Coordinates.__init__`` is not allowed "
                    "(this constructor does not support merging them)"
                )
            variables = {k: v.copy() for k, v in coords.variables.items()}
            coords_obj_indexes = dict(coords.xindexes)
        else:
            variables = {}
            for name, data in coords.items():
                var = as_variable(data, name=name, auto_convert=False)
                if var.dims == (name,) and indexes is None:
                    index, index_vars = create_default_index_implicit(var, list(coords))
                    default_indexes.update({k: index for k in index_vars})
                    variables.update(index_vars)
                else:
                    variables[name] = var

        if indexes is None:
            indexes = {}
        else:
            indexes = dict(indexes)

        indexes.update(default_indexes)
        indexes.update(coords_obj_indexes)

        no_coord_index = set(indexes) - set(variables)
        if no_coord_index:
            raise ValueError(
                f"no coordinate variables found for these indexes: {no_coord_index}"
            )

        for k, idx in indexes.items():
            if not isinstance(idx, Index):
                raise TypeError(f"'{k}' is not an `xarray.indexes.Index` object")

        # maybe convert to base variable
        for k, v in variables.items():
            if k not in indexes:
                variables[k] = v.to_base_variable()

        self._data = Dataset._construct_direct(
            coord_names=set(variables), variables=variables, indexes=indexes
        )

    @classmethod
    def _construct_direct(
        cls,
        coords: dict[Any, Variable],
        indexes: dict[Any, Index],
        dims: dict[Any, int] | None = None,
    ) -> Self:
        from xarray.core.dataset import Dataset

        obj = object.__new__(cls)
        obj._data = Dataset._construct_direct(
            coord_names=set(coords),
            variables=coords,
            indexes=indexes,
            dims=dims,
        )
        return obj

    @classmethod
    def from_pandas_multiindex(cls, midx: pd.MultiIndex, dim: str) -> Self:
        """Wrap a pandas multi-index as Xarray coordinates (dimension + levels).

        The returned coordinates can be directly assigned to a
        :py:class:`~xarray.Dataset` or :py:class:`~xarray.DataArray` via the
        ``coords`` argument of their constructor.

        Parameters
        ----------
        midx : :py:class:`pandas.MultiIndex`
            Pandas multi-index object.
        dim : str
            Dimension name.

        Returns
        -------
        coords : Coordinates
            A collection of Xarray indexed coordinates created from the multi-index.

        """
        xr_idx = PandasMultiIndex(midx, dim)

        variables = xr_idx.create_variables()
        indexes = {k: xr_idx for k in variables}

        return cls(coords=variables, indexes=indexes)

    @property
    def _names(self) -> set[Hashable]:
        return self._data._coord_names

    @property
    def dims(self) -> Frozen[Hashable, int] | tuple[Hashable, ...]:
        """Mapping from dimension names to lengths or tuple of dimension names."""
        return self._data.dims

    @property
    def sizes(self) -> Frozen[Hashable, int]:
        """Mapping from dimension names to lengths."""
        return self._data.sizes

    @property
    def dtypes(self) -> Frozen[Hashable, np.dtype]:
        """Mapping from coordinate names to dtypes.

        Cannot be modified directly.

        See Also
        --------
        Dataset.dtypes
        """
        return Frozen({n: v.dtype for n, v in self._data.variables.items()})

    @property
    def variables(self) -> Mapping[Hashable, Variable]:
        """Low level interface to Coordinates contents as dict of Variable objects.

        This dictionary is frozen to prevent mutation.
        """
        return self._data.variables

    def to_dataset(self) -> Dataset:
        """Convert these coordinates into a new Dataset."""
        names = [name for name in self._data._variables if name in self._names]
        return self._data._copy_listed(names)

    def __getitem__(self, key: Hashable) -> DataArray:
        return self._data[key]

    def __delitem__(self, key: Hashable) -> None:
        # redirect to DatasetCoordinates.__delitem__
        del self._data.coords[key]

    def equals(self, other: Self) -> bool:
        """Two Coordinates objects are equal if they have matching variables,
        all of which are equal.

        See Also
        --------
        Coordinates.identical
        """
        if not isinstance(other, Coordinates):
            return False
        return self.to_dataset().equals(other.to_dataset())

    def identical(self, other: Self) -> bool:
        """Like equals, but also checks all variable attributes.

        See Also
        --------
        Coordinates.equals
        """
        if not isinstance(other, Coordinates):
            return False
        return self.to_dataset().identical(other.to_dataset())

    def _update_coords(
        self, coords: dict[Hashable, Variable], indexes: Mapping[Any, Index]
    ) -> None:
        # redirect to DatasetCoordinates._update_coords
        self._data.coords._update_coords(coords, indexes)

    def _drop_coords(self, coord_names):
        # redirect to DatasetCoordinates._drop_coords
        self._data.coords._drop_coords(coord_names)

    def _merge_raw(self, other, reflexive):
        """For use with binary arithmetic."""
        if other is None:
            variables = dict(self.variables)
            indexes = dict(self.xindexes)
        else:
            coord_list = [self, other] if not reflexive else [other, self]
            variables, indexes = merge_coordinates_without_align(coord_list)
        return variables, indexes

    @contextmanager
    def _merge_inplace(self, other):
        """For use with in-place binary arithmetic."""
        if other is None:
            yield
        else:
            # don't include indexes in prioritized, because we didn't align
            # first and we want indexes to be checked
            prioritized = {
                k: (v, None)
                for k, v in self.variables.items()
                if k not in self.xindexes
            }
            variables, indexes = merge_coordinates_without_align(
                [self, other], prioritized
            )
            yield
            self._update_coords(variables, indexes)

    def merge(self, other: Mapping[Any, Any] | None) -> Dataset:
        """Merge two sets of coordinates to create a new Dataset

        The method implements the logic used for joining coordinates in the
        result of a binary operation performed on xarray objects:

        - If two index coordinates conflict (are not equal), an exception is
          raised. You must align your data before passing it to this method.
        - If an index coordinate and a non-index coordinate conflict, the non-
          index coordinate is dropped.
        - If two non-index coordinates conflict, both are dropped.

        Parameters
        ----------
        other : dict-like, optional
            A :py:class:`Coordinates` object or any mapping that can be turned
            into coordinates.

        Returns
        -------
        merged : Dataset
            A new Dataset with merged coordinates.
        """
        from xarray.core.dataset import Dataset

        if other is None:
            return self.to_dataset()

        if not isinstance(other, Coordinates):
            other = Dataset(coords=other).coords

        coords, indexes = merge_coordinates_without_align([self, other])
        coord_names = set(coords)
        return Dataset._construct_direct(
            variables=coords, coord_names=coord_names, indexes=indexes
        )

    def __setitem__(self, key: Hashable, value: Any) -> None:
        self.update({key: value})

    def update(self, other: Mapping[Any, Any]) -> None:
        """Update this Coordinates variables with other coordinate variables."""

        if not len(other):
            return

        other_coords: Coordinates

        if isinstance(other, Coordinates):
            # Coordinates object: just pass it (default indexes won't be created)
            other_coords = other
        else:
            other_coords = create_coords_with_default_indexes(
                getattr(other, "variables", other)
            )

        # Discard original indexed coordinates prior to merge allows to:
        # - fail early if the new coordinates don't preserve the integrity of existing
        #   multi-coordinate indexes
        # - drop & replace coordinates without alignment (note: we must keep indexed
        #   coordinates extracted from the DataArray objects passed as values to
        #   `other` - if any - as those are still used for aligning the old/new coordinates)
        coords_to_align = drop_indexed_coords(set(other_coords) & set(other), self)

        coords, indexes = merge_coords(
            [coords_to_align, other_coords],
            priority_arg=1,
            indexes=coords_to_align.xindexes,
        )

        # special case for PandasMultiIndex: updating only its dimension coordinate
        # is still allowed but depreciated.
        # It is the only case where we need to actually drop coordinates here (multi-index levels)
        # TODO: remove when removing PandasMultiIndex's dimension coordinate.
        self._drop_coords(self._names - coords_to_align._names)

        self._update_coords(coords, indexes)

    def assign(self, coords: Mapping | None = None, **coords_kwargs: Any) -> Self:
        """Assign new coordinates (and indexes) to a Coordinates object, returning
        a new object with all the original coordinates in addition to the new ones.

        Parameters
        ----------
        coords : mapping of dim to coord, optional
            A mapping whose keys are the names of the coordinates and values are the
            coordinates to assign. The mapping will generally be a dict or
            :class:`Coordinates`.

            * If a value is a standard data value — for example, a ``DataArray``,
              scalar, or array — the data is simply assigned as a coordinate.

            * A coordinate can also be defined and attached to an existing dimension
              using a tuple with the first element the dimension name and the second
              element the values for this new coordinate.

        **coords_kwargs
            The keyword arguments form of ``coords``.
            One of ``coords`` or ``coords_kwargs`` must be provided.

        Returns
        -------
        new_coords : Coordinates
            A new Coordinates object with the new coordinates (and indexes)
            in addition to all the existing coordinates.

        Examples
        --------
        >>> coords = xr.Coordinates()
        >>> coords
        Coordinates:
            *empty*

        >>> coords.assign(x=[1, 2])
        Coordinates:
          * x        (x) int64 16B 1 2

        >>> midx = pd.MultiIndex.from_product([["a", "b"], [0, 1]])
        >>> coords.assign(xr.Coordinates.from_pandas_multiindex(midx, "y"))
        Coordinates:
          * y          (y) object 32B MultiIndex
          * y_level_0  (y) object 32B 'a' 'a' 'b' 'b'
          * y_level_1  (y) int64 32B 0 1 0 1

        """
        # TODO: this doesn't support a callable, which is inconsistent with `DataArray.assign_coords`
        coords = either_dict_or_kwargs(coords, coords_kwargs, "assign")
        new_coords = self.copy()
        new_coords.update(coords)
        return new_coords

    def _overwrite_indexes(
        self,
        indexes: Mapping[Any, Index],
        variables: Mapping[Any, Variable] | None = None,
    ) -> Self:
        results = self.to_dataset()._overwrite_indexes(indexes, variables)

        # TODO: remove cast once we get rid of DatasetCoordinates
        # and DataArrayCoordinates (i.e., Dataset and DataArray encapsulate Coordinates)
        return cast(Self, results.coords)

    def _reindex_callback(
        self,
        aligner: Aligner,
        dim_pos_indexers: dict[Hashable, Any],
        variables: dict[Hashable, Variable],
        indexes: dict[Hashable, Index],
        fill_value: Any,
        exclude_dims: frozenset[Hashable],
        exclude_vars: frozenset[Hashable],
    ) -> Self:
        """Callback called from ``Aligner`` to create a new reindexed Coordinate."""
        aligned = self.to_dataset()._reindex_callback(
            aligner,
            dim_pos_indexers,
            variables,
            indexes,
            fill_value,
            exclude_dims,
            exclude_vars,
        )

        # TODO: remove cast once we get rid of DatasetCoordinates
        # and DataArrayCoordinates (i.e., Dataset and DataArray encapsulate Coordinates)
        return cast(Self, aligned.coords)

    def _ipython_key_completions_(self):
        """Provide method for the key-autocompletions in IPython."""
        return self._data._ipython_key_completions_()

    def copy(
        self,
        deep: bool = False,
        memo: dict[int, Any] | None = None,
    ) -> Self:
        """Return a copy of this Coordinates object."""
        # do not copy indexes (may corrupt multi-coordinate indexes)
        # TODO: disable variables deepcopy? it may also be problematic when they
        # encapsulate index objects like pd.Index
        variables = {
            k: v._copy(deep=deep, memo=memo) for k, v in self.variables.items()
        }

        # TODO: getting an error with `self._construct_direct`, possibly because of how
        # a subclass implements `_construct_direct`. (This was originally the same
        # runtime code, but we switched the type definitions in #8216, which
        # necessitates the cast.)
        return cast(
            Self,
            Coordinates._construct_direct(
                coords=variables, indexes=dict(self.xindexes), dims=dict(self.sizes)
            ),
        )


class DatasetCoordinates(Coordinates):
    """Dictionary like container for Dataset coordinates (variables + indexes).

    This collection can be passed directly to the :py:class:`~xarray.Dataset`
    and :py:class:`~xarray.DataArray` constructors via their `coords` argument.
    This will add both the coordinates variables and their index.
    """

    _data: Dataset

    __slots__ = ("_data",)

    def __init__(self, dataset: Dataset):
        self._data = dataset

    @property
    def _names(self) -> set[Hashable]:
        return self._data._coord_names

    @property
    def dims(self) -> Frozen[Hashable, int]:
        return self._data.dims

    @property
    def dtypes(self) -> Frozen[Hashable, np.dtype]:
        """Mapping from coordinate names to dtypes.

        Cannot be modified directly, but is updated when adding new variables.

        See Also
        --------
        Dataset.dtypes
        """
        return Frozen(
            {
                n: v.dtype
                for n, v in self._data._variables.items()
                if n in self._data._coord_names
            }
        )

    @property
    def variables(self) -> Mapping[Hashable, Variable]:
        return Frozen(
            {k: v for k, v in self._data.variables.items() if k in self._names}
        )

    def __getitem__(self, key: Hashable) -> DataArray:
        if key in self._data.data_vars:
            raise KeyError(key)
        return self._data[key]

    def to_dataset(self) -> Dataset:
        """Convert these coordinates into a new Dataset"""

        names = [name for name in self._data._variables if name in self._names]
        return self._data._copy_listed(names)

    def _update_coords(
        self, coords: dict[Hashable, Variable], indexes: Mapping[Any, Index]
    ) -> None:
        variables = self._data._variables.copy()
        variables.update(coords)

        # check for inconsistent state *before* modifying anything in-place
        dims = calculate_dimensions(variables)
        new_coord_names = set(coords)
        for dim, size in dims.items():
            if dim in variables:
                new_coord_names.add(dim)

        self._data._variables = variables
        self._data._coord_names.update(new_coord_names)
        self._data._dims = dims

        # TODO(shoyer): once ._indexes is always populated by a dict, modify
        # it to update inplace instead.
        original_indexes = dict(self._data.xindexes)
        original_indexes.update(indexes)
        self._data._indexes = original_indexes

    def _drop_coords(self, coord_names):
        # should drop indexed coordinates only
        for name in coord_names:
            del self._data._variables[name]
            del self._data._indexes[name]
        self._data._coord_names.difference_update(coord_names)

    def _drop_indexed_coords(self, coords_to_drop: set[Hashable]) -> None:
        assert self._data.xindexes is not None
        new_coords = drop_indexed_coords(coords_to_drop, self)
        for name in self._data._coord_names - new_coords._names:
            del self._data._variables[name]
        self._data._indexes = dict(new_coords.xindexes)
        self._data._coord_names.intersection_update(new_coords._names)

    def __delitem__(self, key: Hashable) -> None:
        if key in self:
            del self._data[key]
        else:
            raise KeyError(
                f"{key!r} is not in coordinate variables {tuple(self.keys())}"
            )

    def _ipython_key_completions_(self):
        """Provide method for the key-autocompletions in IPython."""
        return [
            key
            for key in self._data._ipython_key_completions_()
            if key not in self._data.data_vars
        ]


class DataArrayCoordinates(Coordinates, Generic[T_DataArray]):
    """Dictionary like container for DataArray coordinates (variables + indexes).

    This collection can be passed directly to the :py:class:`~xarray.Dataset`
    and :py:class:`~xarray.DataArray` constructors via their `coords` argument.
    This will add both the coordinates variables and their index.
    """

    _data: T_DataArray

    __slots__ = ("_data",)

    def __init__(self, dataarray: T_DataArray) -> None:
        self._data = dataarray

    @property
    def dims(self) -> tuple[Hashable, ...]:
        return self._data.dims

    @property
    def dtypes(self) -> Frozen[Hashable, np.dtype]:
        """Mapping from coordinate names to dtypes.

        Cannot be modified directly, but is updated when adding new variables.

        See Also
        --------
        DataArray.dtype
        """
        return Frozen({n: v.dtype for n, v in self._data._coords.items()})

    @property
    def _names(self) -> set[Hashable]:
        return set(self._data._coords)

    def __getitem__(self, key: Hashable) -> T_DataArray:
        return self._data._getitem_coord(key)

    def _update_coords(
        self, coords: dict[Hashable, Variable], indexes: Mapping[Any, Index]
    ) -> None:
        coords_plus_data = coords.copy()
        coords_plus_data[_THIS_ARRAY] = self._data.variable
        dims = calculate_dimensions(coords_plus_data)
        if not set(dims) <= set(self.dims):
            raise ValueError(
                "cannot add coordinates with new dimensions to a DataArray"
            )
        self._data._coords = coords

        # TODO(shoyer): once ._indexes is always populated by a dict, modify
        # it to update inplace instead.
        original_indexes = dict(self._data.xindexes)
        original_indexes.update(indexes)
        self._data._indexes = original_indexes

    def _drop_coords(self, coord_names):
        # should drop indexed coordinates only
        for name in coord_names:
            del self._data._coords[name]
            del self._data._indexes[name]

    @property
    def variables(self):
        return Frozen(self._data._coords)

    def to_dataset(self) -> Dataset:
        from xarray.core.dataset import Dataset

        coords = {k: v.copy(deep=False) for k, v in self._data._coords.items()}
        indexes = dict(self._data.xindexes)
        return Dataset._construct_direct(coords, set(coords), indexes=indexes)

    def __delitem__(self, key: Hashable) -> None:
        if key not in self:
            raise KeyError(
                f"{key!r} is not in coordinate variables {tuple(self.keys())}"
            )
        assert_no_index_corrupted(self._data.xindexes, {key})

        del self._data._coords[key]
        if self._data._indexes is not None and key in self._data._indexes:
            del self._data._indexes[key]

    def _ipython_key_completions_(self):
        """Provide method for the key-autocompletions in IPython."""
        return self._data._ipython_key_completions_()


def drop_indexed_coords(
    coords_to_drop: set[Hashable], coords: Coordinates
) -> Coordinates:
    """Drop indexed coordinates associated with coordinates in coords_to_drop.

    This will raise an error in case it corrupts any passed index and its
    coordinate variables.

    """
    new_variables = dict(coords.variables)
    new_indexes = dict(coords.xindexes)

    for idx, idx_coords in coords.xindexes.group_by_index():
        idx_drop_coords = set(idx_coords) & coords_to_drop

        # special case for pandas multi-index: still allow but deprecate
        # dropping only its dimension coordinate.
        # TODO: remove when removing PandasMultiIndex's dimension coordinate.
        if isinstance(idx, PandasMultiIndex) and idx_drop_coords == {idx.dim}:
            idx_drop_coords.update(idx.index.names)
            emit_user_level_warning(
                f"updating coordinate {idx.dim!r} with a PandasMultiIndex would leave "
                f"the multi-index level coordinates {list(idx.index.names)!r} in an inconsistent state. "
                f"This will raise an error in the future. Use `.drop_vars({list(idx_coords)!r})` before "
                "assigning new coordinate values.",
                FutureWarning,
            )

        elif idx_drop_coords and len(idx_drop_coords) != len(idx_coords):
            idx_drop_coords_str = ", ".join(f"{k!r}" for k in idx_drop_coords)
            idx_coords_str = ", ".join(f"{k!r}" for k in idx_coords)
            raise ValueError(
                f"cannot drop or update coordinate(s) {idx_drop_coords_str}, which would corrupt "
                f"the following index built from coordinates {idx_coords_str}:\n"
                f"{idx}"
            )

        for k in idx_drop_coords:
            del new_variables[k]
            del new_indexes[k]

    return Coordinates._construct_direct(coords=new_variables, indexes=new_indexes)


def assert_coordinate_consistent(obj: T_Xarray, coords: Mapping[Any, Variable]) -> None:
    """Make sure the dimension coordinate of obj is consistent with coords.

    obj: DataArray or Dataset
    coords: Dict-like of variables
    """
    for k in obj.dims:
        # make sure there are no conflict in dimension coordinates
        if k in coords and k in obj.coords and not coords[k].equals(obj[k].variable):
            raise IndexError(
                f"dimension coordinate {k!r} conflicts between "
                f"indexed and indexing objects:\n{obj[k]}\nvs.\n{coords[k]}"
            )


def create_coords_with_default_indexes(
    coords: Mapping[Any, Any], data_vars: DataVars | None = None
) -> Coordinates:
    """Returns a Coordinates object from a mapping of coordinates (arbitrary objects).

    Create default (pandas) indexes for each of the input dimension coordinates.
    Extract coordinates from each input DataArray.

    """
    # Note: data_vars is needed here only because a pd.MultiIndex object
    # can be promoted as coordinates.
    # TODO: It won't be relevant anymore when this behavior will be dropped
    # in favor of the more explicit ``Coordinates.from_pandas_multiindex()``.

    from xarray.core.dataarray import DataArray

    all_variables = dict(coords)
    if data_vars is not None:
        all_variables.update(data_vars)

    indexes: dict[Hashable, Index] = {}
    variables: dict[Hashable, Variable] = {}

    # promote any pandas multi-index in data_vars as coordinates
    coords_promoted: dict[Hashable, Any] = {}
    pd_mindex_keys: list[Hashable] = []

    for k, v in all_variables.items():
        if isinstance(v, pd.MultiIndex):
            coords_promoted[k] = v
            pd_mindex_keys.append(k)
        elif k in coords:
            coords_promoted[k] = v

    if pd_mindex_keys:
        pd_mindex_keys_fmt = ",".join([f"'{k}'" for k in pd_mindex_keys])
        emit_user_level_warning(
            f"the `pandas.MultiIndex` object(s) passed as {pd_mindex_keys_fmt} coordinate(s) or "
            "data variable(s) will no longer be implicitly promoted and wrapped into "
            "multiple indexed coordinates in the future "
            "(i.e., one coordinate for each multi-index level + one dimension coordinate). "
            "If you want to keep this behavior, you need to first wrap it explicitly using "
            "`mindex_coords = xarray.Coordinates.from_pandas_multiindex(mindex_obj, 'dim')` "
            "and pass it as coordinates, e.g., `xarray.Dataset(coords=mindex_coords)`, "
            "`dataset.assign_coords(mindex_coords)` or `dataarray.assign_coords(mindex_coords)`.",
            FutureWarning,
        )

    dataarray_coords: list[DataArrayCoordinates] = []

    for name, obj in coords_promoted.items():
        if isinstance(obj, DataArray):
            dataarray_coords.append(obj.coords)

        variable = as_variable(obj, name=name, auto_convert=False)

        if variable.dims == (name,):
            # still needed to convert to IndexVariable first due to some
            # pandas multi-index edge cases.
            variable = variable.to_index_variable()
            idx, idx_vars = create_default_index_implicit(variable, all_variables)
            indexes.update({k: idx for k in idx_vars})
            variables.update(idx_vars)
            all_variables.update(idx_vars)
        else:
            variables[name] = variable

    new_coords = Coordinates._construct_direct(coords=variables, indexes=indexes)

    # extract and merge coordinates and indexes from input DataArrays
    if dataarray_coords:
        prioritized = {k: (v, indexes.get(k, None)) for k, v in variables.items()}
        variables, indexes = merge_coordinates_without_align(
            dataarray_coords + [new_coords],
            prioritized=prioritized,
        )
        new_coords = Coordinates._construct_direct(coords=variables, indexes=indexes)

    return new_coords
