from __future__ import annotations

import functools
import operator
from collections import defaultdict
from collections.abc import Callable, Hashable, Iterable, Mapping
from contextlib import suppress
from itertools import starmap
from typing import TYPE_CHECKING, Any, Final, Generic, TypeVar, get_args, overload

import numpy as np
import pandas as pd

from xarray.core import dtypes
from xarray.core.indexes import (
    Index,
    Indexes,
    PandasIndex,
    PandasMultiIndex,
    indexes_all_equal,
    safe_cast_to_index,
)
from xarray.core.types import JoinOptions, T_Alignable
from xarray.core.utils import emit_user_level_warning, is_dict_like, is_full_slice
from xarray.core.variable import Variable, as_compatible_data, calculate_dimensions
from xarray.util.deprecation_helpers import CombineKwargDefault

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import (
        Alignable,
        T_DataArray,
        T_Dataset,
        T_DuckArray,
    )


class AlignmentError(ValueError):
    """Error class for alignment failures due to incompatible arguments."""


def reindex_variables(
    variables: Mapping[Any, Variable],
    dim_pos_indexers: Mapping[Any, Any],
    copy: bool = True,
    fill_value: Any = dtypes.NA,
    sparse: bool = False,
) -> dict[Hashable, Variable]:
    """Conform a dictionary of variables onto a new set of variables reindexed
    with dimension positional indexers and possibly filled with missing values.

    Not public API.

    """
    new_variables = {}
    dim_sizes = calculate_dimensions(variables)

    masked_dims = set()
    unchanged_dims = set()
    for dim, indxr in dim_pos_indexers.items():
        # Negative values in dim_pos_indexers mean values missing in the new index
        # See ``Index.reindex_like``.
        if (indxr < 0).any():
            masked_dims.add(dim)
        elif np.array_equal(indxr, np.arange(dim_sizes.get(dim, 0))):
            unchanged_dims.add(dim)

    for name, var in variables.items():
        if isinstance(fill_value, dict):
            fill_value_ = fill_value.get(name, dtypes.NA)
        else:
            fill_value_ = fill_value

        if sparse:
            var = var._as_sparse(fill_value=fill_value_)
        indxr = tuple(
            slice(None) if d in unchanged_dims else dim_pos_indexers.get(d, slice(None))
            for d in var.dims
        )
        needs_masking = any(d in masked_dims for d in var.dims)

        if needs_masking:
            new_var = var._getitem_with_mask(indxr, fill_value=fill_value_)
        elif all(is_full_slice(k) for k in indxr):
            # no reindexing necessary
            # here we need to manually deal with copying data, since
            # we neither created a new ndarray nor used fancy indexing
            new_var = var.copy(deep=copy)
        else:
            new_var = var[indxr]

        new_variables[name] = new_var

    return new_variables


def _normalize_indexes(
    indexes: Mapping[Any, Any | T_DuckArray],
) -> Indexes:
    """Normalize the indexes/indexers given for re-indexing or alignment.

    Wrap any arbitrary array or `pandas.Index` as an Xarray `PandasIndex`
    associated with its corresponding dimension coordinate variable.

    """
    xr_indexes: dict[Hashable, Index] = {}
    xr_variables: dict[Hashable, Variable]

    if isinstance(indexes, Indexes):
        xr_variables = dict(indexes.variables)
    else:
        xr_variables = {}

    for k, idx in indexes.items():
        if not isinstance(idx, Index):
            if getattr(idx, "dims", (k,)) != (k,):
                raise AlignmentError(
                    f"Indexer has dimensions {idx.dims} that are different "
                    f"from that to be indexed along '{k}'"
                )
            data: T_DuckArray = as_compatible_data(idx)
            pd_idx = safe_cast_to_index(data)
            if pd_idx.name != k:
                pd_idx = pd_idx.copy()
                pd_idx.name = k
            if isinstance(pd_idx, pd.MultiIndex):
                idx = PandasMultiIndex(pd_idx, k)
            else:
                idx = PandasIndex(pd_idx, k, coord_dtype=data.dtype)
            xr_variables.update(idx.create_variables())
        xr_indexes[k] = idx

    return Indexes(xr_indexes, xr_variables)


CoordNamesAndDims = tuple[tuple[Hashable, tuple[Hashable, ...]], ...]
MatchingIndexKey = tuple[CoordNamesAndDims, type[Index]]
IndexesToAlign = dict[MatchingIndexKey, Index]
IndexVarsToAlign = dict[MatchingIndexKey, dict[Hashable, Variable]]


class Aligner(Generic[T_Alignable]):
    """Implements all the complex logic for the re-indexing and alignment of Xarray
    objects.

    For internal use only, not public API.
    Usage:

    aligner = Aligner(*objects, **kwargs)
    aligner.align()
    aligned_objects = aligner.results

    """

    objects: tuple[T_Alignable, ...]
    results: tuple[T_Alignable, ...]
    objects_matching_index_vars: tuple[
        dict[MatchingIndexKey, dict[Hashable, Variable]], ...
    ]
    join: JoinOptions | CombineKwargDefault
    exclude_dims: frozenset[Hashable]
    exclude_vars: frozenset[Hashable]
    copy: bool
    fill_value: Any
    sparse: bool
    indexes: dict[MatchingIndexKey, Index]
    index_vars: dict[MatchingIndexKey, dict[Hashable, Variable]]
    all_indexes: dict[MatchingIndexKey, list[Index]]
    all_index_vars: dict[MatchingIndexKey, list[dict[Hashable, Variable]]]
    aligned_indexes: dict[MatchingIndexKey, Index]
    aligned_index_vars: dict[MatchingIndexKey, dict[Hashable, Variable]]
    reindex: dict[MatchingIndexKey, bool]
    keep_original_indexes: set[MatchingIndexKey]
    reindex_kwargs: dict[str, Any]
    unindexed_dim_sizes: dict[Hashable, set]
    new_indexes: Indexes[Index]

    def __init__(
        self,
        objects: Iterable[T_Alignable],
        join: JoinOptions | CombineKwargDefault = "inner",
        indexes: Mapping[Any, Any] | None = None,
        exclude_dims: str | Iterable[Hashable] = frozenset(),
        exclude_vars: Iterable[Hashable] = frozenset(),
        method: str | None = None,
        tolerance: float | Iterable[float] | str | None = None,
        copy: bool = True,
        fill_value: Any = dtypes.NA,
        sparse: bool = False,
    ):
        self.objects = tuple(objects)
        self.objects_matching_indexes: tuple[Any, ...] = ()
        self.objects_matching_index_vars = ()

        if not isinstance(join, CombineKwargDefault) and join not in get_args(
            JoinOptions
        ):
            raise ValueError(f"invalid value for join: {join}")
        self.join = join

        self.copy = copy
        self.fill_value = fill_value
        self.sparse = sparse

        if method is None and tolerance is None:
            self.reindex_kwargs = {}
        else:
            self.reindex_kwargs = {"method": method, "tolerance": tolerance}

        if isinstance(exclude_dims, str):
            exclude_dims = [exclude_dims]
        self.exclude_dims = frozenset(exclude_dims)
        self.exclude_vars = frozenset(exclude_vars)

        if indexes is None:
            indexes = {}
        self.indexes, self.index_vars = self._collect_indexes(
            _normalize_indexes(indexes)
        )

        self.all_indexes = {}
        self.all_index_vars = {}
        self.unindexed_dim_sizes = {}

        self.aligned_indexes = {}
        self.aligned_index_vars = {}
        self.reindex = {}
        self.keep_original_indexes = set()

        self.results = tuple()

    def _collect_indexes(
        self, indexes: Indexes
    ) -> tuple[IndexesToAlign, IndexVarsToAlign]:
        """Collect input and/or object indexes for alignment.

        Return new dictionaries of xarray Index objects and coordinate
        variables, whose keys are used to later retrieve all the indexes to
        compare with each other (based on the name and dimensions of their
        associated coordinate variables as well as the Index type).

        """
        collected_indexes = {}
        collected_index_vars = {}

        for idx, idx_vars in indexes.group_by_index():
            idx_coord_names_and_dims = []
            idx_all_dims: set[Hashable] = set()

            for name, var in idx_vars.items():
                dims = var.dims
                idx_coord_names_and_dims.append((name, dims))
                idx_all_dims.update(dims)

            key: MatchingIndexKey = (tuple(idx_coord_names_and_dims), type(idx))

            if idx_all_dims:
                exclude_dims = idx_all_dims & self.exclude_dims
                if exclude_dims == idx_all_dims:
                    # Do not collect an index if all the dimensions it uses are
                    # also excluded from the alignment
                    continue
                elif exclude_dims:
                    # If the dimensions used by index partially overlap with the dimensions
                    # excluded from alignment, it is possible to check index equality along
                    # non-excluded dimensions only. However, in this case each of the aligned
                    # objects must retain (a copy of) their original index. Re-indexing and
                    # overriding the index are not supported.
                    if self.join == "override":
                        excl_dims_str = ", ".join(str(d) for d in exclude_dims)
                        incl_dims_str = ", ".join(
                            str(d) for d in idx_all_dims - exclude_dims
                        )
                        raise AlignmentError(
                            f"cannot exclude dimension(s) {excl_dims_str} from alignment "
                            "with `join='override` because these are used by an index "
                            f"together with non-excluded dimensions {incl_dims_str}"
                            "(cannot safely override the index)."
                        )
                    else:
                        self.keep_original_indexes.add(key)

            collected_indexes[key] = idx
            collected_index_vars[key] = idx_vars

        return collected_indexes, collected_index_vars

    def find_matching_indexes(self) -> None:
        all_indexes: dict[MatchingIndexKey, list[Index]]
        all_index_vars: dict[MatchingIndexKey, list[dict[Hashable, Variable]]]
        all_indexes_dim_sizes: dict[MatchingIndexKey, dict[Hashable, set]]
        objects_matching_indexes: list[dict[MatchingIndexKey, Index]]
        objects_matching_index_vars: list[
            dict[MatchingIndexKey, dict[Hashable, Variable]]
        ]

        all_indexes = defaultdict(list)
        all_index_vars = defaultdict(list)
        all_indexes_dim_sizes = defaultdict(lambda: defaultdict(set))
        objects_matching_indexes = []
        objects_matching_index_vars = []

        for obj in self.objects:
            obj_indexes, obj_index_vars = self._collect_indexes(obj.xindexes)
            objects_matching_indexes.append(obj_indexes)
            objects_matching_index_vars.append(obj_index_vars)
            for key, idx in obj_indexes.items():
                all_indexes[key].append(idx)
            for key, index_vars in obj_index_vars.items():
                all_index_vars[key].append(index_vars)
                for dim, size in calculate_dimensions(index_vars).items():
                    all_indexes_dim_sizes[key][dim].add(size)

        self.objects_matching_indexes = tuple(objects_matching_indexes)
        self.objects_matching_index_vars = tuple(objects_matching_index_vars)
        self.all_indexes = all_indexes
        self.all_index_vars = all_index_vars

        if self.join == "override":
            for dim_sizes in all_indexes_dim_sizes.values():
                for dim, sizes in dim_sizes.items():
                    if len(sizes) > 1:
                        raise AlignmentError(
                            "cannot align objects with join='override' with matching indexes "
                            f"along dimension {dim!r} that don't have the same size"
                        )

    def find_matching_unindexed_dims(self) -> None:
        unindexed_dim_sizes = defaultdict(set)

        for obj in self.objects:
            for dim in obj.dims:
                if dim not in self.exclude_dims and dim not in obj.xindexes.dims:
                    unindexed_dim_sizes[dim].add(obj.sizes[dim])

        self.unindexed_dim_sizes = unindexed_dim_sizes

    def _need_reindex(self, dim, cmp_indexes) -> bool:
        """Whether or not we need to reindex variables for a set of
        matching indexes.

        We don't reindex when all matching indexes are equal for two reasons:
        - It's faster for the usual case (already aligned objects).
        - It ensures it's possible to do operations that don't require alignment
          on indexes with duplicate values (which cannot be reindexed with
          pandas). This is useful, e.g., for overwriting such duplicate indexes.

        """
        if not indexes_all_equal(cmp_indexes, self.exclude_dims):
            # always reindex when matching indexes are not equal
            return True

        unindexed_dims_sizes = {}
        for d in dim:
            if d in self.unindexed_dim_sizes:
                sizes = self.unindexed_dim_sizes[d]
                if len(sizes) > 1:
                    # reindex if different sizes are found for unindexed dims
                    return True
                else:
                    unindexed_dims_sizes[d] = next(iter(sizes))

        if unindexed_dims_sizes:
            indexed_dims_sizes = {}
            for cmp in cmp_indexes:
                index_vars = cmp[1]
                for var in index_vars.values():
                    indexed_dims_sizes.update(var.sizes)

            for d, size in unindexed_dims_sizes.items():
                if indexed_dims_sizes.get(d, -1) != size:
                    # reindex if unindexed dimension size doesn't match
                    return True

        return False

    def _get_index_joiner(self, index_cls) -> Callable:
        if self.join in ["outer", "inner"]:
            return functools.partial(
                functools.reduce,
                functools.partial(index_cls.join, how=self.join),
            )
        elif self.join == "left":
            return operator.itemgetter(0)
        elif self.join == "right":
            return operator.itemgetter(-1)
        elif self.join == "override":
            # We rewrite all indexes and then use join='left'
            return operator.itemgetter(0)
        else:
            # join='exact' return dummy lambda (error is raised)
            return lambda _: None

    def align_indexes(self) -> None:
        """Compute all aligned indexes and their corresponding coordinate variables."""

        aligned_indexes: dict[MatchingIndexKey, Index] = {}
        aligned_index_vars: dict[MatchingIndexKey, dict[Hashable, Variable]] = {}
        reindex: dict[MatchingIndexKey, bool] = {}
        new_indexes: dict[Hashable, Index] = {}
        new_index_vars: dict[Hashable, Variable] = {}

        def update_dicts(
            key: MatchingIndexKey,
            idx: Index,
            idx_vars: dict[Hashable, Variable],
            need_reindex: bool,
        ):
            reindex[key] = need_reindex
            aligned_indexes[key] = idx
            aligned_index_vars[key] = idx_vars

            for name, var in idx_vars.items():
                if name in new_indexes:
                    other_idx = new_indexes[name]
                    other_var = new_index_vars[name]
                    raise AlignmentError(
                        f"cannot align objects on coordinate {name!r} because of conflicting indexes\n"
                        f"first index: {idx!r}\nsecond index: {other_idx!r}\n"
                        f"first variable: {var!r}\nsecond variable: {other_var!r}\n"
                    )
                new_indexes[name] = idx
                new_index_vars[name] = var

        for key, matching_indexes in self.all_indexes.items():
            matching_index_vars = self.all_index_vars[key]
            dims = {d for coord in matching_index_vars[0].values() for d in coord.dims}
            index_cls = key[1]

            if self.join == "override":
                joined_index = matching_indexes[0]
                joined_index_vars = matching_index_vars[0]
                need_reindex = False
            elif key in self.indexes:
                joined_index = self.indexes[key]
                joined_index_vars = self.index_vars[key]
                cmp_indexes = list(
                    zip(
                        [joined_index] + matching_indexes,
                        [joined_index_vars] + matching_index_vars,
                        strict=True,
                    )
                )
                need_reindex = self._need_reindex(dims, cmp_indexes)
            else:
                if len(matching_indexes) > 1:
                    need_reindex = self._need_reindex(
                        dims,
                        list(zip(matching_indexes, matching_index_vars, strict=True)),
                    )
                else:
                    need_reindex = False
                if need_reindex:
                    if (
                        isinstance(self.join, CombineKwargDefault)
                        and self.join != "exact"
                    ):
                        emit_user_level_warning(
                            self.join.warning_message(
                                "This change will result in the following ValueError: "
                                "cannot be aligned with join='exact' because "
                                "index/labels/sizes are not equal along "
                                "these coordinates (dimensions): "
                                + ", ".join(
                                    f"{name!r} {dims!r}" for name, dims in key[0]
                                ),
                                recommend_set_options=False,
                            ),
                            FutureWarning,
                        )
                    if self.join == "exact":
                        raise AlignmentError(
                            "cannot align objects with join='exact' where "
                            "index/labels/sizes are not equal along "
                            "these coordinates (dimensions): "
                            + ", ".join(f"{name!r} {dims!r}" for name, dims in key[0])
                            + (
                                self.join.error_message()
                                if isinstance(self.join, CombineKwargDefault)
                                else ""
                            )
                        )
                    joiner = self._get_index_joiner(index_cls)
                    joined_index = joiner(matching_indexes)
                    if self.join == "left":
                        joined_index_vars = matching_index_vars[0]
                    elif self.join == "right":
                        joined_index_vars = matching_index_vars[-1]
                    else:
                        joined_index_vars = joined_index.create_variables()
                else:
                    joined_index = matching_indexes[0]
                    joined_index_vars = matching_index_vars[0]

            update_dicts(key, joined_index, joined_index_vars, need_reindex)

        # Explicitly provided indexes that are not found in objects to align
        # may relate to unindexed dimensions so we add them too
        for key, idx in self.indexes.items():
            if key not in aligned_indexes:
                index_vars = self.index_vars[key]
                update_dicts(key, idx, index_vars, False)

        self.aligned_indexes = aligned_indexes
        self.aligned_index_vars = aligned_index_vars
        self.reindex = reindex
        self.new_indexes = Indexes(new_indexes, new_index_vars)

    def assert_unindexed_dim_sizes_equal(self) -> None:
        for dim, sizes in self.unindexed_dim_sizes.items():
            index_size = self.new_indexes.dims.get(dim)
            if index_size is not None:
                sizes.add(index_size)
                add_err_msg = (
                    f" (note: an index is found along that dimension "
                    f"with size={index_size!r})"
                )
            else:
                add_err_msg = ""
            if len(sizes) > 1:
                raise AlignmentError(
                    f"cannot reindex or align along dimension {dim!r} "
                    f"because of conflicting dimension sizes: {sizes!r}" + add_err_msg
                )

    def override_indexes(self) -> None:
        objects = list(self.objects)

        for i, obj in enumerate(objects[1:]):
            new_indexes = {}
            new_variables = {}
            matching_indexes = self.objects_matching_indexes[i + 1]

            for key, aligned_idx in self.aligned_indexes.items():
                obj_idx = matching_indexes.get(key)
                if obj_idx is not None:
                    for name, var in self.aligned_index_vars[key].items():
                        new_indexes[name] = aligned_idx
                        new_variables[name] = var.copy(deep=self.copy)

            objects[i + 1] = obj._overwrite_indexes(new_indexes, new_variables)

        self.results = tuple(objects)

    def _get_dim_pos_indexers(
        self,
        matching_indexes: dict[MatchingIndexKey, Index],
    ) -> dict[Hashable, Any]:
        dim_pos_indexers: dict[Hashable, Any] = {}
        dim_index: dict[Hashable, Index] = {}

        for key, aligned_idx in self.aligned_indexes.items():
            obj_idx = matching_indexes.get(key)
            if obj_idx is not None and self.reindex[key]:
                indexers = obj_idx.reindex_like(aligned_idx, **self.reindex_kwargs)
                for dim, idxer in indexers.items():
                    if dim in self.exclude_dims:
                        raise AlignmentError(
                            f"cannot reindex or align along dimension {dim!r} because "
                            "it is explicitly excluded from alignment. This is likely caused by "
                            "wrong results returned by the `reindex_like` method of this index:\n"
                            f"{obj_idx!r}"
                        )
                    if dim in dim_pos_indexers and not np.array_equal(
                        idxer, dim_pos_indexers[dim]
                    ):
                        raise AlignmentError(
                            f"cannot reindex or align along dimension {dim!r} because "
                            "of conflicting re-indexers returned by multiple indexes\n"
                            f"first index: {obj_idx!r}\nsecond index: {dim_index[dim]!r}\n"
                        )
                    dim_pos_indexers[dim] = idxer
                    dim_index[dim] = obj_idx

        return dim_pos_indexers

    def _get_indexes_and_vars(
        self,
        obj: T_Alignable,
        matching_indexes: dict[MatchingIndexKey, Index],
        matching_index_vars: dict[MatchingIndexKey, dict[Hashable, Variable]],
    ) -> tuple[dict[Hashable, Index], dict[Hashable, Variable]]:
        new_indexes = {}
        new_variables = {}

        for key, aligned_idx in self.aligned_indexes.items():
            aligned_idx_vars = self.aligned_index_vars[key]
            obj_idx = matching_indexes.get(key)
            obj_idx_vars = matching_index_vars.get(key)

            if obj_idx is None:
                # add the aligned index if it relates to unindexed dimensions in obj
                dims = {d for var in aligned_idx_vars.values() for d in var.dims}
                if dims <= set(obj.dims):
                    obj_idx = aligned_idx

            if obj_idx is not None:
                # TODO: always copy object's index when no re-indexing is required?
                # (instead of assigning the aligned index)
                # (need performance assessment)
                if key in self.keep_original_indexes:
                    assert self.reindex[key] is False
                    new_idx = obj_idx.copy(deep=self.copy)
                    new_idx_vars = new_idx.create_variables(obj_idx_vars)
                else:
                    new_idx = aligned_idx
                    new_idx_vars = {
                        k: v.copy(deep=self.copy) for k, v in aligned_idx_vars.items()
                    }
                new_indexes.update(dict.fromkeys(new_idx_vars, new_idx))
                new_variables.update(new_idx_vars)

        return new_indexes, new_variables

    def _reindex_one(
        self,
        obj: T_Alignable,
        matching_indexes: dict[MatchingIndexKey, Index],
        matching_index_vars: dict[MatchingIndexKey, dict[Hashable, Variable]],
    ) -> T_Alignable:
        new_indexes, new_variables = self._get_indexes_and_vars(
            obj, matching_indexes, matching_index_vars
        )
        dim_pos_indexers = self._get_dim_pos_indexers(matching_indexes)

        return obj._reindex_callback(
            self,
            dim_pos_indexers,
            new_variables,
            new_indexes,
            self.fill_value,
            self.exclude_dims,
            self.exclude_vars,
        )

    def reindex_all(self) -> None:
        self.results = tuple(
            starmap(
                self._reindex_one,
                zip(
                    self.objects,
                    self.objects_matching_indexes,
                    self.objects_matching_index_vars,
                    strict=True,
                ),
            )
        )

    def align(self) -> None:
        if not self.indexes and len(self.objects) == 1:
            # fast path for the trivial case
            (obj,) = self.objects
            self.results = (obj.copy(deep=self.copy),)
            return

        self.find_matching_indexes()
        self.find_matching_unindexed_dims()
        self.align_indexes()
        self.assert_unindexed_dim_sizes_equal()

        if self.join == "override":
            self.override_indexes()
        elif self.join == "exact" and not self.copy:
            self.results = self.objects
        else:
            self.reindex_all()


T_Obj1 = TypeVar("T_Obj1", bound="Alignable")
T_Obj2 = TypeVar("T_Obj2", bound="Alignable")
T_Obj3 = TypeVar("T_Obj3", bound="Alignable")
T_Obj4 = TypeVar("T_Obj4", bound="Alignable")
T_Obj5 = TypeVar("T_Obj5", bound="Alignable")


@overload
def align(
    obj1: T_Obj1,
    /,
    *,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Obj1]: ...


@overload
def align(
    obj1: T_Obj1,
    obj2: T_Obj2,
    /,
    *,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Obj1, T_Obj2]: ...


@overload
def align(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    /,
    *,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Obj1, T_Obj2, T_Obj3]: ...


@overload
def align(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    obj4: T_Obj4,
    /,
    *,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Obj1, T_Obj2, T_Obj3, T_Obj4]: ...


@overload
def align(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    obj4: T_Obj4,
    obj5: T_Obj5,
    /,
    *,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Obj1, T_Obj2, T_Obj3, T_Obj4, T_Obj5]: ...


@overload
def align(
    *objects: T_Alignable,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Alignable, ...]: ...


def align(
    *objects: T_Alignable,
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    fill_value=dtypes.NA,
) -> tuple[T_Alignable, ...]:
    """
    Given any number of Dataset and/or DataArray objects, returns new
    objects with aligned indexes and dimension sizes.

    Array from the aligned objects are suitable as input to mathematical
    operators, because along each dimension they have the same index and size.

    Missing values (if ``join != 'inner'``) are filled with ``fill_value``.
    The default fill value is NaN.

    Parameters
    ----------
    *objects : Dataset or DataArray
        Objects to align.
    join : {"outer", "inner", "left", "right", "exact", "override"}, optional
        Method for joining the indexes of the passed objects along each
        dimension:

        - "outer": use the union of object indexes
        - "inner": use the intersection of object indexes
        - "left": use indexes from the first object with each dimension
        - "right": use indexes from the last object with each dimension
        - "exact": instead of aligning, raise `ValueError` when indexes to be
          aligned are not equal
        - "override": if indexes are of same size, rewrite indexes to be
          those of the first object with that dimension. Indexes for the same
          dimension must have the same size in all objects.

    copy : bool, default: True
        If ``copy=True``, data in the return values is always copied. If
        ``copy=False`` and reindexing is unnecessary, or can be performed with
        only slice operations, then the output may share memory with the input.
        In either case, new xarray objects are always returned.
    indexes : dict-like, optional
        Any indexes explicitly provided with the `indexes` argument should be
        used in preference to the aligned indexes.
    exclude : str, iterable of hashable or None, optional
        Dimensions that must be excluded from alignment
    fill_value : scalar or dict-like, optional
        Value to use for newly missing values. If a dict-like, maps
        variable names to fill values. Use a data array's name to
        refer to its values.

    Returns
    -------
    aligned : tuple of DataArray or Dataset
        Tuple of objects with the same type as `*objects` with aligned
        coordinates.

    Raises
    ------
    AlignmentError
        If any dimensions without labels on the arguments have different sizes,
        or a different size than the size of the aligned dimension labels.

    Examples
    --------
    >>> x = xr.DataArray(
    ...     [[25, 35], [10, 24]],
    ...     dims=("lat", "lon"),
    ...     coords={"lat": [35.0, 40.0], "lon": [100.0, 120.0]},
    ... )
    >>> y = xr.DataArray(
    ...     [[20, 5], [7, 13]],
    ...     dims=("lat", "lon"),
    ...     coords={"lat": [35.0, 42.0], "lon": [100.0, 120.0]},
    ... )

    >>> x
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[25, 35],
           [10, 24]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> y
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[20,  5],
           [ 7, 13]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 42.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y)
    >>> a
    <xarray.DataArray (lat: 1, lon: 2)> Size: 16B
    array([[25, 35]])
    Coordinates:
      * lat      (lat) float64 8B 35.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 1, lon: 2)> Size: 16B
    array([[20,  5]])
    Coordinates:
      * lat      (lat) float64 8B 35.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y, join="outer")
    >>> a
    <xarray.DataArray (lat: 3, lon: 2)> Size: 48B
    array([[25., 35.],
           [10., 24.],
           [nan, nan]])
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 3, lon: 2)> Size: 48B
    array([[20.,  5.],
           [nan, nan],
           [ 7., 13.]])
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y, join="outer", fill_value=-999)
    >>> a
    <xarray.DataArray (lat: 3, lon: 2)> Size: 48B
    array([[  25,   35],
           [  10,   24],
           [-999, -999]])
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 3, lon: 2)> Size: 48B
    array([[  20,    5],
           [-999, -999],
           [   7,   13]])
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y, join="left")
    >>> a
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[25, 35],
           [10, 24]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[20.,  5.],
           [nan, nan]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y, join="right")
    >>> a
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[25., 35.],
           [nan, nan]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 42.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[20,  5],
           [ 7, 13]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 42.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> a, b = xr.align(x, y, join="exact")
    Traceback (most recent call last):
    ...
    xarray.structure.alignment.AlignmentError: cannot align objects with join='exact' ...

    >>> a, b = xr.align(x, y, join="override")
    >>> a
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[25, 35],
           [10, 24]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0
    >>> b
    <xarray.DataArray (lat: 2, lon: 2)> Size: 32B
    array([[20,  5],
           [ 7, 13]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0

    """
    aligner = Aligner(
        objects,
        join=join,
        copy=copy,
        indexes=indexes,
        exclude_dims=exclude,
        fill_value=fill_value,
    )
    aligner.align()
    return aligner.results


def deep_align(
    objects: Iterable[Any],
    join: JoinOptions | CombineKwargDefault = "inner",
    copy: bool = True,
    indexes=None,
    exclude: str | Iterable[Hashable] = frozenset(),
    raise_on_invalid: bool = True,
    fill_value=dtypes.NA,
) -> list[Any]:
    """Align objects for merging, recursing into dictionary values.

    This function is not public API.
    """
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    if indexes is None:
        indexes = {}

    def is_alignable(obj):
        return isinstance(obj, Coordinates | DataArray | Dataset)

    positions: list[int] = []
    keys: list[type[object] | Hashable] = []
    out: list[Any] = []
    targets: list[Alignable] = []
    no_key: Final = object()
    not_replaced: Final = object()
    for position, variables in enumerate(objects):
        if is_alignable(variables):
            positions.append(position)
            keys.append(no_key)
            targets.append(variables)
            out.append(not_replaced)
        elif is_dict_like(variables):
            current_out = {}
            for k, v in variables.items():
                if is_alignable(v) and k not in indexes:
                    # Skip variables in indexes for alignment, because these
                    # should to be overwritten instead:
                    # https://github.com/pydata/xarray/issues/725
                    # https://github.com/pydata/xarray/issues/3377
                    # TODO(shoyer): doing this here feels super-hacky -- can we
                    # move it explicitly into merge instead?
                    positions.append(position)
                    keys.append(k)
                    targets.append(v)
                    current_out[k] = not_replaced
                else:
                    current_out[k] = v
            out.append(current_out)
        elif raise_on_invalid:
            raise ValueError(
                "object to align is neither an xarray.Dataset, "
                f"an xarray.DataArray nor a dictionary: {variables!r}"
            )
        else:
            out.append(variables)

    aligned = align(
        *targets,
        join=join,
        copy=copy,
        indexes=indexes,
        exclude=exclude,
        fill_value=fill_value,
    )

    for position, key, aligned_obj in zip(positions, keys, aligned, strict=True):
        if key is no_key:
            out[position] = aligned_obj
        else:
            out[position][key] = aligned_obj

    return out


def reindex(
    obj: T_Alignable,
    indexers: Mapping[Any, Any],
    method: str | None = None,
    tolerance: float | Iterable[float] | str | None = None,
    copy: bool = True,
    fill_value: Any = dtypes.NA,
    sparse: bool = False,
    exclude_vars: Iterable[Hashable] = frozenset(),
) -> T_Alignable:
    """Re-index either a Dataset or a DataArray.

    Not public API.

    """

    # TODO: (benbovy - explicit indexes): uncomment?
    # --> from reindex docstrings: "any mismatched dimension is simply ignored"
    # bad_keys = [k for k in indexers if k not in obj._indexes and k not in obj.dims]
    # if bad_keys:
    #     raise ValueError(
    #         f"indexer keys {bad_keys} do not correspond to any indexed coordinate "
    #         "or unindexed dimension in the object to reindex"
    #     )

    aligner = Aligner(
        (obj,),
        indexes=indexers,
        method=method,
        tolerance=tolerance,
        copy=copy,
        fill_value=fill_value,
        sparse=sparse,
        exclude_vars=exclude_vars,
    )
    aligner.align()
    return aligner.results[0]


def reindex_like(
    obj: T_Alignable,
    other: Dataset | DataArray,
    method: str | None = None,
    tolerance: float | Iterable[float] | str | None = None,
    copy: bool = True,
    fill_value: Any = dtypes.NA,
) -> T_Alignable:
    """Re-index either a Dataset or a DataArray like another Dataset/DataArray.

    Not public API.

    """
    if not other._indexes:
        # This check is not performed in Aligner.
        for dim in other.dims:
            if dim in obj.dims:
                other_size = other.sizes[dim]
                obj_size = obj.sizes[dim]
                if other_size != obj_size:
                    raise ValueError(
                        "different size for unlabeled "
                        f"dimension on argument {dim!r}: {other_size!r} vs {obj_size!r}"
                    )

    return reindex(
        obj,
        indexers=other.xindexes,
        method=method,
        tolerance=tolerance,
        copy=copy,
        fill_value=fill_value,
    )


def _get_broadcast_dims_map_common_coords(args, exclude):
    common_coords = {}
    dims_map = {}
    for arg in args:
        for dim in arg.dims:
            if dim not in common_coords and dim not in exclude:
                dims_map[dim] = arg.sizes[dim]
                if dim in arg._indexes:
                    common_coords.update(arg.xindexes.get_all_coords(dim))

    return dims_map, common_coords


def _broadcast_helper(
    arg: T_Alignable, exclude, dims_map, common_coords
) -> T_Alignable:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    def _set_dims(var):
        # Add excluded dims to a copy of dims_map
        var_dims_map = dims_map.copy()
        for dim in exclude:
            with suppress(ValueError):
                # ignore dim not in var.dims
                var_dims_map[dim] = var.shape[var.dims.index(dim)]

        return var.set_dims(var_dims_map)

    def _broadcast_array(array: T_DataArray) -> T_DataArray:
        data = _set_dims(array.variable)
        coords = dict(array.coords)
        coords.update(common_coords)
        return array.__class__(
            data, coords, data.dims, name=array.name, attrs=array.attrs
        )

    def _broadcast_dataset(ds: T_Dataset) -> T_Dataset:
        data_vars = {k: _set_dims(ds.variables[k]) for k in ds.data_vars}
        coords = dict(ds.coords)
        coords.update(common_coords)
        return ds.__class__(data_vars, coords, ds.attrs)

    # remove casts once https://github.com/python/mypy/issues/12800 is resolved
    if isinstance(arg, DataArray):
        return _broadcast_array(arg)  # type: ignore[return-value,unused-ignore]
    elif isinstance(arg, Dataset):
        return _broadcast_dataset(arg)  # type: ignore[return-value,unused-ignore]
    else:
        raise ValueError("all input must be Dataset or DataArray objects")


@overload
def broadcast(
    obj1: T_Obj1, /, *, exclude: str | Iterable[Hashable] | None = None
) -> tuple[T_Obj1]: ...


@overload
def broadcast(
    obj1: T_Obj1, obj2: T_Obj2, /, *, exclude: str | Iterable[Hashable] | None = None
) -> tuple[T_Obj1, T_Obj2]: ...


@overload
def broadcast(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    /,
    *,
    exclude: str | Iterable[Hashable] | None = None,
) -> tuple[T_Obj1, T_Obj2, T_Obj3]: ...


@overload
def broadcast(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    obj4: T_Obj4,
    /,
    *,
    exclude: str | Iterable[Hashable] | None = None,
) -> tuple[T_Obj1, T_Obj2, T_Obj3, T_Obj4]: ...


@overload
def broadcast(
    obj1: T_Obj1,
    obj2: T_Obj2,
    obj3: T_Obj3,
    obj4: T_Obj4,
    obj5: T_Obj5,
    /,
    *,
    exclude: str | Iterable[Hashable] | None = None,
) -> tuple[T_Obj1, T_Obj2, T_Obj3, T_Obj4, T_Obj5]: ...


@overload
def broadcast(
    *args: T_Alignable, exclude: str | Iterable[Hashable] | None = None
) -> tuple[T_Alignable, ...]: ...


def broadcast(
    *args: T_Alignable, exclude: str | Iterable[Hashable] | None = None
) -> tuple[T_Alignable, ...]:
    """Explicitly broadcast any number of DataArray or Dataset objects against
    one another.

    xarray objects automatically broadcast against each other in arithmetic
    operations, so this function should not be necessary for normal use.

    If no change is needed, the input data is returned to the output without
    being copied.

    Parameters
    ----------
    *args : DataArray or Dataset
        Arrays to broadcast against each other.
    exclude : str, iterable of hashable or None, optional
        Dimensions that must not be broadcasted

    Returns
    -------
    broadcast : tuple of DataArray or tuple of Dataset
        The same data as the input arrays, but with additional dimensions
        inserted so that all data arrays have the same dimensions and shape.

    Examples
    --------
    Broadcast two data arrays against one another to fill out their dimensions:

    >>> a = xr.DataArray([1, 2, 3], dims="x")
    >>> b = xr.DataArray([5, 6], dims="y")
    >>> a
    <xarray.DataArray (x: 3)> Size: 24B
    array([1, 2, 3])
    Dimensions without coordinates: x
    >>> b
    <xarray.DataArray (y: 2)> Size: 16B
    array([5, 6])
    Dimensions without coordinates: y
    >>> a2, b2 = xr.broadcast(a, b)
    >>> a2
    <xarray.DataArray (x: 3, y: 2)> Size: 48B
    array([[1, 1],
           [2, 2],
           [3, 3]])
    Dimensions without coordinates: x, y
    >>> b2
    <xarray.DataArray (x: 3, y: 2)> Size: 48B
    array([[5, 6],
           [5, 6],
           [5, 6]])
    Dimensions without coordinates: x, y

    Fill out the dimensions of all data variables in a dataset:

    >>> ds = xr.Dataset({"a": a, "b": b})
    >>> (ds2,) = xr.broadcast(ds)  # use tuple unpacking to extract one dataset
    >>> ds2
    <xarray.Dataset> Size: 96B
    Dimensions:  (x: 3, y: 2)
    Dimensions without coordinates: x, y
    Data variables:
        a        (x, y) int64 48B 1 1 2 2 3 3
        b        (x, y) int64 48B 5 6 5 6 5 6
    """

    if exclude is None:
        exclude = set()
    args = align(*args, join="outer", copy=False, exclude=exclude)

    dims_map, common_coords = _get_broadcast_dims_map_common_coords(args, exclude)
    result = [_broadcast_helper(arg, exclude, dims_map, common_coords) for arg in args]

    return tuple(result)
