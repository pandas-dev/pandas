from __future__ import annotations

import copy
import itertools
from typing import (
    TYPE_CHECKING,
    Dict,
    List,
    Sequence,
)

import numpy as np

from pandas._libs import internals as libinternals
from pandas._typing import (
    ArrayLike,
    DtypeObj,
    Manager,
    Shape,
)
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.cast import (
    ensure_dtype_can_hold_na,
    find_common_type,
)
from pandas.core.dtypes.common import (
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_sparse,
)
from pandas.core.dtypes.concat import concat_compat
from pandas.core.dtypes.missing import (
    is_valid_na_for_dtype,
    isna_all,
)

import pandas.core.algorithms as algos
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
)
from pandas.core.internals.array_manager import ArrayManager
from pandas.core.internals.blocks import (
    ensure_block_shape,
    new_block,
)
from pandas.core.internals.managers import BlockManager

if TYPE_CHECKING:
    from pandas import Index


def concatenate_array_managers(
    mgrs_indexers, axes: List[Index], concat_axis: int, copy: bool
) -> Manager:
    """
    Concatenate array managers into one.

    Parameters
    ----------
    mgrs_indexers : list of (ArrayManager, {axis: indexer,...}) tuples
    axes : list of Index
    concat_axis : int
    copy : bool

    Returns
    -------
    ArrayManager
    """
    # reindex all arrays
    mgrs = []
    for mgr, indexers in mgrs_indexers:
        for ax, indexer in indexers.items():
            mgr = mgr.reindex_indexer(axes[ax], indexer, axis=ax, allow_dups=True)
        mgrs.append(mgr)

    if concat_axis == 1:
        # concatting along the rows -> concat the reindexed arrays
        # TODO(ArrayManager) doesn't yet preserve the correct dtype
        arrays = [
            concat_compat([mgrs[i].arrays[j] for i in range(len(mgrs))])
            for j in range(len(mgrs[0].arrays))
        ]
        return ArrayManager(arrays, [axes[1], axes[0]], verify_integrity=False)
    else:
        # concatting along the columns -> combine reindexed arrays in a single manager
        assert concat_axis == 0
        arrays = list(itertools.chain.from_iterable([mgr.arrays for mgr in mgrs]))
        return ArrayManager(arrays, [axes[1], axes[0]], verify_integrity=False)


def concatenate_managers(
    mgrs_indexers, axes: List[Index], concat_axis: int, copy: bool
) -> Manager:
    """
    Concatenate block managers into one.

    Parameters
    ----------
    mgrs_indexers : list of (BlockManager, {axis: indexer,...}) tuples
    axes : list of Index
    concat_axis : int
    copy : bool

    Returns
    -------
    BlockManager
    """
    # TODO(ArrayManager) this assumes that all managers are of the same type
    if isinstance(mgrs_indexers[0][0], ArrayManager):
        return concatenate_array_managers(mgrs_indexers, axes, concat_axis, copy)

    concat_plans = [
        _get_mgr_concatenation_plan(mgr, indexers) for mgr, indexers in mgrs_indexers
    ]
    concat_plan = _combine_concat_plans(concat_plans, concat_axis)
    blocks = []

    for placement, join_units in concat_plan:

        if len(join_units) == 1 and not join_units[0].indexers:
            b = join_units[0].block
            values = b.values
            if copy:
                values = values.copy()
            else:
                values = values.view()
            b = b.make_block_same_class(values, placement=placement)
        elif _is_uniform_join_units(join_units):
            blk = join_units[0].block
            vals = [ju.block.values for ju in join_units]

            if not blk.is_extension:
                # _is_uniform_join_units ensures a single dtype, so
                #  we can use np.concatenate, which is more performant
                #  than concat_compat
                values = np.concatenate(vals, axis=blk.ndim - 1)
            else:
                # TODO(EA2D): special-casing not needed with 2D EAs
                values = concat_compat(vals)
                if not isinstance(values, ExtensionArray):
                    values = values.reshape(1, len(values))

            if blk.values.dtype == values.dtype:
                # Fast-path
                b = blk.make_block_same_class(values, placement=placement)
            else:
                b = new_block(values, placement=placement, ndim=blk.ndim)
        else:
            new_values = _concatenate_join_units(join_units, concat_axis, copy=copy)
            b = new_block(new_values, placement=placement, ndim=len(axes))
        blocks.append(b)

    return BlockManager(blocks, axes)


def _get_mgr_concatenation_plan(mgr: BlockManager, indexers: Dict[int, np.ndarray]):
    """
    Construct concatenation plan for given block manager and indexers.

    Parameters
    ----------
    mgr : BlockManager
    indexers : dict of {axis: indexer}

    Returns
    -------
    plan : list of (BlockPlacement, JoinUnit) tuples

    """
    # Calculate post-reindex shape , save for item axis which will be separate
    # for each block anyway.
    mgr_shape_list = list(mgr.shape)
    for ax, indexer in indexers.items():
        mgr_shape_list[ax] = len(indexer)
    mgr_shape = tuple(mgr_shape_list)

    if 0 in indexers:
        ax0_indexer = indexers.pop(0)
        blknos = algos.take_nd(mgr.blknos, ax0_indexer, fill_value=-1)
        blklocs = algos.take_nd(mgr.blklocs, ax0_indexer, fill_value=-1)
    else:

        if mgr.is_single_block:
            blk = mgr.blocks[0]
            return [(blk.mgr_locs, JoinUnit(blk, mgr_shape, indexers))]

        # error: Incompatible types in assignment (expression has type "None", variable
        # has type "ndarray")
        ax0_indexer = None  # type: ignore[assignment]
        blknos = mgr.blknos
        blklocs = mgr.blklocs

    plan = []
    for blkno, placements in libinternals.get_blkno_placements(blknos, group=False):

        assert placements.is_slice_like

        join_unit_indexers = indexers.copy()

        shape_list = list(mgr_shape)
        shape_list[0] = len(placements)
        shape = tuple(shape_list)

        if blkno == -1:
            unit = JoinUnit(None, shape)
        else:
            blk = mgr.blocks[blkno]
            ax0_blk_indexer = blklocs[placements.indexer]

            unit_no_ax0_reindexing = (
                len(placements) == len(blk.mgr_locs)
                and
                # Fastpath detection of join unit not
                # needing to reindex its block: no ax0
                # reindexing took place and block
                # placement was sequential before.
                (
                    (
                        ax0_indexer is None
                        and blk.mgr_locs.is_slice_like
                        and blk.mgr_locs.as_slice.step == 1
                    )
                    or
                    # Slow-ish detection: all indexer locs
                    # are sequential (and length match is
                    # checked above).
                    (np.diff(ax0_blk_indexer) == 1).all()
                )
            )

            # Omit indexer if no item reindexing is required.
            if unit_no_ax0_reindexing:
                join_unit_indexers.pop(0, None)
            else:
                join_unit_indexers[0] = ax0_blk_indexer

            unit = JoinUnit(blk, shape, join_unit_indexers)

        plan.append((placements, unit))

    return plan


class JoinUnit:
    def __init__(self, block, shape: Shape, indexers=None):
        # Passing shape explicitly is required for cases when block is None.
        if indexers is None:
            indexers = {}
        self.block = block
        self.indexers = indexers
        self.shape = shape

    def __repr__(self) -> str:
        return f"{type(self).__name__}({repr(self.block)}, {self.indexers})"

    @cache_readonly
    def needs_filling(self) -> bool:
        for indexer in self.indexers.values():
            # FIXME: cache results of indexer == -1 checks.
            if (indexer == -1).any():
                return True

        return False

    @cache_readonly
    def dtype(self):
        blk = self.block
        if blk is None:
            raise AssertionError("Block is None, no dtype")

        if not self.needs_filling:
            return blk.dtype
        return ensure_dtype_can_hold_na(blk.dtype)

    def is_valid_na_for(self, dtype: DtypeObj) -> bool:
        """
        Check that we are all-NA of a type/dtype that is compatible with this dtype.
        Augments `self.is_na` with an additional check of the type of NA values.
        """
        if not self.is_na:
            return False
        if self.block is None:
            return True

        if self.dtype == object:
            values = self.block.values
            return all(is_valid_na_for_dtype(x, dtype) for x in values.ravel(order="K"))

        if self.dtype.kind == dtype.kind == "M" and not is_dtype_equal(
            self.dtype, dtype
        ):
            # fill_values match but we should not cast self.block.values to dtype
            return False

        na_value = self.block.fill_value
        return is_valid_na_for_dtype(na_value, dtype)

    @cache_readonly
    def is_na(self) -> bool:
        if self.block is None:
            return True

        if not self.block._can_hold_na:
            return False

        # Usually it's enough to check but a small fraction of values to see if
        # a block is NOT null, chunks should help in such cases.  1000 value
        # was chosen rather arbitrarily.
        values = self.block.values
        if is_sparse(self.block.values.dtype):
            return False
        elif self.block.is_extension:
            # TODO(EA2D): no need for special case with 2D EAs
            values_flat = values
        else:
            values_flat = values.ravel(order="K")

        return isna_all(values_flat)

    def get_reindexed_values(self, empty_dtype: DtypeObj, upcasted_na) -> ArrayLike:
        if upcasted_na is None:
            # No upcasting is necessary
            fill_value = self.block.fill_value
            values = self.block.get_values()
        else:
            fill_value = upcasted_na

            if self.is_valid_na_for(empty_dtype):
                blk_dtype = getattr(self.block, "dtype", None)

                # error: Value of type variable "_DTypeScalar" of "dtype" cannot be
                # "object"
                if blk_dtype == np.dtype(object):  # type: ignore[type-var]
                    # we want to avoid filling with np.nan if we are
                    # using None; we already know that we are all
                    # nulls
                    values = self.block.values.ravel(order="K")
                    if len(values) and values[0] is None:
                        fill_value = None

                if is_datetime64tz_dtype(empty_dtype):
                    # TODO(EA2D): special case unneeded with 2D EAs
                    i8values = np.full(self.shape[1], fill_value.value)
                    return DatetimeArray(i8values, dtype=empty_dtype)
                elif is_extension_array_dtype(blk_dtype):
                    pass
                elif is_extension_array_dtype(empty_dtype):
                    # error: Item "dtype[Any]" of "Union[dtype[Any], ExtensionDtype]"
                    # has no attribute "construct_array_type"
                    cls = empty_dtype.construct_array_type()  # type: ignore[union-attr]
                    missing_arr = cls._from_sequence([], dtype=empty_dtype)
                    ncols, nrows = self.shape
                    assert ncols == 1, ncols
                    empty_arr = -1 * np.ones((nrows,), dtype=np.intp)
                    return missing_arr.take(
                        empty_arr, allow_fill=True, fill_value=fill_value
                    )
                else:
                    # NB: we should never get here with empty_dtype integer or bool;
                    #  if we did, the missing_arr.fill would cast to gibberish

                    # error: Argument "dtype" to "empty" has incompatible type
                    # "Union[dtype[Any], ExtensionDtype]"; expected "Union[dtype[Any],
                    # None, type, _SupportsDType, str, Union[Tuple[Any, int], Tuple[Any,
                    # Union[int, Sequence[int]]], List[Any], _DTypeDict, Tuple[Any,
                    # Any]]]"
                    missing_arr = np.empty(
                        self.shape, dtype=empty_dtype  # type: ignore[arg-type]
                    )
                    missing_arr.fill(fill_value)
                    return missing_arr

            if (not self.indexers) and (not self.block._can_consolidate):
                # preserve these for validation in concat_compat
                return self.block.values

            if self.block.is_bool and not self.block.is_categorical:
                # External code requested filling/upcasting, bool values must
                # be upcasted to object to avoid being upcasted to numeric.
                values = self.block.astype(np.object_).values
            elif self.block.is_extension:
                values = self.block.values
            else:
                # No dtype upcasting is done here, it will be performed during
                # concatenation itself.
                values = self.block.values

        if not self.indexers:
            # If there's no indexing to be done, we want to signal outside
            # code that this array must be copied explicitly.  This is done
            # by returning a view and checking `retval.base`.
            values = values.view()

        else:
            for ax, indexer in self.indexers.items():
                values = algos.take_nd(values, indexer, axis=ax)

        return values


def _concatenate_join_units(
    join_units: List[JoinUnit], concat_axis: int, copy: bool
) -> ArrayLike:
    """
    Concatenate values from several join units along selected axis.
    """
    if concat_axis == 0 and len(join_units) > 1:
        # Concatenating join units along ax0 is handled in _merge_blocks.
        raise AssertionError("Concatenating join units along axis0")

    empty_dtype = _get_empty_dtype(join_units)

    has_none_blocks = any(unit.block is None for unit in join_units)
    upcasted_na = _dtype_to_na_value(empty_dtype, has_none_blocks)

    to_concat = [
        ju.get_reindexed_values(empty_dtype=empty_dtype, upcasted_na=upcasted_na)
        for ju in join_units
    ]

    if len(to_concat) == 1:
        # Only one block, nothing to concatenate.
        concat_values = to_concat[0]
        if copy:
            if isinstance(concat_values, np.ndarray):
                # non-reindexed (=not yet copied) arrays are made into a view
                # in JoinUnit.get_reindexed_values
                if concat_values.base is not None:
                    concat_values = concat_values.copy()
            else:
                concat_values = concat_values.copy()
    elif any(isinstance(t, ExtensionArray) for t in to_concat):
        # concatting with at least one EA means we are concatting a single column
        # the non-EA values are 2D arrays with shape (1, n)
        to_concat = [t if isinstance(t, ExtensionArray) else t[0, :] for t in to_concat]
        concat_values = concat_compat(to_concat, axis=0, ea_compat_axis=True)
        concat_values = ensure_block_shape(concat_values, 2)

    else:
        concat_values = concat_compat(to_concat, axis=concat_axis)

    return concat_values


def _dtype_to_na_value(dtype: DtypeObj, has_none_blocks: bool):
    """
    Find the NA value to go with this dtype.
    """
    if is_extension_array_dtype(dtype):
        # error: Item "dtype[Any]" of "Union[dtype[Any], ExtensionDtype]" has no
        # attribute "na_value"
        return dtype.na_value  # type: ignore[union-attr]
    elif dtype.kind in ["m", "M"]:
        return dtype.type("NaT")
    elif dtype.kind in ["f", "c"]:
        return dtype.type("NaN")
    elif dtype.kind == "b":
        return None
    elif dtype.kind in ["i", "u"]:
        if not has_none_blocks:
            return None
        return np.nan
    elif dtype.kind == "O":
        return np.nan
    raise NotImplementedError


def _get_empty_dtype(join_units: Sequence[JoinUnit]) -> DtypeObj:
    """
    Return dtype and N/A values to use when concatenating specified units.

    Returned N/A value may be None which means there was no casting involved.

    Returns
    -------
    dtype
    """
    if len(join_units) == 1:
        blk = join_units[0].block
        if blk is None:
            return np.dtype(np.float64)

    if _is_uniform_reindex(join_units):
        # FIXME: integrate property
        empty_dtype = join_units[0].block.dtype
        return empty_dtype

    has_none_blocks = any(unit.block is None for unit in join_units)

    dtypes = [
        unit.dtype for unit in join_units if unit.block is not None and not unit.is_na
    ]
    if not len(dtypes):
        dtypes = [unit.dtype for unit in join_units if unit.block is not None]

    dtype = find_common_type(dtypes)
    if has_none_blocks:
        dtype = ensure_dtype_can_hold_na(dtype)
    return dtype


def _is_uniform_join_units(join_units: List[JoinUnit]) -> bool:
    """
    Check if the join units consist of blocks of uniform type that can
    be concatenated using Block.concat_same_type instead of the generic
    _concatenate_join_units (which uses `concat_compat`).

    """
    # TODO: require dtype match in addition to same type?  e.g. DatetimeTZBlock
    #  cannot necessarily join
    return (
        # all blocks need to have the same type
        all(type(ju.block) is type(join_units[0].block) for ju in join_units)  # noqa
        and
        # no blocks that would get missing values (can lead to type upcasts)
        # unless we're an extension dtype.
        all(not ju.is_na or ju.block.is_extension for ju in join_units)
        and
        # no blocks with indexers (as then the dimensions do not fit)
        all(not ju.indexers for ju in join_units)
        and
        # only use this path when there is something to concatenate
        len(join_units) > 1
    )


def _is_uniform_reindex(join_units) -> bool:
    return (
        # TODO: should this be ju.block._can_hold_na?
        all(ju.block and ju.block.is_extension for ju in join_units)
        and len({ju.block.dtype.name for ju in join_units}) == 1
    )


def _trim_join_unit(join_unit: JoinUnit, length: int) -> JoinUnit:
    """
    Reduce join_unit's shape along item axis to length.

    Extra items that didn't fit are returned as a separate block.
    """
    if 0 not in join_unit.indexers:
        extra_indexers = join_unit.indexers

        if join_unit.block is None:
            extra_block = None
        else:
            extra_block = join_unit.block.getitem_block(slice(length, None))
            join_unit.block = join_unit.block.getitem_block(slice(length))
    else:
        extra_block = join_unit.block

        extra_indexers = copy.copy(join_unit.indexers)
        extra_indexers[0] = extra_indexers[0][length:]
        join_unit.indexers[0] = join_unit.indexers[0][:length]

    extra_shape = (join_unit.shape[0] - length,) + join_unit.shape[1:]
    join_unit.shape = (length,) + join_unit.shape[1:]

    return JoinUnit(block=extra_block, indexers=extra_indexers, shape=extra_shape)


def _combine_concat_plans(plans, concat_axis: int):
    """
    Combine multiple concatenation plans into one.

    existing_plan is updated in-place.
    """
    if len(plans) == 1:
        for p in plans[0]:
            yield p[0], [p[1]]

    elif concat_axis == 0:
        offset = 0
        for plan in plans:
            last_plc = None

            for plc, unit in plan:
                yield plc.add(offset), [unit]
                last_plc = plc

            if last_plc is not None:
                offset += last_plc.as_slice.stop

    else:
        num_ended = [0]

        def _next_or_none(seq):
            retval = next(seq, None)
            if retval is None:
                num_ended[0] += 1
            return retval

        plans = list(map(iter, plans))
        next_items = list(map(_next_or_none, plans))

        while num_ended[0] != len(next_items):
            if num_ended[0] > 0:
                raise ValueError("Plan shapes are not aligned")

            placements, units = zip(*next_items)

            lengths = list(map(len, placements))
            min_len, max_len = min(lengths), max(lengths)

            if min_len == max_len:
                yield placements[0], units
                next_items[:] = map(_next_or_none, plans)
            else:
                yielded_placement = None
                yielded_units = [None] * len(next_items)
                for i, (plc, unit) in enumerate(next_items):
                    yielded_units[i] = unit
                    if len(plc) > min_len:
                        # _trim_join_unit updates unit in place, so only
                        # placement needs to be sliced to skip min_len.
                        next_items[i] = (plc[min_len:], _trim_join_unit(unit, min_len))
                    else:
                        yielded_placement = plc
                        next_items[i] = _next_or_none(plans[i])

                yield yielded_placement, yielded_units
