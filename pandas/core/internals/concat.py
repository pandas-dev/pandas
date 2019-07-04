# TODO: Needs a better name; too many modules are already called "concat"
from collections import defaultdict
import copy

import numpy as np

from pandas._libs import internals as libinternals, tslibs
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.cast import maybe_promote
from pandas.core.dtypes.common import (
    _get_dtype,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_extension_array_dtype,
    is_float_dtype,
    is_numeric_dtype,
    is_sparse,
    is_timedelta64_dtype,
)
import pandas.core.dtypes.concat as _concat
from pandas.core.dtypes.missing import isna

import pandas.core.algorithms as algos


def get_mgr_concatenation_plan(mgr, indexers):
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
    mgr_shape = list(mgr.shape)
    for ax, indexer in indexers.items():
        mgr_shape[ax] = len(indexer)
    mgr_shape = tuple(mgr_shape)

    if 0 in indexers:
        ax0_indexer = indexers.pop(0)
        blknos = algos.take_1d(mgr._blknos, ax0_indexer, fill_value=-1)
        blklocs = algos.take_1d(mgr._blklocs, ax0_indexer, fill_value=-1)
    else:

        if mgr._is_single_block:
            blk = mgr.blocks[0]
            return [(blk.mgr_locs, JoinUnit(blk, mgr_shape, indexers))]

        ax0_indexer = None
        blknos = mgr._blknos
        blklocs = mgr._blklocs

    plan = []
    for blkno, placements in libinternals.get_blkno_placements(
        blknos, mgr.nblocks, group=False
    ):

        assert placements.is_slice_like

        join_unit_indexers = indexers.copy()

        shape = list(mgr_shape)
        shape[0] = len(placements)
        shape = tuple(shape)

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
    def __init__(self, block, shape, indexers=None):
        # Passing shape explicitly is required for cases when block is None.
        if indexers is None:
            indexers = {}
        self.block = block
        self.indexers = indexers
        self.shape = shape

    def __repr__(self):
        return "{name}({block!r}, {indexers})".format(
            name=self.__class__.__name__, block=self.block, indexers=self.indexers
        )

    @cache_readonly
    def needs_filling(self):
        for indexer in self.indexers.values():
            # FIXME: cache results of indexer == -1 checks.
            if (indexer == -1).any():
                return True

        return False

    @cache_readonly
    def dtype(self):
        if self.block is None:
            raise AssertionError("Block is None, no dtype")

        if not self.needs_filling:
            return self.block.dtype
        else:
            return _get_dtype(maybe_promote(self.block.dtype, self.block.fill_value)[0])

    @cache_readonly
    def is_na(self):
        if self.block is None:
            return True

        if not self.block._can_hold_na:
            return False

        # Usually it's enough to check but a small fraction of values to see if
        # a block is NOT null, chunks should help in such cases.  1000 value
        # was chosen rather arbitrarily.
        values = self.block.values
        if self.block.is_categorical:
            values_flat = values.categories
        elif is_sparse(self.block.values.dtype):
            return False
        elif self.block.is_extension:
            values_flat = values
        else:
            values_flat = values.ravel(order="K")
        total_len = values_flat.shape[0]
        chunk_len = max(total_len // 40, 1000)
        for i in range(0, total_len, chunk_len):
            if not isna(values_flat[i : i + chunk_len]).all():
                return False

        return True

    def get_reindexed_values(self, empty_dtype, upcasted_na):
        if upcasted_na is None:
            # No upcasting is necessary
            fill_value = self.block.fill_value
            values = self.block.get_values()
        else:
            fill_value = upcasted_na

            if self.is_na:
                if getattr(self.block, "is_object", False):
                    # we want to avoid filling with np.nan if we are
                    # using None; we already know that we are all
                    # nulls
                    values = self.block.values.ravel(order="K")
                    if len(values) and values[0] is None:
                        fill_value = None

                if getattr(self.block, "is_datetimetz", False) or is_datetime64tz_dtype(
                    empty_dtype
                ):
                    if self.block is None:
                        array = empty_dtype.construct_array_type()
                        return array(
                            np.full(self.shape[1], fill_value.value), dtype=empty_dtype
                        )
                    pass
                elif getattr(self.block, "is_categorical", False):
                    pass
                elif getattr(self.block, "is_extension", False):
                    pass
                else:
                    missing_arr = np.empty(self.shape, dtype=empty_dtype)
                    missing_arr.fill(fill_value)
                    return missing_arr

            if not self.indexers:
                if not self.block._can_consolidate:
                    # preserve these for validation in _concat_compat
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
                values = self.block.get_values()

        if not self.indexers:
            # If there's no indexing to be done, we want to signal outside
            # code that this array must be copied explicitly.  This is done
            # by returning a view and checking `retval.base`.
            values = values.view()

        else:
            for ax, indexer in self.indexers.items():
                values = algos.take_nd(values, indexer, axis=ax, fill_value=fill_value)

        return values


def concatenate_join_units(join_units, concat_axis, copy):
    """
    Concatenate values from several join units along selected axis.
    """
    if concat_axis == 0 and len(join_units) > 1:
        # Concatenating join units along ax0 is handled in _merge_blocks.
        raise AssertionError("Concatenating join units along axis0")

    empty_dtype, upcasted_na = get_empty_dtype_and_na(join_units)

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
    else:
        concat_values = _concat._concat_compat(to_concat, axis=concat_axis)

    return concat_values


def get_empty_dtype_and_na(join_units):
    """
    Return dtype and N/A values to use when concatenating specified units.

    Returned N/A value may be None which means there was no casting involved.

    Returns
    -------
    dtype
    na
    """
    if len(join_units) == 1:
        blk = join_units[0].block
        if blk is None:
            return np.float64, np.nan

    if is_uniform_reindex(join_units):
        # XXX: integrate property
        empty_dtype = join_units[0].block.dtype
        upcasted_na = join_units[0].block.fill_value
        return empty_dtype, upcasted_na

    has_none_blocks = False
    dtypes = [None] * len(join_units)
    for i, unit in enumerate(join_units):
        if unit.block is None:
            has_none_blocks = True
        else:
            dtypes[i] = unit.dtype

    upcast_classes = defaultdict(list)
    null_upcast_classes = defaultdict(list)
    for dtype, unit in zip(dtypes, join_units):
        if dtype is None:
            continue

        if is_categorical_dtype(dtype):
            upcast_cls = "category"
        elif is_datetime64tz_dtype(dtype):
            upcast_cls = "datetimetz"
        elif issubclass(dtype.type, np.bool_):
            upcast_cls = "bool"
        elif issubclass(dtype.type, np.object_):
            upcast_cls = "object"
        elif is_datetime64_dtype(dtype):
            upcast_cls = "datetime"
        elif is_timedelta64_dtype(dtype):
            upcast_cls = "timedelta"
        elif is_sparse(dtype):
            upcast_cls = dtype.subtype.name
        elif is_extension_array_dtype(dtype):
            upcast_cls = "object"
        elif is_float_dtype(dtype) or is_numeric_dtype(dtype):
            upcast_cls = dtype.name
        else:
            upcast_cls = "float"

        # Null blocks should not influence upcast class selection, unless there
        # are only null blocks, when same upcasting rules must be applied to
        # null upcast classes.
        if unit.is_na:
            null_upcast_classes[upcast_cls].append(dtype)
        else:
            upcast_classes[upcast_cls].append(dtype)

    if not upcast_classes:
        upcast_classes = null_upcast_classes

    # create the result
    if "object" in upcast_classes:
        return np.dtype(np.object_), np.nan
    elif "bool" in upcast_classes:
        if has_none_blocks:
            return np.dtype(np.object_), np.nan
        else:
            return np.dtype(np.bool_), None
    elif "category" in upcast_classes:
        return np.dtype(np.object_), np.nan
    elif "datetimetz" in upcast_classes:
        # GH-25014. We use NaT instead of iNaT, since this eventually
        # ends up in DatetimeArray.take, which does not allow iNaT.
        dtype = upcast_classes["datetimetz"]
        return dtype[0], tslibs.NaT
    elif "datetime" in upcast_classes:
        return np.dtype("M8[ns]"), tslibs.iNaT
    elif "timedelta" in upcast_classes:
        return np.dtype("m8[ns]"), tslibs.iNaT
    else:  # pragma
        try:
            g = np.find_common_type(upcast_classes, [])
        except TypeError:
            # At least one is an ExtensionArray
            return np.dtype(np.object_), np.nan
        else:
            if is_float_dtype(g):
                return g, g.type(np.nan)
            elif is_numeric_dtype(g):
                if has_none_blocks:
                    return np.float64, np.nan
                else:
                    return g, None

    msg = "invalid dtype determination in get_concat_dtype"
    raise AssertionError(msg)


def is_uniform_join_units(join_units):
    """
    Check if the join units consist of blocks of uniform type that can
    be concatenated using Block.concat_same_type instead of the generic
    concatenate_join_units (which uses `_concat._concat_compat`).

    """
    return (
        # all blocks need to have the same type
        all(type(ju.block) is type(join_units[0].block) for ju in join_units)
        and  # noqa
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


def is_uniform_reindex(join_units):
    return (
        # TODO: should this be ju.block._can_hold_na?
        all(ju.block and ju.block.is_extension for ju in join_units)
        and len({ju.block.dtype.name for ju in join_units}) == 1
    )


def trim_join_unit(join_unit, length):
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


def combine_concat_plans(plans, concat_axis):
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
                        # trim_join_unit updates unit in place, so only
                        # placement needs to be sliced to skip min_len.
                        next_items[i] = (plc[min_len:], trim_join_unit(unit, min_len))
                    else:
                        yielded_placement = plc
                        next_items[i] = _next_or_none(plans[i])

                yield yielded_placement, yielded_units
