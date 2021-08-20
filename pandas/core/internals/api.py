"""
This is a pseudo-public API for downstream libraries.  We ask that downstream
authors

1) Try to avoid using internals directly altogether, and failing that,
2) Use only functions exposed here (or in core.internals)

"""
from __future__ import annotations

from collections import defaultdict
from typing import DefaultDict

import numpy as np

from pandas._libs.internals import BlockPlacement
from pandas._typing import (
    ArrayLike,
    Dtype,
)

from pandas.core.dtypes.common import (
    is_datetime64tz_dtype,
    pandas_dtype,
)

from pandas.core.arrays import DatetimeArray
from pandas.core.construction import extract_array
from pandas.core.indexes.api import Index
from pandas.core.internals.blocks import (
    Block,
    CategoricalBlock,
    DatetimeTZBlock,
    ExtensionBlock,
    check_ndim,
    ensure_block_shape,
    extract_pandas_array,
    get_block_type,
    maybe_coerce_values,
    new_block,
)
from pandas.core.internals.managers import (
    BlockManager,
    construction_error,
    multi_blockify,
    simple_blockify,
)


def make_block(
    values, placement, klass=None, ndim=None, dtype: Dtype | None = None
) -> Block:
    """
    This is a pseudo-public analogue to blocks.new_block.

    We ask that downstream libraries use this rather than any fully-internal
    APIs, including but not limited to:

    - core.internals.blocks.make_block
    - Block.make_block
    - Block.make_block_same_class
    - Block.__init__
    """
    if dtype is not None:
        dtype = pandas_dtype(dtype)

    values, dtype = extract_pandas_array(values, dtype, ndim)

    if klass is None:
        dtype = dtype or values.dtype
        klass = get_block_type(values, dtype)

    elif klass is DatetimeTZBlock and not is_datetime64tz_dtype(values.dtype):
        # pyarrow calls get here
        values = DatetimeArray._simple_new(values, dtype=dtype)

    if not isinstance(placement, BlockPlacement):
        placement = BlockPlacement(placement)

    ndim = maybe_infer_ndim(values, placement, ndim)
    if is_datetime64tz_dtype(values.dtype):
        # GH#41168 ensure we can pass 1D dt64tz values
        values = extract_array(values, extract_numpy=True)
        values = ensure_block_shape(values, ndim)

    check_ndim(values, placement, ndim)
    values = maybe_coerce_values(values)
    return klass(values, ndim=ndim, placement=placement)


def maybe_infer_ndim(values, placement: BlockPlacement, ndim: int | None) -> int:
    """
    If `ndim` is not provided, infer it from placment and values.
    """
    if ndim is None:
        # GH#38134 Block constructor now assumes ndim is not None
        if not isinstance(values.dtype, np.dtype):
            if len(placement) != 1:
                ndim = 1
            else:
                ndim = 2
        else:
            ndim = values.ndim
    return ndim


def create_block_manager_from_arrays(
    arrays,
    names: Index,
    axes: list[Index],
    consolidate: bool = True,
) -> BlockManager:
    # Assertions disabled for performance
    # assert isinstance(names, Index)
    # assert isinstance(axes, list)
    # assert all(isinstance(x, Index) for x in axes)

    arrays = [extract_array(x, extract_numpy=True) for x in arrays]

    try:
        blocks = _form_blocks(arrays, names, axes, consolidate)
        mgr = BlockManager(blocks, axes)
    except ValueError as e:
        raise construction_error(len(arrays), arrays[0].shape, axes, e)
    if consolidate:
        mgr._consolidate_inplace()
    return mgr


def _form_blocks(
    arrays: list[ArrayLike], names: Index, axes: list[Index], consolidate: bool
) -> list[Block]:
    # put "leftover" items in float bucket, where else?
    # generalize?
    items_dict: DefaultDict[str, list] = defaultdict(list)
    extra_locs = []

    names_idx = names
    if names_idx.equals(axes[0]):
        names_indexer = np.arange(len(names_idx))
    else:
        # Assertion disabled for performance
        # assert names_idx.intersection(axes[0]).is_unique
        names_indexer = names_idx.get_indexer_for(axes[0])

    for i, name_idx in enumerate(names_indexer):
        if name_idx == -1:
            extra_locs.append(i)
            continue

        v = arrays[name_idx]

        block_type = get_block_type(v)
        items_dict[block_type.__name__].append((i, v))

    blocks: list[Block] = []
    if len(items_dict["NumericBlock"]):
        numeric_blocks = multi_blockify(
            items_dict["NumericBlock"], consolidate=consolidate
        )
        blocks.extend(numeric_blocks)

    if len(items_dict["DatetimeLikeBlock"]):
        dtlike_blocks = multi_blockify(
            items_dict["DatetimeLikeBlock"], consolidate=consolidate
        )
        blocks.extend(dtlike_blocks)

    if len(items_dict["DatetimeTZBlock"]):
        dttz_blocks = [
            DatetimeTZBlock(
                ensure_block_shape(extract_array(array), 2),
                placement=BlockPlacement(i),
                ndim=2,
            )
            for i, array in items_dict["DatetimeTZBlock"]
        ]
        blocks.extend(dttz_blocks)

    if len(items_dict["ObjectBlock"]) > 0:
        object_blocks = simple_blockify(
            items_dict["ObjectBlock"], np.object_, consolidate=consolidate
        )
        blocks.extend(object_blocks)

    if len(items_dict["CategoricalBlock"]) > 0:
        cat_blocks = [
            CategoricalBlock(array, placement=BlockPlacement(i), ndim=2)
            for i, array in items_dict["CategoricalBlock"]
        ]
        blocks.extend(cat_blocks)

    if len(items_dict["ExtensionBlock"]):
        external_blocks = [
            ExtensionBlock(array, placement=BlockPlacement(i), ndim=2)
            for i, array in items_dict["ExtensionBlock"]
        ]

        blocks.extend(external_blocks)

    if len(extra_locs):
        shape = (len(extra_locs),) + tuple(len(x) for x in axes[1:])

        # empty items -> dtype object
        block_values = np.empty(shape, dtype=object)
        block_values.fill(np.nan)

        na_block = new_block(block_values, placement=extra_locs, ndim=2)
        blocks.append(na_block)

    return blocks
