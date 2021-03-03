"""
This is a pseudo-public API for downstream libraries.  We ask that downstream
authors

1) Try to avoid using internals directly altogether, and failing that,
2) Use only functions exposed here (or in core.internals)

"""
from typing import Optional

import numpy as np

from pandas._typing import Dtype

from pandas.core.dtypes.common import is_datetime64tz_dtype
from pandas.core.dtypes.dtypes import PandasDtype
from pandas.core.dtypes.generic import ABCPandasArray

from pandas.core.arrays import DatetimeArray
from pandas.core.internals.blocks import (
    Block,
    DatetimeTZBlock,
    get_block_type,
)


def make_block(
    values, placement, klass=None, ndim=None, dtype: Optional[Dtype] = None
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
    if isinstance(values, ABCPandasArray):
        # Ensure that we don't allow PandasArray / PandasDtype in internals.
        # For now, blocks should be backed by ndarrays when possible.
        values = values.to_numpy()
        if ndim and ndim > 1:
            # TODO(EA2D): special case not needed with 2D EAs
            values = np.atleast_2d(values)

    if isinstance(dtype, PandasDtype):
        dtype = dtype.numpy_dtype

    if klass is None:
        dtype = dtype or values.dtype
        klass = get_block_type(values, dtype)

    elif klass is DatetimeTZBlock and not is_datetime64tz_dtype(values.dtype):
        # TODO: This is no longer hit internally; does it need to be retained
        #  for e.g. pyarrow?
        values = DatetimeArray._simple_new(values, dtype=dtype)

    return klass(values, ndim=ndim, placement=placement)
