from __future__ import annotations

import dask.array._array_expr._backends
from dask.array._array_expr import random
from dask.array._array_expr._collection import (
    Array,
    array,
    asanyarray,
    asarray,
    blockwise,
    concatenate,
    elemwise,
    from_array,
    rechunk,
    stack,
)
from dask.array._array_expr._creation import arange, empty, linspace, ones, zeros
from dask.array._array_expr._map_blocks import map_blocks
from dask.array._array_expr._reductions import _tree_reduce, reduction
