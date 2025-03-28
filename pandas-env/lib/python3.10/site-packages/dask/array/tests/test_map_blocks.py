from __future__ import annotations

import dask.array as da
from dask.base import collections_to_dsk


def test_map_blocks_block_id_fusion():
    arr = da.ones((20, 10), chunks=(2, 5))

    def dummy(x, block_id=None, block_info=None):
        return x

    result = arr.map_blocks(dummy).astype("f8")
    dsk = collections_to_dsk([result])
    assert len(dsk) == 20
