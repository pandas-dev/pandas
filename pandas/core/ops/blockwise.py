from typing import TYPE_CHECKING, List, Tuple

import numpy as np

from pandas._typing import ArrayLike

if TYPE_CHECKING:
    from pandas.core.internals.blocks import Block  # noqa:F401


def operate_blockwise(left, right, array_op):
    assert right._indexed_same(left)

    def get_same_shape_values(
        lblk: "Block", rblk: "Block", left_ea: bool, right_ea: bool
    ) -> Tuple[ArrayLike, ArrayLike]:
        """
        Slice lblk.values to align with rblk.  Squeeze if we have EAs.
        """
        lvals = lblk.values
        rvals = rblk.values

        # TODO(EA2D): with 2D EAs pnly this first clause would be needed
        if not (left_ea or right_ea):
            lvals = lvals[rblk.mgr_locs.indexer, :]
            assert lvals.shape == rvals.shape, (lvals.shape, rvals.shape)
        elif left_ea and right_ea:
            assert lvals.shape == rvals.shape, (lvals.shape, rvals.shape)
        elif right_ea:
            # lvals are 2D, rvals are 1D
            lvals = lvals[rblk.mgr_locs.indexer, :]
            assert lvals.shape[0] == 1, lvals.shape
            lvals = lvals[0, :]
        else:
            # lvals are 1D, rvals are 2D
            assert rvals.shape[0] == 1, rvals.shape
            rvals = rvals[0, :]

        return lvals, rvals

    res_blks: List["Block"] = []
    rmgr = right._mgr
    for n, blk in enumerate(left._mgr.blocks):
        locs = blk.mgr_locs
        blk_vals = blk.values

        left_ea = not isinstance(blk_vals, np.ndarray)

        rblks = rmgr._slice_take_blocks_ax0(locs.indexer, only_slice=True)

        # Assertions are disabled for performance, but should hold:
        # if left_ea:
        #    assert len(locs) == 1, locs
        #    assert len(rblks) == 1, rblks
        #    assert rblks[0].shape[0] == 1, rblks[0].shape

        for k, rblk in enumerate(rblks):
            right_ea = not isinstance(rblk.values, np.ndarray)

            lvals, rvals = get_same_shape_values(blk, rblk, left_ea, right_ea)

            res_values = array_op(lvals, rvals)
            if left_ea and not right_ea and hasattr(res_values, "reshape"):
                res_values = res_values.reshape(1, -1)
            nbs = rblk._split_op_result(res_values)

            # Assertions are disabled for performance, but should hold:
            # if right_ea or left_ea:
            #    assert len(nbs) == 1
            # else:
            #    assert res_values.shape == lvals.shape, (res_values.shape, lvals.shape)

            for nb in nbs:
                # Reset mgr_locs to correspond to our original DataFrame
                nblocs = locs.as_array[nb.mgr_locs.indexer]
                nb.mgr_locs = nblocs
                # Assertions are disabled for performance, but should hold:
                #  assert len(nblocs) == nb.shape[0], (len(nblocs), nb.shape)
                #  assert all(x in locs.as_array for x in nb.mgr_locs.as_array)

            res_blks.extend(nbs)

    # Assertions are disabled for performance, but should hold:
    #  slocs = {y for nb in res_blks for y in nb.mgr_locs.as_array}
    #  nlocs = sum(len(nb.mgr_locs.as_array) for nb in res_blks)
    #  assert nlocs == len(left.columns), (nlocs, len(left.columns))
    #  assert len(slocs) == nlocs, (len(slocs), nlocs)
    #  assert slocs == set(range(nlocs)), slocs

    new_mgr = type(rmgr)(res_blks, axes=rmgr.axes, do_integrity_check=False)
    return new_mgr
