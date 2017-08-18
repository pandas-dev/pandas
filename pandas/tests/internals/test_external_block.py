# -*- coding: utf-8 -*-
# pylint: disable=W0102

import numpy as np

import pandas as pd
from pandas.core.internals import Block, BlockManager, SingleBlockManager


class CustomBlock(Block):

    def formatting_values(self):
        return np.array(["Val: {}".format(i) for i in self.values])


def test_custom_repr():
    values = np.arange(3, dtype='int64')

    # series
    block = CustomBlock(values, placement=slice(0, 3))

    s = pd.Series(SingleBlockManager(block, pd.RangeIndex(3)))
    assert repr(s) == '0    Val: 0\n1    Val: 1\n2    Val: 2\ndtype: int64'

    # dataframe
    block = CustomBlock(values.reshape(1, -1), placement=slice(0, 1))
    blk_mgr = BlockManager([block], [['col'], range(3)])
    df = pd.DataFrame(blk_mgr)
    assert repr(df) == '      col\n0  Val: 0\n1  Val: 1\n2  Val: 2'
