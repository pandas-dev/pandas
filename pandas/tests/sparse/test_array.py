import numpy as np

import pandas as pd
from pandas.util import testing as tm


class TestSparseArray(object):

    def test_nonzero(self):
        sa = pd.SparseArray([
            float('nan'),
            float('nan'),
            1, 0, 0,
            2, 0, 0, 0,
            3, 0, 0
        ])
        tm.assert_numpy_array_equal(np.array([2, 5, 9], dtype=np.int32),
                                    sa.nonzero()[0])

        sa = pd.SparseArray(
            [0, 0, 1, 0, 0, 2, 0, 0, 0, 3, 0, 0])
        tm.assert_numpy_array_equal(np.array([2, 5, 9], dtype=np.int32),
                                    sa.nonzero()[0])
