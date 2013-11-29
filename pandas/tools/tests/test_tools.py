from pandas import DataFrame
from pandas.tools.describe import value_range

import numpy as np
import pandas.util.testing as tm


class TestTools(tm.TestCase):

    def test_value_range(self):
        df = DataFrame(np.random.randn(5, 5))
        df.ix[0, 2] = -5
        df.ix[2, 0] = 5

        res = value_range(df)

        self.assert_(res['Minimum'] == -5)
        self.assert_(res['Maximum'] == 5)

        df.ix[0, 1] = np.NaN

        self.assert_(res['Minimum'] == -5)
        self.assert_(res['Maximum'] == 5)
