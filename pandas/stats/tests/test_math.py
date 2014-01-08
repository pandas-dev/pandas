import nose

from datetime import datetime
from numpy.random import randn
import numpy as np

from pandas.core.api import Series, DataFrame, date_range
from pandas.util.testing import assert_almost_equal
import pandas.core.datetools as datetools
import pandas.stats.moments as mom
import pandas.util.testing as tm
import pandas.stats.math as pmath
import pandas.tests.test_series as ts
from pandas import ols

N, K = 100, 10

_have_statsmodels = True
try:
    import statsmodels.api as sm
except ImportError:
    try:
        import scikits.statsmodels.api as sm
    except ImportError:
        _have_statsmodels = False


class TestMath(tm.TestCase):

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def setUp(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = date_range(datetime(2009, 1, 1), periods=N)

        self.series = Series(arr.copy(), index=self.rng)

        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))

    def test_rank_1d(self):
        self.assertEqual(1, pmath.rank(self.series))
        self.assertEqual(0, pmath.rank(Series(0, self.series.index)))

    def test_solve_rect(self):
        if not _have_statsmodels:
            raise nose.SkipTest("no statsmodels")

        b = Series(np.random.randn(N), self.frame.index)
        result = pmath.solve(self.frame, b)
        expected = ols(y=b, x=self.frame, intercept=False).beta
        self.assert_(np.allclose(result, expected))

    def test_inv_illformed(self):
        singular = DataFrame(np.array([[1, 1], [2, 2]]))
        rs = pmath.inv(singular)
        expected = np.array([[0.1, 0.2], [0.1, 0.2]])
        self.assert_(np.allclose(rs, expected))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
