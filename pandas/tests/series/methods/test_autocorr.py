import numpy as np

import pandas._testing as tm


class TestAutoCorr:
    def test_autocorr(self, datetime_series):
        # Just run the function
        corr1 = datetime_series.autocorr()

        # Now run it with the lag parameter
        corr2 = datetime_series.autocorr(lag=1)

        # corr() with lag needs Series of at least length 2
        if len(datetime_series) <= 2:
            assert np.isnan(corr1)
            assert np.isnan(corr2)
        else:
            # Not exact equality: the dot-product Pearson kernel (GH#65502) is not
            # bit-reproducible across calls on threaded BLAS (e.g. Accelerate), so
            # compare with a tolerance rather than ``==``.
            tm.assert_almost_equal(corr1, corr2)

        # Choose a random lag between 1 and length of Series - 2
        # and compare the result with the Series corr() function
        n = 1 + np.random.default_rng(2).integers(max(1, len(datetime_series) - 2))
        corr1 = datetime_series.corr(datetime_series.shift(n))
        corr2 = datetime_series.autocorr(lag=n)

        # corr() with lag needs Series of at least length 2
        if len(datetime_series) <= 2:
            assert np.isnan(corr1)
            assert np.isnan(corr2)
        else:
            tm.assert_almost_equal(corr1, corr2)
