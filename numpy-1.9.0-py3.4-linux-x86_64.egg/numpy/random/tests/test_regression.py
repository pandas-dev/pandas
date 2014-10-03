from __future__ import division, absolute_import, print_function

from numpy.testing import (TestCase, run_module_suite, assert_,
                           assert_array_equal)
from numpy import random
from numpy.compat import long
import numpy as np


class TestRegression(TestCase):

    def test_VonMises_range(self):
        # Make sure generated random variables are in [-pi, pi].
        # Regression test for ticket #986.
        for mu in np.linspace(-7., 7., 5):
            r = random.mtrand.vonmises(mu, 1, 50)
            assert_(np.all(r > -np.pi) and np.all(r <= np.pi))

    def test_hypergeometric_range(self):
        # Test for ticket #921
        assert_(np.all(np.random.hypergeometric(3, 18, 11, size=10) < 4))
        assert_(np.all(np.random.hypergeometric(18, 3, 11, size=10) > 0))

    def test_logseries_convergence(self):
        # Test for ticket #923
        N = 1000
        np.random.seed(0)
        rvsn = np.random.logseries(0.8, size=N)
        # these two frequency counts should be close to theoretical
        # numbers with this large sample
        # theoretical large N result is 0.49706795
        freq = np.sum(rvsn == 1) / float(N)
        msg = "Frequency was %f, should be > 0.45" % freq
        assert_(freq > 0.45, msg)
        # theoretical large N result is 0.19882718
        freq = np.sum(rvsn == 2) / float(N)
        msg = "Frequency was %f, should be < 0.23" % freq
        assert_(freq < 0.23, msg)

    def test_permutation_longs(self):
        np.random.seed(1234)
        a = np.random.permutation(12)
        np.random.seed(1234)
        b = np.random.permutation(long(12))
        assert_array_equal(a, b)

    def test_randint_range(self):
        # Test for ticket #1690
        lmax = np.iinfo('l').max
        lmin = np.iinfo('l').min
        try:
            random.randint(lmin, lmax)
        except:
            raise AssertionError

    def test_shuffle_mixed_dimension(self):
        # Test for trac ticket #2074
        for t in [[1, 2, 3, None],
                  [(1, 1), (2, 2), (3, 3), None],
                  [1, (2, 2), (3, 3), None],
                  [(1, 1), 2, 3, None]]:
            np.random.seed(12345)
            shuffled = list(t)
            random.shuffle(shuffled)
            assert_array_equal(shuffled, [t[0], t[3], t[1], t[2]])

    def test_call_within_randomstate(self):
        # Check that custom RandomState does not call into global state
        m = np.random.RandomState()
        res = np.array([0, 8, 7, 2, 1, 9, 4, 7, 0, 3])
        for i in range(3):
            np.random.seed(i)
            m.seed(4321)
            # If m.state is not honored, the result will change
            assert_array_equal(m.choice(10, size=10, p=np.ones(10)/10.), res)

    def test_multivariate_normal_size_types(self):
        # Test for multivariate_normal issue with 'size' argument.
        # Check that the multivariate_normal size argument can be a
        # numpy integer.
        np.random.multivariate_normal([0], [[0]], size=1)
        np.random.multivariate_normal([0], [[0]], size=np.int_(1))
        np.random.multivariate_normal([0], [[0]], size=np.int64(1))

if __name__ == "__main__":
    run_module_suite()
