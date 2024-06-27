import itertools
import unittest
from numba import njit
from numba.core import types


def template(fromty, toty):
    def closure(self):
        def cast(x):
            y = x
            return y

        cfunc = njit(toty(fromty))(cast)
        self.assertAlmostEqual(cfunc(1), 1)

    return closure


class TestNumberConversion(unittest.TestCase):
    """
    Test all int/float numeric conversion to ensure we have all the external
    dependencies to perform these conversions.
    """
    # NOTE: more implicit tests are in test_numberctor

    @classmethod
    def automatic_populate(cls):
        tys = types.integer_domain | types.real_domain
        for fromty, toty in itertools.permutations(tys, r=2):
            test_name = "test_{fromty}_to_{toty}".format(fromty=fromty,
                                                         toty=toty)
            setattr(cls, test_name, template(fromty, toty))


TestNumberConversion.automatic_populate()

if __name__ == '__main__':
    unittest.main()
