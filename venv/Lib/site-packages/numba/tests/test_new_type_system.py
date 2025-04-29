import numpy as np

from numba import njit, config
from numba.tests.support import TestCase


class TestTypes(TestCase):

    def setUp(self) -> None:
        if config.USE_LEGACY_TYPE_SYSTEM:
            self.skipTest("This test is only for the new type system")
        return super().setUp()

    def test_return_types(self):
        @njit
        def foo(x):
            return x

        cases = [
            # Python types
            1,
            1.2,
            (1 + 2j),
            True,
            # NumPy types
            np.int32(1),
            np.float64(1.2),
            np.complex64(1 + 2j),
            np.complex128(1 + 2j),
            np.bool_(True),
            np.datetime64('2020-01-01'),
            np.timedelta64(1, 'D'),
        ]

        for case in cases:
            self.assertEqual(foo(case), case)
            self.assertEqual(type(foo(case)), type(case))
