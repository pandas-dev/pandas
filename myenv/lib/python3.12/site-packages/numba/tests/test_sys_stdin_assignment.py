import sys

import unittest
from numba import njit


@njit
def f0(a, b):
    return a + b


@njit
def f1(begin1, end1, begin2, end2):
    if begin1 > begin2: return f1(begin2, end2, begin1, end1)
    return end1 + 1 >= begin2


@njit
def f0_2(a, b):
    return a + b


@njit
def f1_2(begin1, end1, begin2, end2):
    if begin1 > begin2: return f1_2(begin2, end2, begin1, end1)
    return end1 + 1 >= begin2


class TestSysStdinAssignment(unittest.TestCase):

    def test_no_reassignment_of_stdout(self):
        """
        https://github.com/numba/numba/issues/3027
        Older versions of colorama break stdout/err when recursive functions
        are compiled.

        This test should work irrespective of colorama version, or indeed its
        presence. If the version is too low, it should be disabled and the test
        should work anyway, if it is a sufficiently high version or it is not
        present, it should work anyway.
        """

        originally = sys.stdout, sys.stderr

        try:
            sys.stdout = None
            f0(0, 1)  # Not changed stdout?
            self.assertEqual(sys.stdout, None)
            f1(0, 1, 2, 3)  # Not changed stdout?
            self.assertEqual(sys.stdout, None)

            sys.stderr = None
            f0_2(0, 1)  # Not changed stderr?
            self.assertEqual(sys.stderr, None)
            f1_2(0, 1, 2, 3)  # Not changed stderr?
            self.assertEqual(sys.stderr, None)

        finally:
            sys.stdout, sys.stderr = originally

        self.assertNotEqual(sys.stderr, None)
        self.assertNotEqual(sys.stdout, None)


if __name__ == '__main__':
    unittest.main()

