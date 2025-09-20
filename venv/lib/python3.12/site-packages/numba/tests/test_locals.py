from numba import float32, njit
import unittest


def foo():
    x = 123
    return x


class TestLocals(unittest.TestCase):

    def test_seed_types(self):
        cfunc = njit((), locals={'x': float32})(foo)
        self.assertEqual(cfunc.nopython_signatures[0].return_type, float32)


if __name__ == '__main__':
    unittest.main()
