# "magictoken" is used for markers as beginning and ending of example text.

import unittest
import numpy as np
import numba


class NumpyGeneratorUsageTest(unittest.TestCase):

    def test_numpy_gen_usage(self):
        # magictoken.npgen_usage.begin
        x = np.random.default_rng(1)
        y = np.random.default_rng(1)

        size = 10

        @numba.njit
        def do_stuff(gen):
            return gen.random(size=int(size / 2))

        original = x.random(size=size)
        # [0.51182162 0.9504637  0.14415961 0.94864945 0.31183145
        #  0.42332645 0.82770259 0.40919914 0.54959369 0.02755911]

        numba_func_res = do_stuff(y)
        # [0.51182162 0.9504637  0.14415961 0.94864945 0.31183145]

        after_numba = y.random(size=int(size / 2))
        # [0.42332645 0.82770259 0.40919914 0.54959369 0.02755911]

        # magictoken.npgen_usage.end
        numba_res = np.concatenate((numba_func_res, after_numba))
        for _np_res, _nb_res in zip(original, numba_res):
            self.assertEqual(_np_res, _nb_res)


if __name__ == '__main__':
    unittest.main()
