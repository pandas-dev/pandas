import unittest
import numpy as np

import tables as tb


class TestAuxiliaryFunctions(unittest.TestCase):
    def test_keysort(self):
        N = 1000
        rnd = np.random.randint(N, size=N)
        for dtype1 in ('S6', 'b1', 'i1', 'i8', 'u4', 'u8', 'f4', 'f8'):
            for dtype2 in ('u4', 'i8'):
                a = np.array(rnd, dtype1)
                b = np.array(rnd, dtype2)

                c = a.copy()
                d = c.argsort()
                e = c[d]
                f = b[d]

                tb.indexesextension.keysort(a, b)
                self.assertTrue((a == e).all())
                self.assertTrue((b == f).all())


def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestAuxiliaryFunctions))
    return theSuite


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
