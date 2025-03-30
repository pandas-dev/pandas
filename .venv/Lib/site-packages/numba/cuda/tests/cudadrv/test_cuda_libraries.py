from numba.cuda.testing import unittest
from numba.cuda.testing import skip_on_cudasim, skip_unless_conda_cudatoolkit
from numba.misc.findlib import find_lib


@skip_on_cudasim('Library detection unsupported in the simulator')
@skip_unless_conda_cudatoolkit
class TestLibraryDetection(unittest.TestCase):
    def test_detect(self):
        """
        This test is solely present to ensure that shipped cudatoolkits have
        additional core libraries in locations that Numba scans by default.
        PyCulib (and potentially others) rely on Numba's library finding
        capacity to find and subsequently load these libraries.
        """
        core_libs = ['nvvm']
        for l in core_libs:
            self.assertNotEqual(find_lib(l), [])


if __name__ == '__main__':
    unittest.main()
