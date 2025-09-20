from numba.tests.support import run_in_subprocess
import unittest


class TestImport(unittest.TestCase):
    def test_no_impl_import(self):
        """
        Tests that importing cuda doesn't trigger the import of modules
        containing lowering implementation that would likely install things in
        the builtins registry and have side effects impacting other targets.
        """

        banlist = (
            'numba.cpython.slicing',
            'numba.cpython.tupleobj',
            'numba.cpython.enumimpl',
            'numba.cpython.hashing',
            'numba.cpython.heapq',
            'numba.cpython.iterators',
            'numba.cpython.numbers',
            'numba.cpython.rangeobj',
            'numba.cpython.cmathimpl',
            'numba.cpython.mathimpl',
            'numba.cpython.printimpl',
            'numba.cpython.randomimpl',
            'numba.core.optional',
            'numba.misc.gdb_hook',
            'numba.misc.literal',
            'numba.misc.cffiimpl',
            'numba.np.linalg',
            'numba.np.polynomial',
            'numba.np.arraymath',
            'numba.np.npdatetime',
            'numba.np.npyimpl',
            'numba.typed.typeddict',
            'numba.typed.typedlist',
            'numba.experimental.jitclass.base',
        )

        code = "import sys; from numba import cuda; print(list(sys.modules))"

        out, _ = run_in_subprocess(code)
        modlist = set(eval(out.strip()))
        unexpected = set(banlist) & set(modlist)
        self.assertFalse(unexpected, "some modules unexpectedly imported")


if __name__ == '__main__':
    unittest.main()
