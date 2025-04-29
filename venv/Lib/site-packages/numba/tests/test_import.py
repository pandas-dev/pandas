import unittest
from numba.tests.support import TestCase, run_in_subprocess


class TestNumbaImport(TestCase):
    """
    Test behaviour of importing Numba.
    """

    def test_laziness(self):
        """
        Importing top-level numba features should not import too many modules.
        """
        # A heuristic set of modules that shouldn't be imported immediately
        banlist = ['cffi',
                   'distutils',
                   'numba.cuda',
                   'numba.cpython.mathimpl',
                   'numba.cpython.randomimpl',
                   'numba.tests',
                   'numba.core.typing.collections',
                   'numba.core.typing.listdecl',
                   'numba.core.typing.npdatetime',
                   ]
        # Sanity check the modules still exist...
        for mod in banlist:
            if mod not in ('cffi',):
                __import__(mod)

        code = """if 1:
            from numba import jit, vectorize
            from numba.core import types
            import sys
            print(list(sys.modules))
            """

        out, _ = run_in_subprocess(code)
        modlist = set(eval(out.strip()))
        unexpected = set(banlist) & set(modlist)
        self.assertFalse(unexpected, "some modules unexpectedly imported")

    def test_no_impl_import(self):
        """
        Tests that importing jit does not trigger import of modules containing
        lowering implementations that would likely install things in the
        builtins registry and have side effects impacting other targets
        """
        # None of these modules should be imported through the process of
        # doing 'import numba' or 'from numba import njit'
        banlist = ['numba.cpython.slicing',
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
                   'numba.experimental.jitclass.base',]

        code1 = """if 1:
            import sys
            import numba
            print(list(sys.modules))
            """

        code2 = """if 1:
            import sys
            from numba import njit
            @njit
            def foo():
                pass
            print(list(sys.modules))
            """

        for code in (code1, code2):
            out, _ = run_in_subprocess(code)
            modlist = set(eval(out.strip()))
            unexpected = set(banlist) & set(modlist)
            self.assertFalse(unexpected, "some modules unexpectedly imported")

    def test_no_accidental_warnings(self):
        # checks that importing Numba isn't accidentally triggering warnings due
        # to e.g. deprecated use of import locations from Python's stdlib
        code = "import numba"
        # See: https://github.com/numba/numba/issues/6831
        # bug in setuptools/packaging causing a deprecation warning
        flags = ["-Werror", "-Wignore::DeprecationWarning:packaging.version:"]
        run_in_subprocess(code, flags)

    def test_import_star(self):
        # checks that "from numba import *" works.
        code = "from numba import *"
        run_in_subprocess(code)


if __name__ == '__main__':
    unittest.main()
