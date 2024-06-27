
import unittest

from numba import njit
from numba.core.funcdesc import PythonFunctionDescriptor, default_mangler
from numba.core.compiler import run_frontend
from numba.core.itanium_mangler import mangle_abi_tag


class TestModule(unittest.TestCase):
    def test_module_not_in_namespace(self):
        """ Test of trying to run a compiled function
        where the module from which the function is being compiled
        doesn't exist in the namespace.
        """
        filename = 'test.py'
        name = 'mypackage'
        code = """
def f(x):
    return x
"""

        objs = dict(__file__=filename, __name__=name)
        compiled = compile(code, filename, 'exec')
        exec(compiled, objs)

        compiled_f = njit(objs['f'])
        self.assertEqual(compiled_f(3), 3)


class TestFuncDescMangledName(unittest.TestCase):
    def test_mangling_abi_tags(self):
        """
        This is a minimal test for the abi-tags support in the mangler.
        """
        def udt():
            pass

        # run minimal frontend to create a function descriptor
        func_ir = run_frontend(udt)
        typemap = {}
        restype = None
        calltypes = ()
        mangler = default_mangler
        inline = False
        noalias = False

        abi_tags = ("Shrubbery", "Herring")

        fd = PythonFunctionDescriptor.from_specialized_function(
            func_ir, typemap, restype, calltypes, mangler, inline, noalias,
            abi_tags=abi_tags,
        )

        # mangled tag must exist in the mangled name
        self.assertIn("".join([mangle_abi_tag(x) for x in abi_tags]),
                      fd.mangled_name)


if __name__ == '__main__':
    unittest.main()
