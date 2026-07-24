import re

from numba import njit
from numba.core.extending import overload
from numba.core.targetconfig import ConfigStack
from numba.core.compiler import Flags, DEFAULT_FLAGS
from numba.core import types
from numba.core.funcdesc import default_mangler

from numba.tests.support import TestCase, unittest


class TestCompilerFlags(TestCase):
    def test_setting_invalid_attribute(self):
        flags = Flags()
        msg = "'Flags' object has no attribute 'this_really_does_not_exist'"
        with self.assertRaisesRegex(AttributeError, msg):
            flags.this_really_does_not_exist = True


class TestCompilerFlagCachedOverload(TestCase):
    def test_fastmath_in_overload(self):
        def fastmath_status():
            pass

        @overload(fastmath_status)
        def ov_fastmath_status():
            flags = ConfigStack().top()
            val = "Has fastmath" if flags.fastmath else "No fastmath"

            def codegen():
                return val

            return codegen

        @njit(fastmath=True)
        def set_fastmath():
            return fastmath_status()

        @njit()
        def foo():
            a = fastmath_status()
            b = set_fastmath()
            return (a, b)

        a, b = foo()
        self.assertEqual(a, "No fastmath")
        self.assertEqual(b, "Has fastmath")


class TestFlagMangling(TestCase):

    def test_demangle(self):

        def check(flags):
            mangled = flags.get_mangle_string()
            out = flags.demangle(mangled)
            # Demangle result MUST match summary()
            self.assertEqual(out, flags.summary())

        # test empty flags
        flags = Flags()
        check(flags)

        # test default
        check(DEFAULT_FLAGS)

        # test other
        flags = Flags()
        flags.no_cpython_wrapper = True
        flags.nrt = True
        flags.fastmath = True
        check(flags)

    def test_mangled_flags_is_shorter(self):
        # at least for these control cases
        flags = Flags()
        flags.nrt = True
        flags.auto_parallel = True
        self.assertLess(len(flags.get_mangle_string()), len(flags.summary()))

    def test_mangled_flags_with_fastmath_parfors_inline(self):
        # at least for these control cases
        flags = Flags()
        flags.nrt = True
        flags.auto_parallel = True
        flags.fastmath = True
        flags.inline = "always"
        self.assertLess(len(flags.get_mangle_string()), len(flags.summary()))
        demangled = flags.demangle(flags.get_mangle_string())
        # There should be no pointer value in the demangled string.
        self.assertNotIn("0x", demangled)

    def test_demangling_from_mangled_symbols(self):
        """Test demangling of flags from mangled symbol"""
        # Use default mangler to mangle the string
        fname = 'foo'
        argtypes = types.int32,
        flags = Flags()
        flags.nrt = True
        flags.inline = "always"
        name = default_mangler(
            fname, argtypes, abi_tags=[flags.get_mangle_string()],
        )
        # Find the ABI-tag. Starts with "B"
        prefix = "_Z3fooB"
        # Find the length of the ABI-tag
        m = re.match("[0-9]+", name[len(prefix):])
        size = m.group(0)
        # Extract the ABI tag
        base = len(prefix) + len(size)
        abi_mangled = name[base:base + int(size)]
        # Demangle and check
        demangled = Flags.demangle(abi_mangled)
        self.assertEqual(demangled, flags.summary())


if __name__ == "__main__":
    unittest.main()
