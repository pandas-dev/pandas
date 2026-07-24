# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout


class DocsTypedSetUsageTest(unittest.TestCase):

    def test_ex_typed_set_from_cpython(self):
        with captured_stdout():
            # magictoken.ex_typed_set_from_cpython.begin
            import numpy as np
            from numba import njit
            from numba.core import types
            from numba.typed import Set

            # The Set.empty() constructs a typed set.
            # The key typed must be explicitly declared.
            s = Set.empty(types.unicode_type)

            # The typed-set can be used from the interpreter.
            s.add("One")
            s.add("Two")
            s.add("Three")

            # Here's a function that expects a typed-set as the argument
            @njit
            def remove_elem(s, elem):
                s.discard(elem)

            print(len(s))  # Out: 3

            # Call remove_elem(s) to remove an element in the typed-set.
            remove_elem(s, "Three")

            print(len(s))  # Out: 2
            # magictoken.ex_typed_set_from_cpython.end

        # Test
        np.testing.assert_equal(len(s), 2)


if __name__ == '__main__':
    unittest.main()
