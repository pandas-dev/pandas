"""Tests for moved modules and their redirection from old path
"""
from numba.tests.support import TestCase


class TestMovedModule(TestCase):
    """Testing moved modules in Q1 2020 but were decided to kept as public API
    """
    def tests_numba_types(self):
        import numba.types
        import numba.core.types as types
        # The old module IS NOT the new module
        self.assertIsNot(numba.types, types)
        # Attribute access are there
        self.assertIs(numba.types.intp, types.intp)
        self.assertIs(numba.types.float64, types.float64)
        self.assertIs(numba.types.Array, types.Array)
        # Submodule access through old import path is possible
        import numba.types.misc
        self.assertIs(types.misc, numba.types.misc)
        self.assertIs(types.misc.Optional, numba.types.misc.Optional)
        # Import time code could be executed twice and causes the following to
        # fail.
        self.assertIs(types.StringLiteral, numba.types.misc.StringLiteral)
        # Check numba.types.container
        from numba.types import containers
        self.assertIs(types.containers, containers)
        self.assertIs(types.containers.Sequence, containers.Sequence)
        from numba.types.containers import Sequence
        self.assertIs(Sequence, containers.Sequence)
