from itertools import product, permutations
from collections import defaultdict

import unittest
from numba.core.base import OverloadSelector
from numba.core.registry import cpu_target
from numba.core.imputils import builtin_registry, RegistryLoader
from numba.core import types
from numba.core.errors import NumbaNotImplementedError, NumbaTypeError


class TestOverloadSelector(unittest.TestCase):
    def test_select_and_sort_1(self):
        os = OverloadSelector()
        os.append(1, (types.Any, types.Boolean))
        os.append(2, (types.Boolean, types.Integer))
        os.append(3, (types.Boolean, types.Any))
        os.append(4, (types.Boolean, types.Boolean))
        compats = os._select_compatible((types.boolean, types.boolean))
        self.assertEqual(len(compats), 3)
        ordered, scoring = os._sort_signatures(compats)
        self.assertEqual(len(ordered), 3)
        self.assertEqual(len(scoring), 3)
        self.assertEqual(ordered[0], (types.Boolean, types.Boolean))
        self.assertEqual(scoring[types.Boolean, types.Boolean], 0)
        self.assertEqual(scoring[types.Boolean, types.Any], 1)
        self.assertEqual(scoring[types.Any, types.Boolean], 1)

    def test_select_and_sort_2(self):
        os = OverloadSelector()
        os.append(1, (types.Container,))
        os.append(2, (types.Sequence,))
        os.append(3, (types.MutableSequence,))
        os.append(4, (types.List,))
        compats = os._select_compatible((types.List,))
        self.assertEqual(len(compats), 4)
        ordered, scoring = os._sort_signatures(compats)
        self.assertEqual(len(ordered), 4)
        self.assertEqual(len(scoring), 4)
        self.assertEqual(ordered[0], (types.List,))
        self.assertEqual(scoring[(types.List,)], 0)
        self.assertEqual(scoring[(types.MutableSequence,)], 1)
        self.assertEqual(scoring[(types.Sequence,)], 2)
        self.assertEqual(scoring[(types.Container,)], 3)

    def test_match(self):
        os = OverloadSelector()
        self.assertTrue(os._match(formal=types.Boolean, actual=types.boolean))
        self.assertTrue(os._match(formal=types.Boolean, actual=types.Boolean))
        # test subclass
        self.assertTrue(issubclass(types.Sequence, types.Container))
        self.assertTrue(os._match(formal=types.Container,
                                  actual=types.Sequence))
        self.assertFalse(os._match(formal=types.Sequence,
                                   actual=types.Container))
        # test any
        self.assertTrue(os._match(formal=types.Any, actual=types.Any))
        self.assertTrue(os._match(formal=types.Any, actual=types.Container))
        self.assertFalse(os._match(formal=types.Container, actual=types.Any))

    def test_ambiguous_detection(self):
        os = OverloadSelector()
        # unambiguous signatures
        os.append(1, (types.Any, types.Boolean))
        os.append(2, (types.Integer, types.Boolean))
        self.assertEqual(os.find((types.boolean, types.boolean)), 1)
        # not implemented
        with self.assertRaises(NumbaNotImplementedError) as raises:
            os.find((types.boolean, types.int32))
        # generic
        os.append(3, (types.Any, types.Any))
        self.assertEqual(os.find((types.boolean, types.int32)), 3)
        self.assertEqual(os.find((types.boolean, types.boolean)), 1)
        # add ambiguous signature; can match (bool, any) and (any, bool)
        os.append(4, (types.Boolean, types.Any))
        with self.assertRaises(NumbaTypeError) as raises:
            os.find((types.boolean, types.boolean))
        self.assertIn('2 ambiguous signatures', str(raises.exception))
        # disambiguous
        os.append(5, (types.boolean, types.boolean))
        self.assertEqual(os.find((types.boolean, types.boolean)), 5)

    def test_subclass_specialization(self):
        os = OverloadSelector()
        self.assertTrue(issubclass(types.Sequence, types.Container))
        os.append(1, (types.Container, types.Container,))
        lstty = types.List(types.boolean)
        self.assertEqual(os.find((lstty, lstty)), 1)
        os.append(2, (types.Container, types.Sequence,))
        self.assertEqual(os.find((lstty, lstty)), 2)

    def test_cache(self):
        os = OverloadSelector()
        self.assertEqual(len(os._cache), 0)
        os.append(1, (types.Any,))
        self.assertEqual(os.find((types.int32,)), 1)
        self.assertEqual(len(os._cache), 1)
        os.append(2, (types.Integer,))
        self.assertEqual(len(os._cache), 0)
        self.assertEqual(os.find((types.int32,)), 2)
        self.assertEqual(len(os._cache), 1)


class TestAmbiguousOverloads(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # ensure all impls are loaded
        cpu_target.target_context.refresh()

    def create_overload_selector(self, kind):
        os = OverloadSelector()
        loader = RegistryLoader(builtin_registry)
        for impl, sig in loader.new_registrations(kind):
            os.append(impl, sig)
        return os

    def test_ambiguous_casts(self):
        os = self.create_overload_selector(kind='casts')
        all_types = set(t for sig, impl in os.versions for t in sig)
        # ensure there are no ambiguous cast overloads
        # note: using permutations to avoid testing cast to the same type
        for sig in permutations(all_types, r=2):
            try:
                os.find(sig)
            except NumbaNotImplementedError:
                pass   # ignore not implemented cast

    def test_ambiguous_functions(self):
        loader = RegistryLoader(builtin_registry)
        selectors = defaultdict(OverloadSelector)
        for impl, fn, sig in loader.new_registrations('functions'):
            os = selectors[fn]
            os.append(impl, sig)

        for fn, os in selectors.items():
            all_types = set(t for sig, impl in os.versions for t in sig)
            # ensure there are no ambiguous overloads
            for sig in product(all_types, all_types):
                try:
                    os.find(sig)
                except NumbaNotImplementedError:
                    pass   # ignore not implemented cast



if __name__ == '__main__':
    unittest.main()
