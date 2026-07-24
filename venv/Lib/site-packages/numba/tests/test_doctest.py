import doctest
import unittest
from numba.tests.support import TestCase


class TestDocTest(TestCase):
    def test_basic_decorators(self):
        from . import doctest_usecase

        # Make sure the finder see all the doctest
        finder = doctest.DocTestFinder()
        tests = finder.find(doctest_usecase)
        testnames = {x.name for x in tests}
        expected = {
            'numba.tests.doctest_usecase',
            'numba.tests.doctest_usecase.a',
            'numba.tests.doctest_usecase.b',
            'numba.tests.doctest_usecase.c',
            'numba.tests.doctest_usecase.d',
        }
        self.assertEqual(testnames, expected)

        # Execute the doctest in the module
        doctest.testmod(doctest_usecase)


if __name__ == "__main__":
    unittest.main()
