from numba.tests.support import TestCase, unittest
from numba.misc import findlib


class TestFindlib(TestCase):
    def test_find_file_nonexistent_path(self):
        candidates = findlib.find_file('libirrelevant.so', 'NONEXISTENT')
        self.assertEqual(candidates, [])


if __name__ == '__main__':
    unittest.main()
