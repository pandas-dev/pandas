import unittest

from numba.tests.support import TestCase

from numba.misc.init_utils import version_info, generate_version_info


class TestGenerateVersionInfo(TestCase):

    def test_major_minor_patch(self):
        expected = version_info(0, 1, 0,
                                (0, 1), (0, 1, 0),
                                "0.1.0", ('0', '1', '0'), None)
        received = generate_version_info("0.1.0")
        self.assertEqual(received, expected)

    def test_unknown(self):
        expected = version_info(None, None, None,
                                (None, None), (None, None, None),
                                '0+unknown', ('0+unknown',), None)
        received = generate_version_info('0+unknown')
        self.assertEqual(received, expected)

    def test_dev(self):
        expected = version_info(0, 1, None,
                                (0, 1), (0, 1, None),
                                '0.1.0dev0', ('0', '1', '0dev0'), None)
        received = generate_version_info('0.1.0dev0')
        self.assertEqual(received, expected)

    def test_full_rev(self):
        expected = version_info(0, 1, None,
                                (0, 1), (0, 1, None),
                                '0.1.0dev0+1.g0123456789abcdef',
                                ('0', '1', '0dev0+1', 'g0123456789abcdef'),
                                'g0123456789abcdef')
        received = generate_version_info('0.1.0dev0+1.g0123456789abcdef')
        self.assertEqual(received, expected)


if __name__ == '__main__':
    unittest.main()
