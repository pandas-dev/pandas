"""
    Tests for the pandas.io.common functionalities
"""
from pandas.compat import StringIO
import os
from os.path import isabs

import pandas.util.testing as tm

from pandas.io import common


class TestCommonIOCapabilities(tm.TestCase):

    def test_expand_user(self):
        filename = '~/sometest'
        expanded_name = common._expand_user(filename)

        self.assertNotEqual(expanded_name, filename)
        self.assertTrue(isabs(expanded_name))
        self.assertEqual(os.path.expanduser(filename), expanded_name)

    def test_expand_user_normal_path(self):
        filename = '/somefolder/sometest'
        expanded_name = common._expand_user(filename)

        self.assertEqual(expanded_name, filename)
        self.assertEqual(os.path.expanduser(filename), expanded_name)

    def test_get_filepath_or_buffer_with_path(self):
        filename = '~/sometest'
        filepath_or_buffer, _, _ = common.get_filepath_or_buffer(filename)
        self.assertNotEqual(filepath_or_buffer, filename)
        self.assertTrue(isabs(filepath_or_buffer))
        self.assertEqual(os.path.expanduser(filename), filepath_or_buffer)

    def test_get_filepath_or_buffer_with_buffer(self):
        input_buffer = StringIO()
        filepath_or_buffer, _, _ = common.get_filepath_or_buffer(input_buffer)
        self.assertEqual(filepath_or_buffer, input_buffer)
