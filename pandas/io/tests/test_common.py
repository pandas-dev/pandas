"""
    Tests for the pandas.io.common functionalities
"""
from pandas.compat import StringIO
import os
from os.path import isabs

import pandas.util.testing as tm

from pandas.io import common

from pandas import read_csv, concat

try:
    from pathlib import Path
except ImportError:
    pass

try:
    from py.path import local as LocalPath
except ImportError:
    pass


class TestCommonIOCapabilities(tm.TestCase):
    data1 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

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

    def test_stringify_path_pathlib(self):
        tm._skip_if_no_pathlib()

        rel_path = common._stringify_path(Path('.'))
        self.assertEqual(rel_path, '.')
        redundant_path = common._stringify_path(Path('foo//bar'))
        self.assertEqual(redundant_path, os.path.join('foo', 'bar'))

    def test_stringify_path_localpath(self):
        tm._skip_if_no_localpath()

        path = os.path.join('foo', 'bar')
        abs_path = os.path.abspath(path)
        lpath = LocalPath(path)
        self.assertEqual(common._stringify_path(lpath), abs_path)

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

    def test_iterator(self):
        reader = read_csv(StringIO(self.data1), chunksize=1)
        result = concat(reader, ignore_index=True)
        expected = read_csv(StringIO(self.data1))
        tm.assert_frame_equal(result, expected)

        # GH12153
        it = read_csv(StringIO(self.data1), chunksize=1)
        first = next(it)
        tm.assert_frame_equal(first, expected.iloc[[0]])
        expected.index = [0 for i in range(len(expected))]
        tm.assert_frame_equal(concat(it), expected.iloc[1:])
