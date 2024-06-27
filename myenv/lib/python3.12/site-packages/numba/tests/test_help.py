import sys
import subprocess
import types as pytypes
import os.path

import numpy as np

import builtins
from numba.core import types
from numba.tests.support import TestCase, temp_directory
from numba.misc.help.inspector import inspect_function, inspect_module


class TestInspector(TestCase):
    def check_function_descriptor(self, info, must_be_defined=False):
        self.assertIsInstance(info, dict)
        self.assertIn('numba_type', info)
        numba_type = info['numba_type']
        if numba_type is None:
            self.assertFalse(must_be_defined)
        else:
            self.assertIsInstance(numba_type, types.Type)
            self.assertIn('explained', info)
            self.assertIsInstance(info['explained'], str)
            self.assertIn('source_infos', info)
            self.assertIsInstance(info['source_infos'], dict)

    def test_inspect_function_on_range(self):
        info = inspect_function(range)
        self.check_function_descriptor(info, must_be_defined=True)

    def test_inspect_function_on_np_all(self):
        info = inspect_function(np.all)
        self.check_function_descriptor(info, must_be_defined=True)
        source_infos = info['source_infos']
        self.assertGreater(len(source_infos), 0)
        c = 0
        for srcinfo in source_infos.values():
            self.assertIsInstance(srcinfo['kind'], str)
            self.assertIsInstance(srcinfo['name'], str)
            self.assertIsInstance(srcinfo['sig'], str)
            self.assertIsInstance(srcinfo['filename'], str)
            self.assertIsInstance(srcinfo['lines'], tuple)
            self.assertIn('docstring', srcinfo)
            c += 1
        self.assertEqual(c, len(source_infos))

    def test_inspect_module(self):
        c = 0
        for it in inspect_module(builtins):
            self.assertIsInstance(it['module'], pytypes.ModuleType)
            self.assertIsInstance(it['name'], str)
            self.assertTrue(callable(it['obj']))
            self.check_function_descriptor(it)
            c += 1
        self.assertGreater(c, 0)

    def test_inspect_cli(self):
        # Try CLI on math module
        cmdbase = [sys.executable, '-m', 'numba.misc.help.inspector']

        dirpath = temp_directory('{}.{}'.format(__name__,
                                                self.__class__.__name__))
        filename = os.path.join(dirpath, 'out')

        # Try default format "html"
        expected_file = filename + '.html'
        cmds = cmdbase + ['--file', filename, 'math']
        # File shouldn't exist yet
        self.assertFalse(os.path.isfile(expected_file))
        # Run CLI
        subprocess.check_output(cmds)
        # File should exist now
        self.assertTrue(os.path.isfile(expected_file))

        # Try changing the format to "rst"
        cmds = cmdbase + ['--file', filename, '--format', 'rst', 'math']
        expected_file = filename + '.rst'
        # File shouldn't exist yet
        self.assertFalse(os.path.isfile(expected_file))
        # Run CLI
        subprocess.check_output(cmds)
        # File should exist now
        self.assertTrue(os.path.isfile(expected_file))

        # Try unsupported format
        cmds = cmdbase + ['--file', filename, '--format', 'foo', 'math']
        # Run CLI
        with self.assertRaises(subprocess.CalledProcessError) as raises:
            subprocess.check_output(cmds, stderr=subprocess.STDOUT)
        self.assertIn("\'foo\' is not supported",
                      raises.exception.stdout.decode())
