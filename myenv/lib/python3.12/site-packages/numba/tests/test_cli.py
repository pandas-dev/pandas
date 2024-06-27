# Tests for the CLI

import os
import subprocess
import sys
import threading
import json
from subprocess import CompletedProcess
from tempfile import TemporaryDirectory
from unittest import mock

import unittest
from numba.tests.support import TestCase, linux_only
import numba.misc.numba_sysinfo as nsi
from numba.tests.gdb_support import needs_gdb
from numba.misc.numba_gdbinfo import collect_gdbinfo
# Going to mock parts of this in testing
from numba.misc.numba_gdbinfo import _GDBTestWrapper


def run_cmd(cmdline, env=os.environ, timeout=60):
    popen = subprocess.Popen(cmdline,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             env=env)

    timeout_timer = threading.Timer(timeout, popen.kill)
    try:
        timeout_timer.start()
        out, err = popen.communicate()
        if popen.returncode != 0:
            raise AssertionError(
                "process failed with code %s: stderr follows\n%s\n" %
                (popen.returncode, err.decode()))
        return out.decode(), err.decode()
    finally:
        timeout_timer.cancel()
    return None, None


class TestCLI(TestCase):

    def test_as_module_exit_code(self):
        cmdline = [sys.executable, "-m", "numba"]
        with self.assertRaises(AssertionError) as raises:
            run_cmd(cmdline)

        self.assertIn("process failed with code 1", str(raises.exception))

    def test_sysinfo_from_module(self):
        cmdline = [sys.executable, "-m", "numba", "-s"]
        o, _ = run_cmd(cmdline)
        self.assertIn("System info", o)

    def test_json_sysinfo_from_module(self):
        with TemporaryDirectory() as d:
            path = os.path.join(d, "test_json_sysinfo.json")
            cmdline = [sys.executable, "-m", "numba", "--sys-json", path]
            run_cmd(cmdline)
            with self.subTest(msg=f"{path} exists"):
                self.assertTrue(os.path.exists(path))
            with self.subTest(msg="json load"):
                with open(path, 'r') as f:
                    info = json.load(f)
            safe_contents = {
                int: (
                    nsi._cpu_count,
                ),
                float: (
                    nsi._runtime,
                ),
                str: (
                    nsi._start,
                    nsi._start_utc,
                    nsi._machine,
                    nsi._cpu_name,
                    nsi._platform_name,
                    nsi._os_name,
                    nsi._os_version,
                    nsi._python_comp,
                    nsi._python_impl,
                    nsi._python_version,
                    nsi._llvm_version,
                ),
                bool: (
                    nsi._cu_dev_init,
                    nsi._svml_state,
                    nsi._svml_loaded,
                    nsi._svml_operational,
                    nsi._llvm_svml_patched,
                    nsi._tbb_thread,
                    nsi._openmp_thread,
                    nsi._wkq_thread,
                ),
                list: (
                    nsi._errors,
                    nsi._warnings,
                ),
                dict: (
                    nsi._numba_env_vars,
                ),
            }
            for t, keys in safe_contents.items():
                for k in keys:
                    with self.subTest(k=k):
                        self.assertIsInstance(info[k], t)

    @needs_gdb
    def test_gdb_status_from_module(self):
        # Check that the `python -m numba -g` works ok
        cmdline = [sys.executable, "-m", "numba", "-g"]
        o, _ = run_cmd(cmdline)
        self.assertIn("GDB info", o)
        # It's not known a priori whether the extension is supported, this just
        # checks that the last logical item in the output is printed.
        self.assertIn("Numba printing extension support", o)


class TestGDBCLIInfo(TestCase):

    def setUp(self):
        # Mock the entire class, to report valid things,
        # then override bits of it locally to check failures etc.

        self._patches = []

        mock_init = lambda self: None
        self._patches.append(mock.patch.object(_GDBTestWrapper, '__init__',
                                               mock_init))

        bpath = 'numba.misc.numba_gdbinfo._GDBTestWrapper.gdb_binary'
        self._patches.append(mock.patch(bpath, 'PATH_TO_GDB'))

        def _patch(fnstr, func):
            self._patches.append(mock.patch.object(_GDBTestWrapper, fnstr,
                                                   func))

        def mock_check_launch(self):
            return CompletedProcess('COMMAND STRING', 0)

        _patch('check_launch', mock_check_launch)

        # NOTE: The Python and NumPy versions are set to something unsupported!
        def mock_check_python(self):
            return CompletedProcess('COMMAND STRING', 0,
                                    stdout='(3, 2)',
                                    stderr='')

        _patch('check_python', mock_check_python)

        def mock_check_numpy(self):
            return CompletedProcess('COMMAND STRING', 0, stdout='True',
                                    stderr='')

        _patch('check_numpy', mock_check_numpy)

        def mock_check_numpy_version(self):
            return CompletedProcess('COMMAND STRING', 0, stdout='1.15',
                                    stderr='')

        _patch('check_numpy_version', mock_check_numpy_version)

        # start the patching
        for p in self._patches:
            p.start()

    def tearDown(self):
        # stop the patching
        for p in self._patches:
            p.stop()

    def test_valid(self):
        collected = collect_gdbinfo()
        self.assertEqual(collected.binary_loc, 'PATH_TO_GDB')
        extp = os.path.exists(os.path.abspath(collected.extension_loc))
        self.assertTrue(extp)
        self.assertEqual(collected.py_ver, '3.2')
        self.assertEqual(collected.np_ver, '1.15')
        self.assertIn('Full', collected.supported)

    def test_invalid_binary(self):

        def mock_fn(self):
            return CompletedProcess('INVALID_BINARY', 1)

        with mock.patch.object(_GDBTestWrapper, 'check_launch', mock_fn):
            info = collect_gdbinfo()
            self.assertIn("Testing gdb binary failed.", info.binary_loc)
            self.assertIn("gdb at 'PATH_TO_GDB' does not appear to work",
                          info.binary_loc)

    def test_no_python(self):
        def mock_fn(self):
            return CompletedProcess('NO PYTHON', 1)

        with mock.patch.object(_GDBTestWrapper, 'check_python', mock_fn):
            collected = collect_gdbinfo()
            self.assertEqual(collected.py_ver, 'No Python support')
            self.assertEqual(collected.supported, 'None')

    def test_unparsable_python_version(self):
        def mock_fn(self):
            return CompletedProcess('NO PYTHON', 0, stdout='(NOT A VERSION)')

        with mock.patch.object(_GDBTestWrapper, 'check_python', mock_fn):
            collected = collect_gdbinfo()
            self.assertEqual(collected.py_ver, 'No Python support')

    def test_no_numpy(self):
        def mock_fn(self):
            return CompletedProcess('NO NUMPY', 1)

        with mock.patch.object(_GDBTestWrapper, 'check_numpy', mock_fn):
            collected = collect_gdbinfo()
            self.assertEqual(collected.np_ver, 'No NumPy support')
            self.assertEqual(collected.py_ver, '3.2')
            self.assertIn('Partial', collected.supported)

    def test_no_numpy_version(self):
        def mock_fn(self):
            return CompletedProcess('NO NUMPY VERSION', 1)

        with mock.patch.object(_GDBTestWrapper, 'check_numpy_version', mock_fn):
            collected = collect_gdbinfo()
            self.assertEqual(collected.np_ver, 'Unknown')

    def test_traceback_in_numpy_version(self):
        def mock_fn(self):
            return CompletedProcess('NO NUMPY VERSION', 0,
                                    stdout='(NOT A VERSION)',
                                    stderr='Traceback')

        with mock.patch.object(_GDBTestWrapper, 'check_numpy_version', mock_fn):
            collected = collect_gdbinfo()
            self.assertEqual(collected.np_ver, 'Unknown')


@linux_only
class TestGDBCLIInfoBrokenGdbs(TestCase):
    # In these tests, it's necessary to actually run the entry point in an env
    # with the env var set as it's in part testing the processing of the
    # handling of consumption of variables by numba.core.config.

    def test_cannot_find_gdb_from_name(self):
        # Tests that supplying an invalid gdb binary name is handled ok
        env = os.environ.copy()
        env['NUMBA_GDB_BINARY'] = 'THIS_IS_NOT_A_VALID_GDB_BINARY_NAME'
        cmdline = [sys.executable, "-m", "numba", "-g"]
        stdout, stderr = run_cmd(cmdline, env=env)
        self.assertIn("Testing gdb binary failed", stdout)
        self.assertIn("No such file or directory", stdout)
        self.assertIn("'THIS_IS_NOT_A_VALID_GDB_BINARY_NAME'", stdout)

    def test_cannot_find_gdb_from_path(self):
        # Tests that an invalid path is handled ok
        env = os.environ.copy()
        with TemporaryDirectory() as d:
            # This wont exist as a path...
            path = os.path.join(d, 'CANNOT_EXIST')
            env['NUMBA_GDB_BINARY'] = path
            cmdline = [sys.executable, "-m", "numba", "-g"]
            stdout, stderr = run_cmd(cmdline, env=env)
            self.assertIn("Testing gdb binary failed", stdout)
            self.assertIn("No such file or directory", stdout)
            self.assertIn(path, stdout)

    def test_nonsense_gdb_binary(self):
        # Tests that a nonsense binary specified as gdb it picked up ok
        env = os.environ.copy()
        env['NUMBA_GDB_BINARY'] = 'python' # 'python' isn't gdb!
        cmdline = [sys.executable, "-m", "numba", "-g"]
        stdout, stderr = run_cmd(cmdline, env=env)
        self.assertIn("Testing gdb binary failed", stdout)
        # NOTE: should 'python' ever add support for the same flags as the gdb
        # commands used in the information gathering code in `numba_gdbinfo`
        # this test will fail, it's reasonably unlikely.
        self.assertIn("Unknown option", stdout)


if __name__ == '__main__':
    unittest.main()
