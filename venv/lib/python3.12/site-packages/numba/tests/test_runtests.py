import os
import sys
import subprocess

from numba import cuda
import unittest
import itertools

try:
    import git  # noqa: F401 from gitpython package
except ImportError:
    has_gitpython = False
else:
    has_gitpython = True

try:
    import yaml  # from pyyaml package
except ImportError:
    has_pyyaml = False
else:
    has_pyyaml = True


class TestCase(unittest.TestCase):
    """These test cases are meant to test the Numba test infrastructure itself.
    Therefore, the logic used here shouldn't use numba.testing, but only the
    upstream unittest, and run the numba test suite only in a subprocess."""

    def get_testsuite_listing(self, args, *, subp_kwargs=None):
        """
        Use `subp_kwargs` to pass extra argument to `subprocess.check_output`.
        """
        subp_kwargs = subp_kwargs or {}
        cmd = [sys.executable, '-m', 'numba.runtests', '-l'] + list(args)
        out_bytes = subprocess.check_output(cmd, **subp_kwargs)
        lines = out_bytes.decode('UTF-8').splitlines()
        lines = [line for line in lines if line.strip()]
        return lines

    def check_listing_prefix(self, prefix):
        listing = self.get_testsuite_listing([prefix])
        for ln in listing[:-1]:
            errmsg = '{!r} not startswith {!r}'.format(ln, prefix)
            self.assertTrue(ln.startswith(prefix), msg=errmsg)

    def check_testsuite_size(self, args, minsize):
        """
        Check that the reported numbers of tests are at least *minsize*.
        """
        lines = self.get_testsuite_listing(args)
        last_line = lines[-1]
        self.assertTrue('tests found' in last_line)
        number = int(last_line.split(' ')[0])
        # There may be some "skipped" messages at the beginning,
        # so do an approximate check.
        self.assertIn(len(lines), range(number + 1, number + 20))
        self.assertGreaterEqual(number, minsize)
        return lines

    def check_all(self, ids):
        lines = self.check_testsuite_size(ids, 5000)
        # CUDA should be included by default
        self.assertTrue(any('numba.cuda.tests.' in line for line in lines))
        # As well as subpackage
        self.assertTrue(any('numba.tests.npyufunc.test_' in line
                            for line in lines),)

    def _get_numba_tests_from_listing(self, listing):
        """returns a filter on strings starting with 'numba.', useful for
        selecting the 'numba' test names from a test listing."""
        return filter(lambda x: x.startswith('numba.'), listing)

    def test_default(self):
        self.check_all([])

    def test_all(self):
        self.check_all(['numba.tests'])

    def test_cuda(self):
        # Even without CUDA enabled, there is at least one test
        # (in numba.cuda.tests.nocuda)
        minsize = 100 if cuda.is_available() else 1
        self.check_testsuite_size(['numba.cuda.tests'], minsize)

    @unittest.skipIf(not cuda.is_available(), "NO CUDA")
    def test_cuda_submodules(self):
        self.check_listing_prefix('numba.cuda.tests.cudadrv')
        self.check_listing_prefix('numba.cuda.tests.cudapy')
        self.check_listing_prefix('numba.cuda.tests.nocuda')
        self.check_listing_prefix('numba.cuda.tests.cudasim')

    def test_module(self):
        self.check_testsuite_size(['numba.tests.test_storeslice'], 2)
        self.check_testsuite_size(['numba.tests.test_nested_calls'], 10)
        # Several modules
        self.check_testsuite_size(['numba.tests.test_nested_calls',
                                   'numba.tests.test_storeslice'], 12)

    def test_subpackage(self):
        self.check_testsuite_size(['numba.tests.npyufunc'], 50)

    def test_random(self):
        self.check_testsuite_size(
            ['--random', '0.1', 'numba.tests.npyufunc'], 5)

    def test_include_exclude_tags(self):
        def get_count(arg_list):
            lines = self.get_testsuite_listing(arg_list)
            self.assertIn('tests found', lines[-1])
            count = int(lines[-1].split()[0])
            self.assertTrue(count > 0)
            return count

        tags = ['long_running', 'long_running, important']

        total = get_count(['numba.tests'])

        for tag in tags:
            included = get_count(['--tags', tag, 'numba.tests'])
            excluded = get_count(['--exclude-tags', tag, 'numba.tests'])
            self.assertEqual(total, included + excluded)

            # check syntax with `=` sign in
            included = get_count(['--tags=%s' % tag, 'numba.tests'])
            excluded = get_count(['--exclude-tags=%s' % tag, 'numba.tests'])
            self.assertEqual(total, included + excluded)

    def test_check_shard(self):
        tmpAll = self.get_testsuite_listing([])
        tmp1 = self.get_testsuite_listing(['-j', '0:2'])
        tmp2 = self.get_testsuite_listing(['-j', '1:2'])

        lAll = set(self._get_numba_tests_from_listing(tmpAll))
        l1 = set(self._get_numba_tests_from_listing(tmp1))
        l2 = set(self._get_numba_tests_from_listing(tmp2))

        # The difference between two adjacent shards should be less than 5% of
        # the total
        self.assertLess(abs(len(l2) - len(l1)), len(lAll) / 20)
        self.assertLess(len(l1), len(lAll))
        self.assertLess(len(l2), len(lAll))

    def test_check_sharding_equivalent(self):
        # get some shards
        sharded = list()
        for i in range(3):
            subset = self.get_testsuite_listing(['-j', '{}:3'.format(i)])
            slist = [*self._get_numba_tests_from_listing(subset)]
            sharded.append(slist)

        # get the always running tests
        tmp = self.get_testsuite_listing(['--tag', 'always_test'])
        always_running = set(self._get_numba_tests_from_listing(tmp))

        # make sure there is at least one test that always runs
        self.assertGreaterEqual(len(always_running), 1)

        # check that each shard contains no repeats
        sharded_sets = [set(x) for x in sharded]
        for i in range(len(sharded)):
            self.assertEqual(len(sharded_sets[i]), len(sharded[i]))

        # check that the always running tests are in every shard, and then
        # remove them from the shards
        for shard in sharded_sets:
            for test in always_running:
                self.assertIn(test, shard)
                shard.remove(test)
                self.assertNotIn(test, shard)

        # check that there is no overlap between the shards
        for a, b in itertools.combinations(sharded_sets, 2):
            self.assertFalse(a & b)

        # check that the sum of the shards and the always running tests is the
        # same as the full listing

        sum_of_parts = set()
        for x in sharded_sets:
            sum_of_parts.update(x)
        sum_of_parts.update(always_running)

        full_listing = set(self._get_numba_tests_from_listing(
            self.get_testsuite_listing([])))

        self.assertEqual(sum_of_parts, full_listing)

    @unittest.skipUnless(has_gitpython, "Requires gitpython")
    def test_gitdiff(self):
        # Check for git
        try:
            subprocess.call("git",
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            self.skipTest("no git available")

        # default
        outs = self.get_testsuite_listing(['-g'])
        self.assertNotIn("Git diff by common ancestor", outs)
        # using ancestor
        outs = self.get_testsuite_listing(['-g=ancestor'])
        self.assertIn("Git diff by common ancestor", outs)
        # misspelled ancestor
        subp_kwargs = dict(stderr=subprocess.DEVNULL)
        with self.assertRaises(subprocess.CalledProcessError):
            self.get_testsuite_listing(['-g=ancest'], subp_kwargs=subp_kwargs)

    @unittest.skipUnless(has_pyyaml, "Requires pyyaml")
    def test_azure_config(self):
        from yaml import Loader
        base_path = os.path.dirname(os.path.abspath(__file__))
        azure_pipe = os.path.join(base_path, '..', '..', 'azure-pipelines.yml')
        if not os.path.isfile(azure_pipe):
            self.skipTest("'azure-pipelines.yml' is not available")
        with open(os.path.abspath(azure_pipe), 'rt') as f:
            data = f.read()
        pipe_yml = yaml.load(data, Loader=Loader)

        templates = pipe_yml['jobs']
        # first look at the items in the first two templates, this is osx/linux
        start_indexes = []
        for tmplt in templates[:2]:
            matrix = tmplt['parameters']['matrix']
            for setup in matrix.values():
                start_indexes.append(setup['TEST_START_INDEX'])

        # next look at the items in the windows only template
        winpath = ['..', '..', 'buildscripts', 'azure', 'azure-windows.yml']
        azure_windows = os.path.join(base_path, *winpath)
        if not os.path.isfile(azure_windows):
            self.skipTest("'azure-windows.yml' is not available")
        with open(os.path.abspath(azure_windows), 'rt') as f:
            data = f.read()
        windows_yml = yaml.load(data, Loader=Loader)

        # There's only one template in windows and its keyed differently to the
        # above, get its matrix.
        matrix = windows_yml['jobs'][0]['strategy']['matrix']
        for setup in matrix.values():
            start_indexes.append(setup['TEST_START_INDEX'])

        # sanity checks
        # 1. That the TEST_START_INDEX is unique
        self.assertEqual(len(start_indexes), len(set(start_indexes)))
        # 2. That the TEST_START_INDEX is a complete range
        lim_start_index = max(start_indexes) + 1
        expected = [*range(lim_start_index)]
        self.assertEqual(sorted(start_indexes), expected)
        # 3. That the number of indexes matches the declared test count
        self.assertEqual(lim_start_index, pipe_yml['variables']['TEST_COUNT'])

    def test_no_compilation_on_list(self):
        # Checks that the test suite doesn't do any CPU-side compilation on
        # listing of tests.
        code = """if 1:
        from unittest import mock
        from llvmlite import binding as llvm
        error = RuntimeError("Detected compilation during test listing")
        with mock.patch.object(llvm.ExecutionEngine, 'finalize_object',
                               side_effect=error):
            import numba
            {0}
        """

        # Run with a jit function in the test to demonstrate failure
        with self.assertRaises(subprocess.CalledProcessError) as raises:
            cmd = [sys.executable, "-c", code.format("numba.njit(lambda:0)()")]
            subprocess.check_output(cmd,
                                    stderr=subprocess.STDOUT,
                                    timeout=60)
        self.assertIn("Detected compilation during test listing",
                      raises.exception.stdout.decode('UTF-8'))

        # Run to validate the test suite does not trigger compilation during
        # listing.
        cmd = [sys.executable, "-c", code.format("numba.test('-l')")]
        subprocess.check_call(cmd,
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL,)


if __name__ == '__main__':
    unittest.main()
