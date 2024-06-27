import collections
import contextlib
import cProfile
import inspect
import gc
import multiprocessing
import os
import random
import sys
import time
import unittest
import warnings
import zlib
import re

from functools import lru_cache
from io import StringIO
from unittest import result, runner, signals, suite, loader, case
from datetime import datetime

from .loader import TestLoader

try:
    from multiprocessing import TimeoutError
except ImportError:
    from Queue import Empty as TimeoutError


def make_tag_decorator(known_tags):
    """
    Create a decorator allowing tests to be tagged with the *known_tags*.
    """

    def tag(*tags):
        """
        Tag a test method with the given tags.
        Can be used in conjunction with the --tags command-line argument
        for runtests.py.
        """
        for t in tags:
            if t not in known_tags:
                raise ValueError("unknown tag: %r" % (t,))

        def decorate(func):
            if (not callable(func) or isinstance(func, type)
                or not func.__name__.startswith('test_')):
                raise TypeError("@tag(...) should be used on test methods")
            try:
                s = func.tags
            except AttributeError:
                s = func.tags = set()
            s.update(tags)
            return func
        return decorate

    return tag


# Chances are the next queried class is the same as the previous, locally 128
# entries seems to be fastest.
# Current number of test classes can be found with:
# $ ./runtests.py -l|sed -e 's/\(.*\)\..*/\1/'|grep ^numba|sort|uniq|wc -l
# as of writing it's 658.
@lru_cache(maxsize=128)
def _get_mtime(cls):
    """
    Gets the mtime of the file in which a test class is defined.
    """
    return str(os.path.getmtime(inspect.getfile(cls)))


def cuda_sensitive_mtime(x):
    """
    Return a key for sorting tests bases on mtime and test name. For CUDA
    tests, interleaving tests from different classes is dangerous as the CUDA
    context might get reset unexpectedly between methods of a class, so for
    CUDA tests the key prioritises the test module and class ahead of the
    mtime.
    """
    cls = x.__class__
    key = _get_mtime(cls) + str(x)

    from numba.cuda.testing import CUDATestCase
    if CUDATestCase in cls.mro():
        key = "%s.%s %s" % (str(cls.__module__), str(cls.__name__), key)

    return key


def parse_slice(useslice):
    """Parses the argument string "useslice" as a shard index and number and
    returns a function that filters on those arguments. i.e. input
    useslice="1:3" leads to output something like `lambda x: zlib.crc32(x) % 3
    == 1`.
    """
    if callable(useslice):
        return useslice
    if not useslice:
        return lambda x: True
    try:
        (index, count) = useslice.split(":")
        index = int(index)
        count = int(count)
    except Exception:
        msg = (
                    "Expected arguments shard index and count to follow "
                    "option `-j i:t`, where i is the shard number and t "
                    "is the total number of shards, found '%s'" % useslice)
        raise ValueError(msg)
    if count == 0:
        return lambda x: True
    elif count < 0 or index < 0 or index >= count:
        raise ValueError("Sharding out of range")
    else:
        def decide(test):
            func = getattr(test, test._testMethodName)
            if "always_test" in getattr(func, 'tags', {}):
                return True
            return abs(zlib.crc32(test.id().encode('utf-8'))) % count == index
        return decide


class TestLister(object):
    """Simply list available tests rather than running them."""
    def __init__(self, useslice):
        self.useslice = parse_slice(useslice)

    def run(self, test):
        result = runner.TextTestResult(sys.stderr, descriptions=True,
                                       verbosity=1)
        self._test_list = _flatten_suite(test)
        masked_list = list(filter(self.useslice, self._test_list))
        self._test_list.sort(key=cuda_sensitive_mtime)
        for t in masked_list:
            print(t.id())
        print('%d tests found. %s selected' % (len(self._test_list),
                                               len(masked_list)))
        return result


class SerialSuite(unittest.TestSuite):
    """A simple marker to make sure tests in this suite are run serially.

    Note: As the suite is going through internals of unittest,
          it may get unpacked and stuffed into a plain TestSuite.
          We need to set an attribute on the TestCase objects to
          remember they should not be run in parallel.
    """

    def addTest(self, test):
        if not isinstance(test, unittest.TestCase):
            # It's a sub-suite, recurse
            for t in test:
                self.addTest(t)
        else:
            # It's a test case, mark it serial
            test._numba_parallel_test_ = False
            super(SerialSuite, self).addTest(test)


class BasicTestRunner(runner.TextTestRunner):
    def __init__(self, useslice, **kwargs):
        runner.TextTestRunner.__init__(self, **kwargs)
        self.useslice = parse_slice(useslice)

    def run(self, test):
        run = list(filter(self.useslice, _flatten_suite(test)))
        run.sort(key=cuda_sensitive_mtime)
        wrapped = unittest.TestSuite(run)
        return super(BasicTestRunner, self).run(wrapped)


# "unittest.main" is really the TestProgram class!
# (defined in a module named itself "unittest.main"...)

class NumbaTestProgram(unittest.main):
    """
    A TestProgram subclass adding the following options:
    * a -R option to enable reference leak detection
    * a --profile option to enable profiling of the test run
    * a -m option for parallel execution
    * a -l option to (only) list tests

    Currently the options are only added in 3.4+.
    """

    refleak = False
    profile = False
    multiprocess = False
    useslice = None
    list = False
    tags = None
    exclude_tags = None
    random_select = None
    random_seed = 42
    show_timing = False
    write_junit = False

    def __init__(self, *args, **kwargs):
        topleveldir = kwargs.pop('topleveldir', None)
        kwargs['testLoader'] = TestLoader(topleveldir)

        # HACK to force unittest not to change warning display options
        # (so that NumbaWarnings don't appear all over the place)
        sys.warnoptions.append(':x')
        self.nomultiproc = kwargs.pop('nomultiproc', False)
        super(NumbaTestProgram, self).__init__(*args, **kwargs)

    def _getParentArgParser(self):
        # NOTE: this hook only exists on Python 3.4+. The options won't be
        # added in earlier versions (which use optparse - 3.3 - or getopt()
        # - 2.x).
        parser = super(NumbaTestProgram, self)._getParentArgParser()
        if self.testRunner is None:
            parser.add_argument('-R', '--refleak', dest='refleak',
                                action='store_true',
                                help='Detect reference / memory leaks')
        parser.add_argument('-m', '--multiprocess', dest='multiprocess',
                            nargs='?',
                            type=int,
                            const=multiprocessing.cpu_count(),
                            help='Parallelize tests')
        parser.add_argument('-l', '--list', dest='list',
                            action='store_true',
                            help='List tests without running them')
        parser.add_argument('--tags', dest='tags', type=str,
                            help='Comma-separated list of tags to select '
                                 'a subset of the test suite')
        parser.add_argument('--exclude-tags', dest='exclude_tags', type=str,
                            help='Comma-separated list of tags to de-select '
                                 'a subset of the test suite')
        parser.add_argument('--random', dest='random_select', type=float,
                            help='Random proportion of tests to select')
        parser.add_argument('--profile', dest='profile',
                            action='store_true',
                            help='Profile the test run')
        parser.add_argument('-j', '--slice', dest='useslice', nargs='?',
                            type=str, const="None",
                            help='Shard the test sequence')
        parser.add_argument('--show-timing', dest='show_timing', action='store_true',
                            help=("Print recorded timing for each test and "
                                  "the whole suite."))
        parser.add_argument('--junit', dest='write_junit', action='store_true',
                            help=("Write junit xml output"))

        def git_diff_str(x):
            if x != 'ancestor':
                raise ValueError("invalid option for --gitdiff")
            return x

        parser.add_argument('-g', '--gitdiff', dest='gitdiff', type=git_diff_str,
                            default=False, nargs='?',
                            help=('Run tests from changes made against '
                                  'origin/release0.60 as identified by `git diff`. '
                                  'If set to "ancestor", the diff compares '
                                  'against the common ancestor.'))
        return parser

    def _handle_tags(self, argv, tagstr):
        found = None
        for x in argv:
            if tagstr in x:
                if found is None:
                    found = x
                else:
                    raise ValueError("argument %s supplied repeatedly" % tagstr)

        if found is not None:
            posn = argv.index(found)
            try:
                if found == tagstr: # --tagstr <arg>
                    tag_args = argv[posn + 1].strip()
                    argv.remove(tag_args)
                else: # --tagstr=<arg>
                    if '=' in found:
                        tag_args =  found.split('=')[1].strip()
                    else:
                        raise AssertionError('unreachable')
            except IndexError:
                # at end of arg list, raise
                msg = "%s requires at least one tag to be specified"
                raise ValueError(msg % tagstr)
            # see if next arg is "end options" or some other flag
            if tag_args.startswith('-'):
                raise ValueError("tag starts with '-', probably a syntax error")
            # see if tag is something like "=<tagname>" which is likely a syntax
            # error of form `--tags =<tagname>`, note the space prior to `=`.
            if '=' in tag_args:
                msg = "%s argument contains '=', probably a syntax error"
                raise ValueError(msg % tagstr)
            attr = tagstr[2:].replace('-', '_')
            setattr(self, attr, tag_args)
            argv.remove(found)


    def parseArgs(self, argv):
        if '-l' in argv:
            argv.remove('-l')
            self.list = True

        super(NumbaTestProgram, self).parseArgs(argv)

        # If at this point self.test doesn't exist, it is because
        # no test ID was given in argv. Use the default instead.
        if not hasattr(self, 'test') or not self.test.countTestCases():
            self.testNames = (self.defaultTest,)
            self.createTests()

        if self.tags:
            tags = [s.strip() for s in self.tags.split(',')]
            self.test = _choose_tagged_tests(self.test, tags, mode='include')

        if self.exclude_tags:
            tags = [s.strip() for s in self.exclude_tags.split(',')]
            self.test = _choose_tagged_tests(self.test, tags, mode='exclude')

        if self.random_select:
            self.test = _choose_random_tests(self.test, self.random_select,
                                             self.random_seed)

        if self.gitdiff is not False:
            self.test = _choose_gitdiff_tests(
                self.test,
                use_common_ancestor=(self.gitdiff == 'ancestor'),
            )

        if self.verbosity <= 0:
            # We aren't interested in informational messages / warnings when
            # running with '-q'.
            self.buffer = True

    def _do_discovery(self, argv, Loader=None):
        # Disable unittest's implicit test discovery when parsing
        # CLI arguments, as it can select other tests than Numba's
        # (e.g. some test_xxx module that may happen to be directly
        #  reachable from sys.path)
        return

    def runTests(self):
        if self.refleak:
            self.testRunner = RefleakTestRunner

            if not hasattr(sys, "gettotalrefcount"):
                warnings.warn("detecting reference leaks requires a debug build "
                              "of Python, only memory leaks will be detected")

        elif self.list:
            self.testRunner = TestLister(self.useslice)

        elif self.testRunner is None:
            self.testRunner = BasicTestRunner(self.useslice,
                                              verbosity=self.verbosity,
                                              failfast=self.failfast,
                                              buffer=self.buffer)

        if self.multiprocess and not self.nomultiproc:
            if self.multiprocess < 1:
                msg = ("Value specified for the number of processes to use in "
                    "running the suite must be > 0")
                raise ValueError(msg)
            self.testRunner = ParallelTestRunner(runner.TextTestRunner,
                                                 self.multiprocess,
                                                 self.useslice,
                                                 show_timing=self.show_timing,
                                                 write_junit=self.write_junit,
                                                 verbosity=self.verbosity,
                                                 failfast=self.failfast,
                                                 buffer=self.buffer)

        def run_tests_real():
            super(NumbaTestProgram, self).runTests()

        if self.profile:
            filename = os.path.splitext(
                os.path.basename(sys.modules['__main__'].__file__)
                )[0] + '.prof'
            p = cProfile.Profile(timer=time.perf_counter)  # 3.3+
            p.enable()
            try:
                p.runcall(run_tests_real)
            finally:
                p.disable()
                print("Writing test profile data into %r" % (filename,))
                p.dump_stats(filename)
        else:
            run_tests_real()


# These are tests which are generated and injected into the test suite, what
# gets injected depends on features of the test environment, e.g. TBB presence
# it's important for doing the CI "slice tests" that these are run at the end
# See notes in `_flatten_suite` for why. Simple substring matching is used to
# determine a match.
_GENERATED = (
    "numba.cuda.tests.cudapy.test_libdevice.TestLibdeviceCompilation",
    "numba.tests.test_num_threads",
    "numba.tests.test_parallel_backend",
    "numba.tests.test_svml",
    "numba.tests.test_ufuncs",
)


def _flatten_suite_inner(test):
    """
    Workhorse for _flatten_suite
    """
    tests = []
    if isinstance(test, (unittest.TestSuite, list, tuple)):
        for x in test:
            tests.extend(_flatten_suite_inner(x))
    else:
        tests.append(test)
    return tests


def _flatten_suite(test):
    """
    Expand nested suite into list of test cases.
    """
    tests = _flatten_suite_inner(test)
    # Strip out generated tests and stick them at the end, this is to make sure
    # that tests appear in a consistent order regardless of features available.
    # This is so that a slice through the test suite e.g. (1::N) would likely be
    # consistent up to the point of the generated tests, which rely on specific
    # features.
    generated = set()
    for t in tests:
        for g in _GENERATED:
            if g in str(t):
                generated.add(t)
    normal = set(tests) - generated
    def key(x):
        return x.__module__, type(x).__name__, x._testMethodName
    tests = sorted(normal, key=key)
    tests.extend(sorted(list(generated), key=key))
    return tests


def _choose_gitdiff_tests(tests, *, use_common_ancestor=False):
    try:
        from git import Repo
    except ImportError:
        raise ValueError("gitpython needed for git functionality")
    repo = Repo('.')
    path = os.path.join('numba', 'tests')
    if use_common_ancestor:
        print(f"Git diff by common ancestor")
        target = 'origin/release0.60...HEAD'
    else:
        target = 'origin/release0.60..HEAD'
    gdiff_paths = repo.git.diff(target, path, name_only=True).split()
    # normalise the paths as they are unix style from repo.git.diff
    gdiff_paths = [os.path.normpath(x) for x in gdiff_paths]
    selected = []
    gdiff_paths = [os.path.join(repo.working_dir, x) for x in gdiff_paths]
    for test in _flatten_suite(tests):
        assert isinstance(test, unittest.TestCase)
        fname = inspect.getsourcefile(test.__class__)
        if fname in gdiff_paths:
            selected.append(test)
    print("Git diff identified %s tests" % len(selected))
    return unittest.TestSuite(selected)

def _choose_tagged_tests(tests, tags, mode='include'):
    """
    Select tests that are tagged/not tagged with at least one of the given tags.
    Set mode to 'include' to include the tests with tags, or 'exclude' to
    exclude the tests with the tags.
    """
    selected = []
    tags = set(tags)
    for test in _flatten_suite(tests):
        assert isinstance(test, unittest.TestCase)
        func = getattr(test, test._testMethodName)
        try:
            # Look up the method's underlying function (Python 2)
            func = func.im_func
        except AttributeError:
            pass

        found_tags = getattr(func, 'tags', None)
        # only include the test if the tags *are* present
        if mode == 'include':
            if found_tags is not None and found_tags & tags:
                selected.append(test)
        elif mode == 'exclude':
            # only include the test if the tags *are not* present
            if found_tags is None or not (found_tags & tags):
                selected.append(test)
        else:
            raise ValueError("Invalid 'mode' supplied: %s." % mode)
    return unittest.TestSuite(selected)


def _choose_random_tests(tests, ratio, seed):
    """
    Choose a given proportion of tests at random.
    """
    rnd = random.Random()
    rnd.seed(seed)
    if isinstance(tests, unittest.TestSuite):
        tests = _flatten_suite(tests)
    tests = rnd.sample(tests, int(len(tests) * ratio))
    tests = sorted(tests, key=lambda case: case.id())
    return unittest.TestSuite(tests)


# The reference leak detection code is liberally taken and adapted from
# Python's own Lib/test/regrtest.py.

def _refleak_cleanup():
    # Collect cyclic trash and read memory statistics immediately after.
    func1 = sys.getallocatedblocks
    try:
        func2 = sys.gettotalrefcount
    except AttributeError:
        func2 = lambda: 42

    # Flush standard output, so that buffered data is sent to the OS and
    # associated Python objects are reclaimed.
    for stream in (sys.stdout, sys.stderr, sys.__stdout__, sys.__stderr__):
        if stream is not None:
            stream.flush()

    sys._clear_type_cache()
    # This also clears the various internal CPython freelists.
    gc.collect()
    return func1(), func2()


class ReferenceLeakError(RuntimeError):
    pass


class IntPool(collections.defaultdict):

    def __missing__(self, key):
        return key


class RefleakTestResult(runner.TextTestResult):

    warmup = 3
    repetitions = 6

    def _huntLeaks(self, test):
        self.stream.flush()

        repcount = self.repetitions
        nwarmup = self.warmup
        rc_deltas = [0] * (repcount - nwarmup)
        alloc_deltas = [0] * (repcount - nwarmup)
        # Preallocate ints likely to be stored in rc_deltas and alloc_deltas,
        # to make sys.getallocatedblocks() less flaky.
        _int_pool = IntPool()
        for i in range(-200, 200):
            _int_pool[i]

        for i in range(repcount):
            # Use a pristine, silent result object to avoid recursion
            res = result.TestResult()
            test.run(res)
            # Poorly-written tests may fail when run several times.
            # In this case, abort the refleak run and report the failure.
            if not res.wasSuccessful():
                self.failures.extend(res.failures)
                self.errors.extend(res.errors)
                raise AssertionError
            del res
            alloc_after, rc_after = _refleak_cleanup()
            if i >= nwarmup:
                rc_deltas[i - nwarmup] = _int_pool[rc_after - rc_before]
                alloc_deltas[i - nwarmup] = _int_pool[alloc_after - alloc_before]
            alloc_before, rc_before = alloc_after, rc_after
        return rc_deltas, alloc_deltas

    def addSuccess(self, test):
        try:
            rc_deltas, alloc_deltas = self._huntLeaks(test)
        except AssertionError:
            # Test failed when repeated
            assert not self.wasSuccessful()
            return

        # These checkers return False on success, True on failure
        def check_rc_deltas(deltas):
            return any(deltas)

        def check_alloc_deltas(deltas):
            # At least 1/3rd of 0s
            if 3 * deltas.count(0) < len(deltas):
                return True
            # Nothing else than 1s, 0s and -1s
            if not set(deltas) <= set((1, 0, -1)):
                return True
            return False

        failed = False

        for deltas, item_name, checker in [
            (rc_deltas, 'references', check_rc_deltas),
            (alloc_deltas, 'memory blocks', check_alloc_deltas)]:
            if checker(deltas):
                msg = '%s leaked %s %s, sum=%s' % (
                    test, deltas, item_name, sum(deltas))
                failed = True
                try:
                    raise ReferenceLeakError(msg)
                except Exception:
                    exc_info = sys.exc_info()
                if self.showAll:
                    self.stream.write("%s = %r " % (item_name, deltas))
                self.addFailure(test, exc_info)

        if not failed:
            super(RefleakTestResult, self).addSuccess(test)


class RefleakTestRunner(runner.TextTestRunner):
    resultclass = RefleakTestResult


class ParallelTestResult(runner.TextTestResult):
    """
    A TestResult able to inject results from other results.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.test_runtime = 0.0
        self.records = []

    def add_results(self, result):
        """
        Add the results from the other *result* to this result.
        """
        streamout = result.stream.getvalue()
        self.records.append({
            'test_id': result.test_id,
            'stream': streamout,
            'runtime': result.test_runtime,
            'failures': result.failures,
            'errors': result.errors,
            'skipped': result.skipped,
            'expectedFailures': result.expectedFailures,
            'unexpectedSuccesses': result.unexpectedSuccesses,
            'pid': result.pid,
            'start_time': result.start_time,
        })
        self.stream.write(streamout)
        self.test_runtime += result.test_runtime
        self.stream.flush()
        self.testsRun += result.testsRun
        self.failures.extend(result.failures)
        self.errors.extend(result.errors)
        self.skipped.extend(result.skipped)
        self.expectedFailures.extend(result.expectedFailures)
        self.unexpectedSuccesses.extend(result.unexpectedSuccesses)

    def write_junit_xml(self, output):
        """
        Parameters
        ----------
        output: file or filename
        """
        from xml.etree.ElementTree import Element, SubElement, ElementTree

        suites = Element("testsuites", time=str(self.test_runtime))
        suite = SubElement(suites, "testsuite", name="numba_testsuite")
        status_to_tags = {
            'errors': 'error',
            'failures': 'failure',
            'skipped': 'skipped',
            'expectedFailures': 'skipped',
            'unexpectedSuccesses': 'failure',
        }
        # Process each record
        for rec in self.records:
            filename, classname, testname = rec['test_id'].rsplit('.', 2)

            testcase = SubElement(
                suite, "testcase",
                name='.'.join([classname, testname]),
                classname=filename,
                time=str(rec['runtime']),
            )
            props = SubElement(testcase, 'properties')
            # Write PID and start_time to help isolate the sequence of tests
            # to reach a interrupted (e.g. segfault) test
            SubElement(props, 'property', name='pid', value=str(rec['pid']))
            SubElement(props, 'property', name='start_time',
                       value=str(rec['start_time']))
            # Add error/failure/skipped tags accordingly
            for kind, tag in status_to_tags.items():
                if rec[kind]:
                    [(_id, traceback)] = rec[kind]
                    traceback = _strip_ansi_escape_sequences(traceback)
                    SubElement(testcase, tag).text = traceback
                    break
        # Write XML to output
        ElementTree(suites).write(
            output, xml_declaration=True, encoding='utf-8',
        )


def _strip_ansi_escape_sequences(text):
    """
    Adapting from https://stackoverflow.com/a/14693789
    """
    ansi_escape = re.compile(r'''
        \x1B  # ESC
        (?:   #  7-bit C1 Fe (except CSI)
            [@-Z\\-_]
        |     # or [ for CSI, followed by a control sequence
            \[
            [0-?]*  # Parameter bytes
            [ -/]*  # Intermediate bytes
            [@-~]   # Final byte
        )
    ''', re.VERBOSE)
    return ansi_escape.sub('', text)


class _MinimalResult(object):
    """
    A minimal, picklable TestResult-alike object.
    """

    __slots__ = (
        'failures', 'errors', 'skipped', 'expectedFailures',
        'unexpectedSuccesses', 'stream', 'shouldStop', 'testsRun',
        'test_id', 'test_runtime', 'start_time', 'pid')

    def fixup_case(self, case):
        """
        Remove any unpicklable attributes from TestCase instance *case*.
        """
        # Python 3.3 doesn't reset this one.
        case._outcomeForDoCleanups = None

    def __init__(self, original_result, test_id=None, test_runtime=None,
                 start_time=None, pid=None):
        for attr in self.__slots__:
            setattr(self, attr, getattr(original_result, attr, None))
        for case, _ in self.expectedFailures:
            self.fixup_case(case)
        for case, _ in self.errors:
            self.fixup_case(case)
        for case, _ in self.failures:
            self.fixup_case(case)
        self.test_id = test_id
        self.test_runtime = test_runtime
        self.start_time = start_time
        self.pid = pid


class _FakeStringIO(object):
    """
    A trivial picklable StringIO-alike for Python 2.
    """

    def __init__(self, value):
        self._value = value

    def getvalue(self):
        return self._value


class _MinimalRunner(object):
    """
    A minimal picklable object able to instantiate a runner in a
    child process and run a test case with it.
    """

    def __init__(self, runner_cls, runner_args, *, show_timing=False, write_junit=False):
        self.runner_cls = runner_cls
        self.runner_args = runner_args
        self.show_timing = show_timing
        self.write_junit = write_junit

    # Python 2 doesn't know how to pickle instance methods, so we use __call__
    # instead.

    def __call__(self, test_and_queue):
        [test, queue] = test_and_queue
        # Executed in child process
        kwargs = self.runner_args
        # Force recording of output in a buffer (it will be printed out
        # by the parent).
        kwargs['stream'] = StringIO()
        runner = self.runner_cls(**kwargs)
        result = runner._makeResult()
        # Avoid child tracebacks when Ctrl-C is pressed.
        signals.installHandler()
        signals.registerResult(result)
        result.failfast = runner.failfast
        result.buffer = runner.buffer
        start_time = time.perf_counter()
        pid = os.getpid()
        # Notify testsuite manager which tests are active
        queue.put_nowait({
            'start_time': start_time,
            'pid': pid,
            'test_id': test.id(),
        })
        with self.cleanup_object(test):
            test(result)
        end_time = time.perf_counter()
        runtime = end_time - start_time
        if self.show_timing and self.runner_args.get('verbosity', 0) > 1:
            print(f"    Runtime (seconds): {runtime}", file=result.stream)
        # HACK as cStringIO.StringIO isn't picklable in 2.x
        result.stream = _FakeStringIO(result.stream.getvalue())
        return _MinimalResult(
            result, test.id(), test_runtime=runtime, start_time=start_time,
            pid=pid,
        )

    @contextlib.contextmanager
    def cleanup_object(self, test):
        """
        A context manager which cleans up unwanted attributes on a test case
        (or any other object).
        """
        vanilla_attrs = set(test.__dict__)
        try:
            yield test
        finally:
            spurious_attrs = set(test.__dict__) - vanilla_attrs
            for name in spurious_attrs:
                del test.__dict__[name]


def _split_nonparallel_tests(test, sliced):
    """
    Split test suite into parallel and serial tests.
    """
    ptests = []
    stests = []

    flat = [*filter(sliced, _flatten_suite(test))]

    def is_parallelizable_test_case(test):
        # Guard for the fake test case created by unittest when test
        # discovery fails, as it isn't picklable (e.g. "LoadTestsFailure")
        method_name = test._testMethodName
        method = getattr(test, method_name)
        if method.__name__ != method_name and method.__name__ == "testFailure":
            return False
        # Was parallel execution explicitly disabled?
        return getattr(test, "_numba_parallel_test_", True)

    for t in flat:
        if is_parallelizable_test_case(t):
            ptests.append(t)
        else:
            stests.append(t)

    return ptests, stests

# A test can't run longer than 10 minutes
_TIMEOUT = 600

class ParallelTestRunner(runner.TextTestRunner):
    """
    A test runner which delegates the actual running to a pool of child
    processes.
    """

    resultclass = ParallelTestResult
    timeout = _TIMEOUT

    def __init__(self, runner_cls, nprocs, useslice, *, show_timing, write_junit, **kwargs):
        runner.TextTestRunner.__init__(self, **kwargs)
        self.runner_cls = runner_cls
        self.nprocs = nprocs
        self.useslice = parse_slice(useslice)
        self.runner_args = kwargs
        self.show_timing = show_timing
        self.write_junit = write_junit

    def _run_inner(self, result):
        # We hijack TextTestRunner.run()'s inner logic by passing this
        # method as if it were a test case.
        child_runner = _MinimalRunner(self.runner_cls, self.runner_args,
                                      show_timing=self.show_timing)

        # Split the tests and recycle the worker process to tame memory usage.
        chunk_size = 100
        splitted_tests = [self._ptests[i:i + chunk_size]
                          for i in range(0, len(self._ptests), chunk_size)]

        # Create the manager
        with multiprocessing.Manager() as manager:
            for tests in splitted_tests:
                pool = multiprocessing.Pool(self.nprocs)
                queue = manager.Queue()
                try:
                    self._run_parallel_tests(result, pool, child_runner, tests,
                                             queue)
                except:
                    # On exception, kill still active workers immediately
                    pool.terminate()
                    # Make sure exception is reported and not ignored
                    raise
                else:
                    # Close the pool cleanly unless asked to early out
                    if result.shouldStop:
                        pool.terminate()
                        break
                    else:
                        pool.close()
                finally:
                    # Always join the pool (this is necessary for coverage.py)
                    pool.join()
        if not result.shouldStop:
            stests = SerialSuite(self._stests)
            stests.run(result)
            return result

    def _run_parallel_tests(self, result, pool, child_runner, tests, queue):
        tests.sort(key=cuda_sensitive_mtime)

        pid_tests = collections.defaultdict(list)
        active_tests = {}

        def process_queue():
            while not queue.empty():
                data = queue.get_nowait()
                active_tests[data['test_id']] = data
                pid_tests[data['pid']].append(data['test_id'])

        it = pool.imap_unordered(child_runner, zip(tests, [queue] * len(tests)))
        while True:
            try:
                try:
                    child_result = it.__next__(self.timeout)
                finally:
                    process_queue()
            except StopIteration:
                return
            except TimeoutError as e:
                # Diagnose the names of unfinished tests
                msgbuf = [
                    "Active tests didn't finish before timeout (or crashed):"
                ]
                active_pids = []
                for tid, info in active_tests.items():
                    pid = info['pid']
                    active_pids.append(pid)
                    msgbuf.append(f"- {tid} [PID {pid}]")
                # Show test sequence in affected processes.
                # These can be used to replicate the problem in case of
                # corrupted state created by earlier tests.
                msgbuf.append("Tests ran by each affected process:")
                for pid in active_pids:
                    msgbuf.append(f"- [PID {pid}]: " + ' '.join(pid_tests[pid]))

                msg = '\n'.join(msgbuf)
                e.args = (msg,) + e.args[1:]
                raise e
            else:
                result.add_results(child_result)
                del active_tests[child_result.test_id]
                if child_result.shouldStop:
                    result.shouldStop = True
                    return

    def run(self, test):
        self._ptests, self._stests = _split_nonparallel_tests(test,
                                                              self.useslice)
        print("Parallel: %s. Serial: %s" % (len(self._ptests),
                                            len(self._stests)))
        # This will call self._run_inner() on the created result object,
        # and print out the detailed test results at the end.
        result = super(ParallelTestRunner, self).run(self._run_inner)
        # Write timing
        if self.show_timing:
            self.stream.write(f"Total test runtime (seconds): {result.test_runtime}\n")
        # Write junit xml
        if self.write_junit:
            tstamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            filename = f"junit_numba_{tstamp}.xml"
            print("Writing junit xml file:", filename)
            result.write_junit_xml(filename)
        return result
