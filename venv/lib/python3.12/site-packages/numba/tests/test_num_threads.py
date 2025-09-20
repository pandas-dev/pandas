# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, division

import sys
import os
import re
import multiprocessing
import unittest

import numpy as np

from numba import (njit, set_num_threads, get_num_threads, prange, config,
                   threading_layer, guvectorize)
from numba.np.ufunc.parallel import get_thread_id
from numba.core.errors import TypingError
from numba.tests.support import TestCase, skip_parfors_unsupported, tag
from numba.tests.test_parallel_backend import TestInSubprocess


class TestNumThreads(TestCase):
    _numba_parallel_test_ = False

    def setUp(self):
        # Make sure the num_threads is set to the max. This also makes sure
        # the threads are launched.
        set_num_threads(config.NUMBA_NUM_THREADS)

    def check_mask(self, expected, result):
        # There's no guarantee that TBB will use a full mask worth of
        # threads if it deems it inefficient to do so
        if threading_layer() == 'tbb':
            self.assertTrue(np.all(result <= expected))
        elif threading_layer() in ('omp', 'workqueue'):
            np.testing.assert_equal(expected, result)
        else:
            assert 0, 'unreachable'

    @skip_parfors_unsupported
    def test_set_num_threads_type(self):

        @njit
        def foo():
            set_num_threads('wrong_type')

        expected = "The number of threads specified must be an integer"
        for fn, errty in ((foo, TypingError), (foo.py_func, TypeError)):
            with self.assertRaises(errty) as raises:
                fn()
            self.assertIn(expected, str(raises.exception))

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_basic(self):
        max_threads = config.NUMBA_NUM_THREADS

        self.assertEqual(get_num_threads(), max_threads)
        set_num_threads(2)
        self.assertEqual(get_num_threads(), 2)
        set_num_threads(max_threads)
        self.assertEqual(get_num_threads(), max_threads)

        with self.assertRaises(ValueError):
            set_num_threads(0)

        with self.assertRaises(ValueError):
            set_num_threads(max_threads + 1)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_basic_jit(self):
        max_threads = config.NUMBA_NUM_THREADS

        @njit
        def get_n():
            return get_num_threads()

        self.assertEqual(get_n(), max_threads)
        set_num_threads(2)
        self.assertEqual(get_n(), 2)
        set_num_threads(max_threads)
        self.assertEqual(get_n(), max_threads)

        @njit
        def set_get_n(n):
            set_num_threads(n)
            return get_num_threads()

        self.assertEqual(set_get_n(2), 2)
        self.assertEqual(set_get_n(max_threads), max_threads)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_basic_guvectorize(self):
        max_threads = config.NUMBA_NUM_THREADS

        @guvectorize(['void(int64[:])'],
                     '(n)',
                     nopython=True,
                     target='parallel')
        def get_n(x):
            x[:] = get_num_threads()

        x = np.zeros((5000000,), dtype=np.int64)
        get_n(x)
        np.testing.assert_equal(x, max_threads)
        set_num_threads(2)
        x = np.zeros((5000000,), dtype=np.int64)
        get_n(x)
        np.testing.assert_equal(x, 2)
        set_num_threads(max_threads)
        x = np.zeros((5000000,), dtype=np.int64)
        get_n(x)
        np.testing.assert_equal(x, max_threads)

        @guvectorize(['void(int64[:])'],
                     '(n)',
                     nopython=True,
                     target='parallel')
        def set_get_n(n):
            set_num_threads(n[0])
            n[:] = get_num_threads()

        x = np.zeros((5000000,), dtype=np.int64)
        x[0] = 2
        set_get_n(x)
        np.testing.assert_equal(x, 2)
        x = np.zeros((5000000,), dtype=np.int64)
        x[0] = max_threads
        set_get_n(x)
        np.testing.assert_equal(x, max_threads)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_outside_jit(self):

        # Test set_num_threads outside a jitted function
        set_num_threads(2)

        @njit(parallel=True)
        def test_func():
            x = 5
            buf = np.empty((x,))
            for i in prange(x):
                buf[i] = get_num_threads()
            return buf

        @guvectorize(['void(int64[:])'],
                     '(n)',
                     nopython=True,
                     target='parallel')
        def test_gufunc(x):
            x[:] = get_num_threads()

        out = test_func()
        np.testing.assert_equal(out, 2)

        x = np.zeros((5000000,), dtype=np.int64)
        test_gufunc(x)
        np.testing.assert_equal(x, 2)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_inside_jit(self):
        # Test set_num_threads inside a jitted function
        @njit(parallel=True)
        def test_func(nthreads):
            x = 5
            buf = np.empty((x,))
            set_num_threads(nthreads)
            for i in prange(x):
                buf[i] = get_num_threads()
            return buf

        mask = 2
        out = test_func(mask)
        np.testing.assert_equal(out, mask)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_set_num_threads_inside_guvectorize(self):
        # Test set_num_threads inside a jitted guvectorize function
        @guvectorize(['void(int64[:])'],
                     '(n)',
                     nopython=True,
                     target='parallel')
        def test_func(x):
            set_num_threads(x[0])
            x[:] = get_num_threads()

        x = np.zeros((5000000,), dtype=np.int64)
        mask = 2
        x[0] = mask
        test_func(x)
        np.testing.assert_equal(x, mask)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_get_num_threads_truth_outside_jit(self):

        for mask in range(2, min(6, config.NUMBA_NUM_THREADS + 1)):
            set_num_threads(mask)

            # a lot of work, hopefully will trigger "mask" count of threads to
            # join the parallel region (for those backends with dynamic threads)
            @njit(parallel=True)
            def test_func():
                x = 5000000
                buf = np.empty((x,))
                for i in prange(x):
                    buf[i] = get_thread_id()
                return len(np.unique(buf)), get_num_threads()

            out = test_func()
            self.check_mask((mask, mask), out)

            @guvectorize(['void(int64[:], int64[:])'],
                         '(n), (m)',
                         nopython=True,
                         target='parallel')
            def test_gufunc(x, out):
                x[:] = get_thread_id()
                out[0] = get_num_threads()

            # Reshape to force parallelism
            x = np.full((5000000,), -1, dtype=np.int64).reshape((100, 50000))
            out = np.zeros((1,), dtype=np.int64)
            test_gufunc(x, out)
            self.check_mask(mask, out)
            self.check_mask(mask, len(np.unique(x)))

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_get_num_threads_truth_inside_jit(self):

        for mask in range(2, min(6, config.NUMBA_NUM_THREADS + 1)):

            # a lot of work, hopefully will trigger "mask" count of threads to
            # join the parallel region (for those backends with dynamic threads)
            @njit(parallel=True)
            def test_func():
                set_num_threads(mask)
                x = 5000000
                buf = np.empty((x,))
                for i in prange(x):
                    buf[i] = get_thread_id()
                return len(np.unique(buf)), get_num_threads()

            out = test_func()
            self.check_mask((mask, mask), out)

            @guvectorize(['void(int64[:], int64[:])'],
                         '(n), (m)',
                         nopython=True,
                         target='parallel')
            def test_gufunc(x, out):
                set_num_threads(mask)
                x[:] = get_thread_id()
                out[0] = get_num_threads()

            # Reshape to force parallelism
            x = np.full((5000000,), -1, dtype=np.int64).reshape((100, 50000))
            out = np.zeros((1,), dtype=np.int64)
            test_gufunc(x, out)
            self.check_mask(mask, out)
            self.check_mask(mask, len(np.unique(x)))

    # this test can only run on OpenMP (providing OMP_MAX_ACTIVE_LEVELS is not
    # set or >= 2) and TBB backends
    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_nested_parallelism_1(self):
        if threading_layer() == 'workqueue':
            self.skipTest("workqueue is not threadsafe")

        # check that get_num_threads is ok in nesting
        mask = config.NUMBA_NUM_THREADS - 1

        N = config.NUMBA_NUM_THREADS
        M = 2 * config.NUMBA_NUM_THREADS

        @njit(parallel=True)
        def child_func(buf, fid):
            M, N = buf.shape
            for i in prange(N):
                buf[fid, i] = get_num_threads()

        def get_test(test_type):
            if test_type == 'njit':
                def test_func(nthreads, py_func=False):
                    @njit(parallel=True)
                    def _test_func(nthreads):
                        acc = 0
                        buf = np.zeros((M, N))
                        set_num_threads(nthreads)
                        for i in prange(M):
                            local_mask = 1 + i % mask
                            # set threads in parent function
                            set_num_threads(local_mask)
                            if local_mask < N:
                                child_func(buf, local_mask)
                            acc += get_num_threads()
                        return acc, buf
                    if py_func:
                        return _test_func.py_func(nthreads)
                    else:
                        return _test_func(nthreads)

            elif test_type == 'guvectorize':
                def test_func(nthreads, py_func=False):
                    def _test_func(acc, buf, local_mask):
                        set_num_threads(nthreads)
                        # set threads in parent function
                        set_num_threads(local_mask[0])
                        if local_mask[0] < N:
                            child_func(buf, local_mask[0])
                        acc[0] += get_num_threads()

                    buf = np.zeros((M, N), dtype=np.int64)
                    acc = np.zeros((M, 1), dtype=np.int64)
                    local_mask = (1 + np.arange(M) % mask).reshape((M, 1))
                    sig = ['void(int64[:], int64[:, :], int64[:])']
                    layout = '(p), (n, m), (p)'
                    if not py_func:
                        _test_func = guvectorize(sig, layout, nopython=True,
                                                 target='parallel')(_test_func)
                    else:
                        _test_func = guvectorize(sig, layout,
                                                 forceobj=True)(_test_func)
                    _test_func(acc, buf, local_mask)
                    return acc, buf

            return test_func

        for test_type in ['njit', 'guvectorize']:
            test_func = get_test(test_type)
            got_acc, got_arr = test_func(mask)
            exp_acc, exp_arr = test_func(mask, py_func=True)
            np.testing.assert_equal(exp_acc, got_acc)
            np.testing.assert_equal(exp_arr, got_arr)

            # check the maths reconciles, guvectorize does not reduce, njit does
            math_acc_exp = 1 + np.arange(M) % mask
            if test_type == 'guvectorize':
                math_acc = math_acc_exp.reshape((M, 1))
            else:
                math_acc = np.sum(math_acc_exp)

            np.testing.assert_equal(math_acc, got_acc)

            math_arr = np.zeros((M, N))
            for i in range(1, N):
                # there's branches on 1, ..., num_threads - 1
                math_arr[i, :] = i
            np.testing.assert_equal(math_arr, got_arr)

    # this test can only run on OpenMP (providing OMP_MAX_ACTIVE_LEVELS is not
    # set or >= 2) and TBB backends
    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    def _test_nested_parallelism_2(self):
        if threading_layer() == 'workqueue':
            self.skipTest("workqueue is not threadsafe")

        # check that get_num_threads is ok in nesting

        N = config.NUMBA_NUM_THREADS + 1
        M = 4 * config.NUMBA_NUM_THREADS + 1

        def get_impl(child_type, test_type):

            if child_type == 'parallel':
                child_dec = njit(parallel=True)
            elif child_type == 'njit':
                child_dec = njit(parallel=False)
            elif child_type == 'none':
                def child_dec(x):
                    return x

            @child_dec
            def child(buf, fid):
                M, N = buf.shape
                set_num_threads(fid)  # set threads in child function
                for i in prange(N):
                    buf[fid, i] = get_num_threads()

            if test_type in ['parallel', 'njit', 'none']:
                if test_type == 'parallel':
                    test_dec = njit(parallel=True)
                elif test_type == 'njit':
                    test_dec = njit(parallel=False)
                elif test_type == 'none':
                    def test_dec(x):
                        return x

                @test_dec
                def test_func(nthreads):
                    buf = np.zeros((M, N))
                    set_num_threads(nthreads)
                    for i in prange(M):
                        local_mask = 1 + i % mask
                        # when the threads exit the child functions they should
                        # have a TLS slot value of the local mask as it was set
                        # in child
                        if local_mask < config.NUMBA_NUM_THREADS:
                            child(buf, local_mask)
                            assert get_num_threads() == local_mask
                    return buf
            else:
                if test_type == 'guvectorize':
                    test_dec = guvectorize(['int64[:,:], int64[:]'],
                                           '(n, m), (k)', nopython=True,
                                           target='parallel')
                elif test_type == 'guvectorize-obj':
                    test_dec = guvectorize(['int64[:,:], int64[:]'],
                                           '(n, m), (k)', forceobj=True)

                def test_func(nthreads):
                    @test_dec
                    def _test_func(buf, local_mask):
                        set_num_threads(nthreads)
                        # when the threads exit the child functions they should
                        # have a TLS slot value of the local mask as it was set
                        # in child
                        if local_mask[0] < config.NUMBA_NUM_THREADS:
                            child(buf, local_mask[0])
                            assert get_num_threads() == local_mask[0]

                    buf = np.zeros((M, N), dtype=np.int64)
                    local_mask = (1 + np.arange(M) % mask).reshape((M, 1))
                    _test_func(buf, local_mask)
                    return buf

            return test_func

        mask = config.NUMBA_NUM_THREADS - 1

        res_arrays = {}
        for test_type in ['parallel', 'njit', 'none',
                          'guvectorize', 'guvectorize-obj']:
            for child_type in ['parallel', 'njit', 'none']:
                if child_type == 'none' and test_type != 'none':
                    continue
                set_num_threads(mask)
                res_arrays[test_type, child_type] = get_impl(
                    child_type, test_type)(mask)

        py_arr = res_arrays['none', 'none']
        for arr in res_arrays.values():
            np.testing.assert_equal(arr, py_arr)

        # check the maths reconciles
        math_arr = np.zeros((M, N))
        # there's branches on modulo mask but only NUMBA_NUM_THREADS funcs
        for i in range(1, config.NUMBA_NUM_THREADS):
            math_arr[i, :] = i

        np.testing.assert_equal(math_arr, py_arr)

    # this test can only run on OpenMP (providing OMP_MAX_ACTIVE_LEVELS is not
    # set or >= 2) and TBB backends
    # This test needs at least 3 threads to run, N>=2 for the launch, M>=N+1 for
    # the nested function
    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 3, "Not enough CPU cores")
    def _test_nested_parallelism_3(self):
        if threading_layer() == 'workqueue':
            self.skipTest("workqueue is not threadsafe")

        # check that the right number of threads are present in nesting
        # this relies on there being a load of cores present
        BIG = 1000000

        @njit(parallel=True)
        def work(local_nt):  # arg is value 3
            tid = np.zeros(BIG)
            acc = 0
            set_num_threads(local_nt)  # set to 3 threads
            for i in prange(BIG):
                acc += 1
                tid[i] = get_thread_id()
            return acc, np.unique(tid)

        @njit(parallel=True)
        def test_func_jit(nthreads):
            set_num_threads(nthreads) # set to 2 threads
            lens = np.zeros(nthreads)
            total = 0
            for i in prange(nthreads):
                my_acc, tids = work(nthreads + 1)  # call with value 3
                lens[i] = len(tids)
                total += my_acc
            return total, np.unique(lens)

        NT = 2
        expected_acc = BIG * NT
        expected_thread_count = NT + 1

        got_acc, got_tc = test_func_jit(NT)
        self.assertEqual(expected_acc, got_acc)
        self.check_mask(expected_thread_count, got_tc)

        def test_guvectorize(nthreads):
            @guvectorize(['int64[:], int64[:]'],
                         '(n), (n)',
                         nopython=True,
                         target='parallel')
            def test_func_guvectorize(total, lens):
                my_acc, tids = work(nthreads + 1)
                lens[0] = len(tids)
                total[0] += my_acc

            total = np.zeros((nthreads, 1), dtype=np.int64)
            lens = np.zeros(nthreads, dtype=np.int64).reshape((nthreads, 1))

            test_func_guvectorize(total, lens)
            # vectorize does not reduce, so total is summed
            return total.sum(), np.unique(lens)

        got_acc, got_tc = test_guvectorize(NT)

        self.assertEqual(expected_acc, got_acc)
        self.check_mask(expected_thread_count, got_tc)

    @skip_parfors_unsupported
    @unittest.skipIf(config.NUMBA_NUM_THREADS < 2, "Not enough CPU cores")
    @unittest.skipIf(not sys.platform.startswith('linux'), "Linux only")
    def _test_threadmask_across_fork(self):
        forkctx = multiprocessing.get_context('fork')

        @njit
        def foo():
            return get_num_threads()

        def wrap(queue):
            queue.put(foo())

        mask = 1
        self.assertEqual(foo(), config.NUMBA_NUM_THREADS)
        set_num_threads(mask)
        self.assertEqual(foo(), mask)
        shared_queue = forkctx.Queue()
        # check TLS slot inheritance in fork
        p = forkctx.Process(target=wrap, args=(shared_queue,))
        p.start()
        p.join()
        self.assertEqual(shared_queue.get(), mask)

    def tearDown(self):
        set_num_threads(config.NUMBA_NUM_THREADS)

    @skip_parfors_unsupported
    def _test_get_thread_id_not_parallel(self):
        python_get_thread_id = get_thread_id()
        check_array_size = 8

        @njit(parallel=False)
        def par_false(size):
            njit_par_false_tid = get_thread_id()
            res = np.ones(size)
            for i in prange(size):
                res[i] = get_thread_id()
            return njit_par_false_tid, res

        @njit(parallel=True)
        def par_true(size):
            njit_par_true_tid = get_thread_id()
            res = np.ones(size)
            for i in range(size):
                res[i] = get_thread_id()
            return njit_par_true_tid, res

        self.assertEqual(python_get_thread_id, 0)
        njit_par_false_tid, njit_par_false_arr = par_false(check_array_size)
        self.assertEqual(njit_par_false_tid, 0)
        np.testing.assert_equal(njit_par_false_arr, 0)
        njit_par_true_tid, njit_par_true_arr = par_true(check_array_size)
        self.assertEqual(njit_par_true_tid, 0)
        np.testing.assert_equal(njit_par_true_arr, 0)


class TestNumThreadsBackends(TestInSubprocess, TestCase):
    _class = TestNumThreads
    _DEBUG = False

    # 1 is mainly here to ensure tests skip correctly
    num_threads = [i for i in [1, 2, 4, 8, 16] if i <= config.NUMBA_NUM_THREADS]

    def run_test_in_separate_process(self, test, threading_layer, num_threads):
        env_copy = os.environ.copy()
        env_copy['NUMBA_THREADING_LAYER'] = str(threading_layer)
        env_copy['NUMBA_NUM_THREADS'] = str(num_threads)
        cmdline = [sys.executable, "-m", "numba.runtests", "-v", test]
        return self.run_cmd(cmdline, env_copy)

    @classmethod
    def _inject(cls, name, backend, backend_guard, num_threads):
        themod = cls.__module__
        thecls = cls._class.__name__
        injected_method = '%s.%s.%s' % (themod, thecls, name)

        def test_template(self):
            o, e = self.run_test_in_separate_process(injected_method, backend,
                                                     num_threads)
            if self._DEBUG:
                print('stdout:\n "%s"\n stderr:\n "%s"' % (o, e))
            # If the test was skipped in the subprocess, then mark this as a
            # skipped test.
            m = re.search(r"\.\.\. skipped '(.*?)'", e)
            if m is not None:
                self.skipTest(m.group(1))
            self.assertIn('OK', e)
            self.assertTrue('FAIL' not in e)
            self.assertTrue('ERROR' not in e)

        injected_test = "%s_%s_%s_threads" % (name[1:], backend, num_threads)
        setattr(cls, injected_test,
                tag('long_running')(backend_guard(test_template)))

    @classmethod
    def generate(cls):
        for name in cls._class.__dict__.copy():
            for backend, backend_guard in cls.backends.items():
                for num_threads in cls.num_threads:
                    if not name.startswith('_test_'):
                        continue
                    cls._inject(name, backend, backend_guard, num_threads)


TestNumThreadsBackends.generate()

if __name__ == '__main__':
    unittest.main()
