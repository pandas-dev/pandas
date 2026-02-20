import copy
import itertools
import math
import random
import sys
import unittest

import numpy as np

from numba import jit, njit
from numba.core import utils, errors
from numba.tests.support import TestCase, MemoryLeakMixin

from numba.misc.quicksort import make_py_quicksort, make_jit_quicksort
from numba.misc.mergesort import make_jit_mergesort
from numba.misc.timsort import make_py_timsort, make_jit_timsort, MergeRun


def make_temp_list(keys, n):
    return [keys[0]] * n

def make_temp_array(keys, n):
    return np.empty(n, keys.dtype)


py_list_timsort = make_py_timsort(make_temp_list)

py_array_timsort = make_py_timsort(make_temp_array)

jit_list_timsort = make_jit_timsort(make_temp_list)

jit_array_timsort = make_jit_timsort(make_temp_array)

py_quicksort = make_py_quicksort()

jit_quicksort = make_jit_quicksort()


def sort_usecase(val):
    val.sort()

def argsort_usecase(val):
    return val.argsort()

def argsort_kind_usecase(val, is_stable=False):
    if is_stable:
        return val.argsort(kind='mergesort')
    else:
        return val.argsort(kind='quicksort')

def sorted_usecase(val):
    return sorted(val)

def sorted_reverse_usecase(val, b):
    return sorted(val, reverse=b)

def np_sort_usecase(val):
    return np.sort(val)

def np_argsort_usecase(val):
    return np.argsort(val)

def np_argsort_kind_usecase(val, is_stable=False):
    if is_stable:
        return np.argsort(val, kind='mergesort')
    else:
        return np.argsort(val, kind='quicksort')

def list_sort_usecase(n):
    np.random.seed(42)
    l = []
    for i in range(n):
        l.append(np.random.random())
    ll = l[:]
    ll.sort()
    return l, ll

def list_sort_reverse_usecase(n, b):
    np.random.seed(42)
    l = []
    for i in range(n):
        l.append(np.random.random())
    ll = l[:]
    ll.sort(reverse=b)
    return l, ll


class BaseSortingTest(object):

    def random_list(self, n, offset=10):
        random.seed(42)
        l = list(range(offset, offset + n))
        random.shuffle(l)
        return l

    def sorted_list(self, n, offset=10):
        return list(range(offset, offset + n))

    def revsorted_list(self, n, offset=10):
        return list(range(offset, offset + n))[::-1]

    def initially_sorted_list(self, n, m=None, offset=10):
        if m is None:
            m = n // 2
        l = self.sorted_list(m, offset)
        l += self.random_list(n - m, offset=l[-1] + offset)
        return l

    def duprandom_list(self, n, factor=None, offset=10):
        random.seed(42)
        if factor is None:
            factor = int(math.sqrt(n))
        l = (list(range(offset, offset + (n // factor) + 1)) * (factor + 1))[:n]
        assert len(l) == n
        random.shuffle(l)
        return l

    def dupsorted_list(self, n, factor=None, offset=10):
        if factor is None:
            factor = int(math.sqrt(n))
        l = (list(range(offset, offset + (n // factor) + 1)) * (factor + 1))[:n]
        assert len(l) == n, (len(l), n)
        l.sort()
        return l

    def assertSorted(self, orig, result):
        self.assertEqual(len(result), len(orig))
        # sorted() returns a list, so make sure we compare to another list
        self.assertEqual(list(result), sorted(orig))

    def assertSortedValues(self, orig, orig_values, result, result_values):
        self.assertEqual(len(result), len(orig))
        self.assertEqual(list(result), sorted(orig))
        zip_sorted = sorted(zip(orig, orig_values), key=lambda x: x[0])
        zip_result = list(zip(result, result_values))
        self.assertEqual(zip_sorted, zip_result)
        # Check stability
        for i in range(len(zip_result) - 1):
            (k1, v1), (k2, v2) = zip_result[i], zip_result[i + 1]
            if k1 == k2:
                # Assuming values are unique, which is enforced by the tests
                self.assertLess(orig_values.index(v1), orig_values.index(v2))

    def fibo(self):
        a = 1
        b = 1
        while True:
            yield a
            a, b = b, a + b

    def make_sample_sorted_lists(self, n):
        lists = []
        for offset in (20, 120):
            lists.append(self.sorted_list(n, offset))
            lists.append(self.dupsorted_list(n, offset))
        return lists

    def make_sample_lists(self, n):
        lists = []
        for offset in (20, 120):
            lists.append(self.sorted_list(n, offset))
            lists.append(self.dupsorted_list(n, offset))
            lists.append(self.revsorted_list(n, offset))
            lists.append(self.duprandom_list(n, offset))
        return lists


class BaseTimsortTest(BaseSortingTest):

    def merge_init(self, keys):
        f = self.timsort.merge_init
        return f(keys)

    def test_binarysort(self):
        n = 20
        def check(l, n, start=0):
            res = self.array_factory(l)
            f(res, res, 0, n, start)
            self.assertSorted(l, res)

        f = self.timsort.binarysort
        l = self.sorted_list(n)
        check(l, n)
        check(l, n, n//2)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.initially_sorted_list(n, n//2)
        check(l, n)
        check(l, n, n//2)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.random_list(n)
        check(l, n)
        l = self.duprandom_list(n)
        check(l, n)

    def test_binarysort_with_values(self):
        n = 20
        v = list(range(100, 100+n))

        def check(l, n, start=0):
            res = self.array_factory(l)
            res_v = self.array_factory(v)
            f(res, res_v, 0, n, start)
            self.assertSortedValues(l, v, res, res_v)

        f = self.timsort.binarysort
        l = self.sorted_list(n)
        check(l, n)
        check(l, n, n//2)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.initially_sorted_list(n, n//2)
        check(l, n)
        check(l, n, n//2)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.random_list(n)
        check(l, n)
        l = self.duprandom_list(n)
        check(l, n)

    def test_count_run(self):
        n = 16
        f = self.timsort.count_run

        def check(l, lo, hi):
            n, desc = f(self.array_factory(l), lo, hi)
            # Fully check invariants
            if desc:
                for k in range(lo, lo + n - 1):
                    a, b = l[k], l[k + 1]
                    self.assertGreater(a, b)
                if lo + n < hi:
                    self.assertLessEqual(l[lo + n - 1], l[lo + n])
            else:
                for k in range(lo, lo + n - 1):
                    a, b = l[k], l[k + 1]
                    self.assertLessEqual(a, b)
                if lo + n < hi:
                    self.assertGreater(l[lo + n - 1], l[lo + n], l)


        l = self.sorted_list(n, offset=100)
        check(l, 0, n)
        check(l, 1, n - 1)
        check(l, 1, 2)
        l = self.revsorted_list(n, offset=100)
        check(l, 0, n)
        check(l, 1, n - 1)
        check(l, 1, 2)
        l = self.random_list(n, offset=100)
        for i in range(len(l) - 1):
            check(l, i, n)
        l = self.duprandom_list(n, offset=100)
        for i in range(len(l) - 1):
            check(l, i, n)

    def test_gallop_left(self):
        n = 20
        f = self.timsort.gallop_left

        def check(l, key, start, stop, hint):
            k = f(key, l, start, stop, hint)
            # Fully check invariants
            self.assertGreaterEqual(k, start)
            self.assertLessEqual(k, stop)
            if k > start:
                self.assertLess(l[k - 1], key)
            if k < stop:
                self.assertGreaterEqual(l[k], key)

        def check_all_hints(l, key, start, stop):
            for hint in range(start, stop):
                check(l, key, start, stop, hint)

        def check_sorted_list(l):
            l = self.array_factory(l)
            for key in (l[5], l[15], l[0], -1000, l[-1], 1000):
                check_all_hints(l, key, 0, n)
                check_all_hints(l, key, 1, n - 1)
                check_all_hints(l, key, 8, n - 8)

        l = self.sorted_list(n, offset=100)
        check_sorted_list(l)
        l = self.dupsorted_list(n, offset=100)
        check_sorted_list(l)

    def test_gallop_right(self):
        n = 20
        f = self.timsort.gallop_right

        def check(l, key, start, stop, hint):
            k = f(key, l, start, stop, hint)
            # Fully check invariants
            self.assertGreaterEqual(k, start)
            self.assertLessEqual(k, stop)
            if k > start:
                self.assertLessEqual(l[k - 1], key)
            if k < stop:
                self.assertGreater(l[k], key)

        def check_all_hints(l, key, start, stop):
            for hint in range(start, stop):
                check(l, key, start, stop, hint)

        def check_sorted_list(l):
            l = self.array_factory(l)
            for key in (l[5], l[15], l[0], -1000, l[-1], 1000):
                check_all_hints(l, key, 0, n)
                check_all_hints(l, key, 1, n - 1)
                check_all_hints(l, key, 8, n - 8)

        l = self.sorted_list(n, offset=100)
        check_sorted_list(l)
        l = self.dupsorted_list(n, offset=100)
        check_sorted_list(l)

    def test_merge_compute_minrun(self):
        f = self.timsort.merge_compute_minrun

        for i in range(0, 64):
            self.assertEqual(f(i), i)
        for i in range(6, 63):
            if 2**i > sys.maxsize:
                break
            self.assertEqual(f(2**i), 32)
        for i in self.fibo():
            if i < 64:
                continue
            if i >= sys.maxsize:
                break
            k = f(i)
            self.assertGreaterEqual(k, 32)
            self.assertLessEqual(k, 64)
            if i > 500:
                # i/k is close to, but strictly less than, an exact power of 2
                quot = i // k
                p = 2 ** utils.bit_length(quot)
                self.assertLess(quot, p)
                self.assertGreaterEqual(quot, 0.9 * p)

    def check_merge_lo_hi(self, func, a, b):
        na = len(a)
        nb = len(b)

        # Add sentinels at start and end, to check they weren't moved
        orig_keys = [42] + a + b + [-42]
        keys = self.array_factory(orig_keys)
        ms = self.merge_init(keys)
        ssa = 1
        ssb = ssa + na

        #new_ms = func(ms, keys, [], ssa, na, ssb, nb)
        new_ms = func(ms, keys, keys, ssa, na, ssb, nb)
        self.assertEqual(keys[0], orig_keys[0])
        self.assertEqual(keys[-1], orig_keys[-1])
        self.assertSorted(orig_keys[1:-1], keys[1:-1])
        # Check the MergeState result
        self.assertGreaterEqual(len(new_ms.keys), len(ms.keys))
        self.assertGreaterEqual(len(new_ms.values), len(ms.values))
        self.assertIs(new_ms.pending, ms.pending)
        self.assertGreaterEqual(new_ms.min_gallop, 1)

    def test_merge_lo_hi(self):
        f_lo = self.timsort.merge_lo
        f_hi = self.timsort.merge_hi

        # The larger sizes exercise galloping
        for (na, nb) in [(12, 16), (40, 40), (100, 110), (1000, 1100)]:
            for a, b in itertools.product(self.make_sample_sorted_lists(na),
                                          self.make_sample_sorted_lists(nb)):
                self.check_merge_lo_hi(f_lo, a, b)
                self.check_merge_lo_hi(f_hi, b, a)

    def check_merge_at(self, a, b):
        f = self.timsort.merge_at
        # Prepare the array to be sorted
        na = len(a)
        nb = len(b)
        # Add sentinels at start and end, to check they weren't moved
        orig_keys = [42] + a + b + [-42]
        ssa = 1
        ssb = ssa + na

        stack_sentinel = MergeRun(-42, -42)

        def run_merge_at(ms, keys, i):
            new_ms = f(ms, keys, keys, i)
            self.assertEqual(keys[0], orig_keys[0])
            self.assertEqual(keys[-1], orig_keys[-1])
            self.assertSorted(orig_keys[1:-1], keys[1:-1])
            # Check stack state
            self.assertIs(new_ms.pending, ms.pending)
            self.assertEqual(ms.pending[i], (ssa, na + nb))
            self.assertEqual(ms.pending[0], stack_sentinel)
            return new_ms

        # First check with i == len(stack) - 2
        keys = self.array_factory(orig_keys)
        ms = self.merge_init(keys)
        # Push sentinel on stack, to check it wasn't touched
        ms = self.timsort.merge_append(ms, stack_sentinel)
        i = ms.n
        ms = self.timsort.merge_append(ms, MergeRun(ssa, na))
        ms = self.timsort.merge_append(ms, MergeRun(ssb, nb))
        ms = run_merge_at(ms, keys, i)
        self.assertEqual(ms.n, i + 1)

        # Now check with i == len(stack) - 3
        keys = self.array_factory(orig_keys)
        ms = self.merge_init(keys)
        # Push sentinel on stack, to check it wasn't touched
        ms = self.timsort.merge_append(ms, stack_sentinel)
        i = ms.n
        ms = self.timsort.merge_append(ms, MergeRun(ssa, na))
        ms = self.timsort.merge_append(ms, MergeRun(ssb, nb))
        # A last run (trivial here)
        last_run = MergeRun(ssb + nb, 1)
        ms = self.timsort.merge_append(ms, last_run)
        ms = run_merge_at(ms, keys, i)
        self.assertEqual(ms.n, i + 2)
        self.assertEqual(ms.pending[ms.n - 1], last_run)

    def test_merge_at(self):
        # The larger sizes exercise galloping
        for (na, nb) in [(12, 16), (40, 40), (100, 110), (500, 510)]:
            for a, b in itertools.product(self.make_sample_sorted_lists(na),
                                          self.make_sample_sorted_lists(nb)):
                self.check_merge_at(a, b)
                self.check_merge_at(b, a)

    def test_merge_force_collapse(self):
        f = self.timsort.merge_force_collapse

        # Test with runs of ascending sizes, then descending sizes
        sizes_list = [(8, 10, 15, 20)]
        sizes_list.append(sizes_list[0][::-1])

        for sizes in sizes_list:
            for chunks in itertools.product(*(self.make_sample_sorted_lists(n)
                                              for n in sizes)):
                # Create runs of the given sizes
                orig_keys = sum(chunks, [])
                keys = self.array_factory(orig_keys)
                ms = self.merge_init(keys)
                pos = 0
                for c in chunks:
                    ms = self.timsort.merge_append(ms, MergeRun(pos, len(c)))
                    pos += len(c)
                # Sanity check
                self.assertEqual(sum(ms.pending[ms.n - 1]), len(keys))
                # Now merge the runs
                ms = f(ms, keys, keys)
                # Remaining run is the whole list
                self.assertEqual(ms.n, 1)
                self.assertEqual(ms.pending[0], MergeRun(0, len(keys)))
                # The list is now sorted
                self.assertSorted(orig_keys, keys)

    def test_run_timsort(self):
        f = self.timsort.run_timsort

        for size_factor in (1, 10):
            # Make lists to be sorted from three chunks of different kinds.
            sizes = (15, 30, 20)

            all_lists = [self.make_sample_lists(n * size_factor) for n in sizes]
            for chunks in itertools.product(*all_lists):
                orig_keys = sum(chunks, [])
                keys = self.array_factory(orig_keys)
                f(keys)
                # The list is now sorted
                self.assertSorted(orig_keys, keys)

    def test_run_timsort_with_values(self):
        # Run timsort, but also with a values array
        f = self.timsort.run_timsort_with_values

        for size_factor in (1, 5):
            chunk_size = 80 * size_factor
            a = self.dupsorted_list(chunk_size)
            b = self.duprandom_list(chunk_size)
            c = self.revsorted_list(chunk_size)
            orig_keys = a + b + c
            orig_values = list(range(1000, 1000 + len(orig_keys)))

            keys = self.array_factory(orig_keys)
            values = self.array_factory(orig_values)
            f(keys, values)
            # This checks sort stability
            self.assertSortedValues(orig_keys, orig_values, keys, values)


class TestTimsortPurePython(BaseTimsortTest, TestCase):

    timsort = py_list_timsort

    # Much faster than a Numpy array in pure Python
    array_factory = list


class TestTimsortArraysPurePython(BaseTimsortTest, TestCase):

    timsort = py_array_timsort

    def array_factory(self, lst):
        return np.array(lst, dtype=np.int32)


class JITTimsortMixin(object):

    timsort = jit_array_timsort

    test_merge_at = None
    test_merge_force_collapse = None

    def wrap_with_mergestate(self, timsort, func, _cache=None):
        """
        Wrap *func* into another compiled function inserting a runtime-created
        mergestate as the first function argument.
        """
        if _cache is None:
            _cache = {}
        key = timsort, func
        if key in _cache:
            return _cache[key]

        merge_init = timsort.merge_init

        @timsort.compile
        def wrapper(keys, values, *args):
            ms = merge_init(keys)
            res = func(ms, keys, values, *args)
            return res

        _cache[key] = wrapper
        return wrapper


class TestTimsortArrays(JITTimsortMixin, BaseTimsortTest, TestCase):

    def array_factory(self, lst):
        return np.array(lst, dtype=np.int32)

    def check_merge_lo_hi(self, func, a, b):
        na = len(a)
        nb = len(b)

        func = self.wrap_with_mergestate(self.timsort, func)

        # Add sentinels at start and end, to check they weren't moved
        orig_keys = [42] + a + b + [-42]
        keys = self.array_factory(orig_keys)
        ssa = 1
        ssb = ssa + na

        new_ms = func(keys, keys, ssa, na, ssb, nb)
        self.assertEqual(keys[0], orig_keys[0])
        self.assertEqual(keys[-1], orig_keys[-1])
        self.assertSorted(orig_keys[1:-1], keys[1:-1])



class BaseQuicksortTest(BaseSortingTest):

    # NOTE these tests assume a non-argsort quicksort.

    def test_insertion_sort(self):
        n = 20
        def check(l, n):
            res = self.array_factory([9999] + l + [-9999])
            f(res, res, 1, n)
            self.assertEqual(res[0], 9999)
            self.assertEqual(res[-1], -9999)
            self.assertSorted(l, res[1:-1])

        f = self.quicksort.insertion_sort
        l = self.sorted_list(n)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.initially_sorted_list(n, n//2)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.random_list(n)
        check(l, n)
        l = self.duprandom_list(n)
        check(l, n)

    def test_partition(self):
        n = 20
        def check(l, n):
            res = self.array_factory([9999] + l + [-9999])
            index = f(res, res, 1, n)
            self.assertEqual(res[0], 9999)
            self.assertEqual(res[-1], -9999)
            pivot = res[index]
            for i in range(1, index):
                self.assertLessEqual(res[i], pivot)
            for i in range(index + 1, n):
                self.assertGreaterEqual(res[i], pivot)

        f = self.quicksort.partition
        l = self.sorted_list(n)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.initially_sorted_list(n, n//2)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.random_list(n)
        check(l, n)
        l = self.duprandom_list(n)
        check(l, n)

    def test_partition3(self):
        # Test the unused partition3() function
        n = 20
        def check(l, n):
            res = self.array_factory([9999] + l + [-9999])
            lt, gt = f(res, 1, n)
            self.assertEqual(res[0], 9999)
            self.assertEqual(res[-1], -9999)
            pivot = res[lt]
            for i in range(1, lt):
                self.assertLessEqual(res[i], pivot)
            for i in range(lt, gt + 1):
                self.assertEqual(res[i], pivot)
            for i in range(gt + 1, n):
                self.assertGreater(res[i], pivot)

        f = self.quicksort.partition3
        l = self.sorted_list(n)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.initially_sorted_list(n, n//2)
        check(l, n)
        l = self.revsorted_list(n)
        check(l, n)
        l = self.random_list(n)
        check(l, n)
        l = self.duprandom_list(n)
        check(l, n)

    def test_run_quicksort(self):
        f = self.quicksort.run_quicksort

        for size_factor in (1, 5):
            # Make lists to be sorted from two chunks of different kinds.
            sizes = (15, 20)

            all_lists = [self.make_sample_lists(n * size_factor) for n in sizes]
            for chunks in itertools.product(*all_lists):
                orig_keys = sum(chunks, [])
                keys = self.array_factory(orig_keys)
                f(keys)
                # The list is now sorted
                self.assertSorted(orig_keys, keys)

    def test_run_quicksort_lt(self):
        def lt(a, b):
            return a > b

        f = self.make_quicksort(lt=lt).run_quicksort

        for size_factor in (1, 5):
            # Make lists to be sorted from two chunks of different kinds.
            sizes = (15, 20)

            all_lists = [self.make_sample_lists(n * size_factor) for n in sizes]
            for chunks in itertools.product(*all_lists):
                orig_keys = sum(chunks, [])
                keys = self.array_factory(orig_keys)
                f(keys)
                # The list is now rev-sorted
                self.assertSorted(orig_keys, keys[::-1])

        # An imperfect comparison function, as LT(a, b) does not imply not LT(b, a).
        # The sort should handle it gracefully.
        def lt_floats(a, b):
            return math.isnan(b) or a < b

        f = self.make_quicksort(lt=lt_floats).run_quicksort

        np.random.seed(42)
        for size in (5, 20, 50, 500):
            orig = np.random.random(size=size) * 100
            orig[np.random.random(size=size) < 0.1] = float('nan')
            orig_keys = list(orig)
            keys = self.array_factory(orig_keys)
            f(keys)
            non_nans = orig[~np.isnan(orig)]
            # Non-NaNs are sorted at the front
            self.assertSorted(non_nans, keys[:len(non_nans)])


class TestQuicksortPurePython(BaseQuicksortTest, TestCase):

    quicksort = py_quicksort
    make_quicksort = staticmethod(make_py_quicksort)

    # Much faster than a Numpy array in pure Python
    array_factory = list


class TestQuicksortArrays(BaseQuicksortTest, TestCase):

    quicksort = jit_quicksort
    make_quicksort = staticmethod(make_jit_quicksort)

    def array_factory(self, lst):
        return np.array(lst, dtype=np.float64)

class TestQuicksortMultidimensionalArrays(BaseSortingTest, TestCase):

    quicksort = make_jit_quicksort(is_np_array=True)
    make_quicksort = staticmethod(make_jit_quicksort)

    def assertSorted(self, orig, result):
        self.assertEqual(orig.shape, result.shape)
        self.assertPreciseEqual(orig, result)

    def array_factory(self, lst, shape=None):
        array = np.array(lst, dtype=np.float64)
        if shape is None:
            return array.reshape(-1, array.shape[0])
        else:
            return array.reshape(shape)

    def get_shapes(self, n):
        shapes = []
        if n == 1:
            return shapes

        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0:
                shapes.append((n // i, i))
                shapes.append((i, n // i))
                _shapes = self.get_shapes(n // i)
                for _shape in _shapes:
                    shapes.append((i,) + _shape)
                    shapes.append(_shape + (i,))

        return shapes

    def test_run_quicksort(self):
        f = self.quicksort.run_quicksort

        for size_factor in (1, 5):
            # Make lists to be sorted from two chunks of different kinds.
            sizes = (15, 20)

            all_lists = [self.make_sample_lists(n * size_factor) for n in sizes]
            for chunks in itertools.product(*all_lists):
                orig_keys = sum(chunks, [])
                shape_list = self.get_shapes(len(orig_keys))
                shape_list.append(None)
                for shape in shape_list:
                    keys = self.array_factory(orig_keys, shape=shape)
                    keys_copy = self.array_factory(orig_keys, shape=shape)
                    f(keys)
                    keys_copy.sort()
                    # The list is now sorted
                    self.assertSorted(keys_copy, keys)

    def test_run_quicksort_lt(self):
        def lt(a, b):
            return a > b

        f = self.make_quicksort(lt=lt, is_np_array=True).run_quicksort

        for size_factor in (1, 5):
            # Make lists to be sorted from two chunks of different kinds.
            sizes = (15, 20)

            all_lists = [self.make_sample_lists(n * size_factor) for n in sizes]
            for chunks in itertools.product(*all_lists):
                orig_keys = sum(chunks, [])
                shape_list = self.get_shapes(len(orig_keys))
                shape_list.append(None)
                for shape in shape_list:
                    keys = self.array_factory(orig_keys, shape=shape)
                    keys_copy = -self.array_factory(orig_keys, shape=shape)
                    f(keys)
                    # The list is now rev-sorted
                    keys_copy.sort()
                    keys_copy = -keys_copy
                    self.assertSorted(keys_copy, keys)

        # An imperfect comparison function, as LT(a, b) does not imply not LT(b, a).
        # The sort should handle it gracefully.
        def lt_floats(a, b):
            return math.isnan(b) or a < b

        f = self.make_quicksort(lt=lt_floats, is_np_array=True).run_quicksort

        np.random.seed(42)
        for size in (5, 20, 50, 500):
            orig = np.random.random(size=size) * 100
            orig[np.random.random(size=size) < 0.1] = float('nan')
            orig_keys = list(orig)
            shape_list = self.get_shapes(len(orig_keys))
            shape_list.append(None)
            for shape in shape_list:
                keys = self.array_factory(orig_keys, shape=shape)
                keys_copy = self.array_factory(orig_keys, shape=shape)
                f(keys)
                keys_copy.sort()
                # Non-NaNs are sorted at the front
                self.assertSorted(keys_copy, keys)

class TestNumpySort(TestCase):

    def setUp(self):
        np.random.seed(42)

    def int_arrays(self):
        for size in (5, 20, 50, 500):
            yield np.random.randint(99, size=size)

    def float_arrays(self):
        for size in (5, 20, 50, 500):
            yield np.random.random(size=size) * 100
        # Now with NaNs.  Numpy sorts them at the end.
        for size in (5, 20, 50, 500):
            orig = np.random.random(size=size) * 100
            orig[np.random.random(size=size) < 0.1] = float('nan')
            yield orig
        # 90% of values are NaNs.
        for size in (50, 500):
            orig = np.random.random(size=size) * 100
            orig[np.random.random(size=size) < 0.9] = float('nan')
            yield orig

    def has_duplicates(self, arr):
        """
        Whether the array has duplicates.  Takes NaNs into account.
        """
        if np.count_nonzero(np.isnan(arr)) > 1:
            return True
        if np.unique(arr).size < arr.size:
            return True
        return False

    def check_sort_inplace(self, pyfunc, cfunc, val):
        expected = copy.copy(val)
        got = copy.copy(val)
        pyfunc(expected)
        cfunc(got)
        self.assertPreciseEqual(got, expected)

    def check_sort_copy(self, pyfunc, cfunc, val):
        orig = copy.copy(val)
        expected = pyfunc(val)
        got = cfunc(val)
        self.assertPreciseEqual(got, expected)
        # The original wasn't mutated
        self.assertPreciseEqual(val, orig)

    def check_argsort(self, pyfunc, cfunc, val, kwargs=None):
        if kwargs is None:
            kwargs = {}
        orig = copy.copy(val)
        expected = pyfunc(val, **kwargs)
        got = cfunc(val, **kwargs)
        self.assertPreciseEqual(orig[got], np.sort(orig),
                                msg="the array wasn't argsorted")
        # Numba and Numpy results may differ if there are duplicates
        # in the array
        if not self.has_duplicates(orig):
            self.assertPreciseEqual(got, expected)
        # The original wasn't mutated
        self.assertPreciseEqual(val, orig)

    def test_array_sort_int(self):
        pyfunc = sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for orig in self.int_arrays():
            self.check_sort_inplace(pyfunc, cfunc, orig)

    def test_array_sort_float(self):
        pyfunc = sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for orig in self.float_arrays():
            self.check_sort_inplace(pyfunc, cfunc, orig)

    def test_array_sort_complex(self):
        pyfunc = sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for real in self.float_arrays():
            imag = real[::]
            np.random.shuffle(imag)
            orig = np.array([complex(*x) for x in zip(real, imag)])
            self.check_sort_inplace(pyfunc, cfunc, orig)

    def test_np_sort_int(self):
        pyfunc = np_sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for orig in self.int_arrays():
            self.check_sort_copy(pyfunc, cfunc, orig)

    def test_np_sort_float(self):
        pyfunc = np_sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for size in (5, 20, 50, 500):
            orig = np.random.random(size=size) * 100
            orig[np.random.random(size=size) < 0.1] = float('nan')
            self.check_sort_copy(pyfunc, cfunc, orig)

    def test_np_sort_complex(self):
        pyfunc = np_sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for size in (5, 20, 50, 500):
            real = np.random.random(size=size) * 100
            imag = np.random.random(size=size) * 100
            real[np.random.random(size=size) < 0.1] = float('nan')
            imag[np.random.random(size=size) < 0.1] = float('nan')
            orig = np.array([complex(*x) for x in zip(real, imag)])
            self.check_sort_copy(pyfunc, cfunc, orig)

    def test_argsort_int(self):
        def check(pyfunc):
            cfunc = jit(nopython=True)(pyfunc)
            for orig in self.int_arrays():
                self.check_argsort(pyfunc, cfunc, orig)

        check(argsort_usecase)
        check(np_argsort_usecase)

    def test_argsort_kind_int(self):
        def check(pyfunc, is_stable):
            cfunc = jit(nopython=True)(pyfunc)
            for orig in self.int_arrays():
                self.check_argsort(pyfunc, cfunc, orig,
                                   dict(is_stable=is_stable))

        check(argsort_kind_usecase, is_stable=True)
        check(np_argsort_kind_usecase, is_stable=True)
        check(argsort_kind_usecase, is_stable=False)
        check(np_argsort_kind_usecase, is_stable=False)

    def test_argsort_float(self):
        def check(pyfunc):
            cfunc = jit(nopython=True)(pyfunc)
            for orig in self.float_arrays():
                self.check_argsort(pyfunc, cfunc, orig)

        check(argsort_usecase)
        check(np_argsort_usecase)

    def test_argsort_float_supplemental(self):
        def check(pyfunc, is_stable):
            cfunc = jit(nopython=True)(pyfunc)
            for orig in self.float_arrays():
                self.check_argsort(pyfunc, cfunc, orig,
                                   dict(is_stable=is_stable))

        check(argsort_kind_usecase, is_stable=True)
        check(np_argsort_kind_usecase, is_stable=True)
        check(argsort_kind_usecase, is_stable=False)
        check(np_argsort_kind_usecase, is_stable=False)

    def test_argsort_complex(self):
        def check(pyfunc):
            cfunc = jit(nopython=True)(pyfunc)
            for real in self.float_arrays():
                imag = real[::]
                np.random.shuffle(imag)
                orig = np.array([complex(*x) for x in zip(real, imag)])
                self.check_argsort(pyfunc, cfunc, orig)

        check(argsort_usecase)
        check(np_argsort_usecase)

    def test_argsort_complex_supplemental(self):
        def check(pyfunc, is_stable):
            cfunc = jit(nopython=True)(pyfunc)
            for real in self.float_arrays():
                imag = real[::]
                np.random.shuffle(imag)
                orig = np.array([complex(*x) for x in zip(real, imag)])
                self.check_argsort(pyfunc, cfunc, orig,
                                   dict(is_stable=is_stable))

        check(argsort_kind_usecase, is_stable=True)
        check(np_argsort_kind_usecase, is_stable=True)
        check(argsort_kind_usecase, is_stable=False)
        check(np_argsort_kind_usecase, is_stable=False)

    def test_bad_array(self):
        cfunc = jit(nopython=True)(np_sort_usecase)
        msg = '.*Argument "a" must be array-like.*'
        with self.assertRaisesRegex(errors.TypingError, msg) as raises:
            cfunc(None)


class TestPythonSort(TestCase):

    def test_list_sort(self):
        pyfunc = list_sort_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for size in (20, 50, 500):
            orig, ret = cfunc(size)
            self.assertEqual(sorted(orig), ret)
            self.assertNotEqual(orig, ret)   # sanity check

    def test_list_sort_reverse(self):
        pyfunc = list_sort_reverse_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for size in (20, 50, 500):
            for b in (False, True):
                orig, ret = cfunc(size, b)
                self.assertEqual(sorted(orig, reverse=b), ret)
                self.assertNotEqual(orig, ret)   # sanity check

    def test_sorted(self):
        pyfunc = sorted_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for size in (20, 50, 500):
            orig = np.random.random(size=size) * 100
            expected = sorted(orig)
            got = cfunc(orig)
            self.assertPreciseEqual(got, expected)
            self.assertNotEqual(list(orig), got)   # sanity check

    def test_sorted_reverse(self):
        pyfunc = sorted_reverse_usecase
        cfunc = jit(nopython=True)(pyfunc)
        size = 20

        orig = np.random.random(size=size) * 100
        for b in (False, True):
            expected = sorted(orig, reverse=b)
            got = cfunc(orig, b)
            self.assertPreciseEqual(got, expected)
            self.assertNotEqual(list(orig), got)   # sanity check


class TestMergeSort(TestCase):
    def setUp(self):
        np.random.seed(321)

    def check_argsort_stable(self, sorter, low, high, count):
        # make data with high possibility of duplicated key
        data = np.random.randint(low, high, count)
        expect = np.argsort(data, kind='mergesort')
        got = sorter(data)
        np.testing.assert_equal(expect, got)

    def test_argsort_stable(self):
        arglist = [
            (-2, 2, 5),
            (-5, 5, 10),
            (0, 10, 101),
            (0, 100, 1003),
        ]
        imp = make_jit_mergesort(is_argsort=True)
        toplevel = imp.run_mergesort
        sorter = njit(lambda arr: toplevel(arr))
        for args in arglist:
            self.check_argsort_stable(sorter, *args)


nop_compiler = lambda x:x


class TestSortSlashSortedWithKey(MemoryLeakMixin, TestCase):

    def test_01(self):

        a = [3, 1, 4, 1, 5, 9]

        @njit
        def external_key(z):
            return 1. / z

        @njit
        def foo(x, key=None):
            new_x = x[:]
            new_x.sort(key=key)
            return sorted(x[:], key=key), new_x

        self.assertPreciseEqual(foo(a[:]), foo.py_func(a[:]))
        self.assertPreciseEqual(foo(a[:], external_key),
                                foo.py_func(a[:], external_key))

    def test_02(self):

        a = [3, 1, 4, 1, 5, 9]

        @njit
        def foo(x):
            def closure_key(z):
                return 1. / z
            new_x = x[:]
            new_x.sort(key=closure_key)
            return sorted(x[:], key=closure_key), new_x

        self.assertPreciseEqual(foo(a[:]), foo.py_func(a[:]))

    def test_03(self):

        a = [3, 1, 4, 1, 5, 9]

        def gen(compiler):

            @compiler
            def bar(x, func):
                new_x = x[:]
                new_x.sort(key=func)
                return sorted(x[:], key=func), new_x

            @compiler
            def foo(x):
                def closure_escapee_key(z):
                    return 1. / z
                return bar(x, closure_escapee_key)

            return foo

        self.assertPreciseEqual(gen(njit)(a[:]), gen(nop_compiler)(a[:]))

    def test_04(self):

        a = ['a','b','B','b','C','A']

        @njit
        def external_key(z):
            return z.upper()

        @njit
        def foo(x, key=None):
            new_x = x[:]
            new_x.sort(key=key)
            return sorted(x[:], key=key), new_x

        self.assertPreciseEqual(foo(a[:]), foo.py_func(a[:]))
        self.assertPreciseEqual(foo(a[:], external_key),
                                foo.py_func(a[:], external_key))

    def test_05(self):

        a = ['a','b','B','b','C','A']

        @njit
        def external_key(z):
            return z.upper()

        @njit
        def foo(x, key=None, reverse=False):
            new_x = x[:]
            new_x.sort(key=key, reverse=reverse)
            return (sorted(x[:], key=key, reverse=reverse), new_x)

        for key, rev in itertools.product((None, external_key),
                                          (True, False, 1, -12, 0)):
            self.assertPreciseEqual(foo(a[:], key, rev),
                                    foo.py_func(a[:], key, rev))

    def test_optional_on_key(self):
        a = [3, 1, 4, 1, 5, 9]

        @njit
        def foo(x, predicate):
            if predicate:
                def closure_key(z):
                    return 1. / z
            else:
                closure_key = None

            new_x = x[:]
            new_x.sort(key=closure_key)

            return (sorted(x[:], key=closure_key), new_x)

        with self.assertRaises(errors.TypingError) as raises:
            TF = True
            foo(a[:], TF)

        msg = "Key must concretely be None or a Numba JIT compiled function"
        self.assertIn(msg, str(raises.exception))

    def test_exceptions_sorted(self):

        @njit
        def foo_sorted(x, key=None, reverse=False):
            return sorted(x[:], key=key, reverse=reverse)

        @njit
        def foo_sort(x, key=None, reverse=False):
            new_x = x[:]
            new_x.sort(key=key, reverse=reverse)
            return new_x

        @njit
        def external_key(z):
            return 1. / z

        a = [3, 1, 4, 1, 5, 9]

        for impl in (foo_sort, foo_sorted):

            # check illegal key
            with self.assertRaises(errors.TypingError) as raises:
                impl(a, key="illegal")

            expect = "Key must be None or a Numba JIT compiled function"
            self.assertIn(expect, str(raises.exception))

            # check illegal reverse
            with self.assertRaises(errors.TypingError) as raises:
                impl(a, key=external_key, reverse="go backwards")

            expect = "an integer is required for 'reverse'"
            self.assertIn(expect, str(raises.exception))


class TestArrayArgsort(MemoryLeakMixin, TestCase):
    """Tests specific to array.argsort"""

    def test_exceptions(self):

        @njit
        def nonliteral_kind(kind):
            np.arange(5).argsort(kind=kind)

        # check non-literal kind
        with self.assertRaises(errors.TypingError) as raises:
            # valid spelling but not literal
            nonliteral_kind('quicksort')

        expect = '"kind" must be a string literal'
        self.assertIn(expect, str(raises.exception))

        @njit
        def unsupported_kwarg():
            np.arange(5).argsort(foo='')

        with self.assertRaises(errors.TypingError) as raises:
            unsupported_kwarg()

        expect = "Unsupported keywords: ['foo']"
        self.assertIn(expect, str(raises.exception))


if __name__ == '__main__':
    unittest.main()
