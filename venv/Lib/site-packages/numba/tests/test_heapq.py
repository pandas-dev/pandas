import heapq as hq
import itertools

import numpy as np

from numba import jit, typed
from numba.tests.support import TestCase, MemoryLeakMixin


def heapify(x):
    return hq.heapify(x)


def heappop(heap):
    return hq.heappop(heap)


def heappush(heap, item):
    return hq.heappush(heap, item)


def heappushpop(heap, item):
    return hq.heappushpop(heap, item)


def heapreplace(heap, item):
    return hq.heapreplace(heap, item)


def nsmallest(n, iterable):
    return hq.nsmallest(n, iterable)


def nlargest(n, iterable):
    return hq.nlargest(n, iterable)


class _TestHeapq(MemoryLeakMixin):

    def setUp(self):
        super(_TestHeapq, self).setUp()
        self.rnd = np.random.RandomState(42)

    def test_heapify_basic_sanity(self):
        pyfunc = heapify
        cfunc = jit(nopython=True)(pyfunc)

        a = [1, 3, 5, 7, 9, 2, 4, 6, 8, 0]
        b = self.listimpl(a)

        pyfunc(a)
        cfunc(b)
        self.assertPreciseEqual(a, list(b))

        # includes non-finite elements
        element_pool = [3.142, -10.0, 5.5, np.nan, -np.inf, np.inf]

        # list which may contain duplicate elements
        for x in itertools.combinations_with_replacement(element_pool, 6):
            a = list(x)
            b = self.listimpl(a)

            pyfunc(a)
            cfunc(b)
            self.assertPreciseEqual(a, list(b))

        # single element list
        for i in range(len(element_pool)):
            a = [element_pool[i]]
            b = self.listimpl(a)

            pyfunc(a)
            cfunc(b)
            self.assertPreciseEqual(a, list(b))

        # elements are tuples
        a = [(3, 33), (1, 11), (2, 22)]
        b = self.listimpl(a)
        pyfunc(a)
        cfunc(b)
        self.assertPreciseEqual(a, list(b))

    def check_invariant(self, heap):
        for pos, item in enumerate(heap):
            if pos:
                parentpos = (pos - 1) >> 1
                self.assertTrue(heap[parentpos] <= item)

    def test_push_pop(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc_heappush = heappush
        cfunc_heappush = jit(nopython=True)(pyfunc_heappush)

        pyfunc_heappop = heappop
        cfunc_heappop = jit(nopython=True)(pyfunc_heappop)

        heap = self.listimpl([-1.0])
        data = self.listimpl([-1.0])
        self.check_invariant(heap)
        for i in range(256):
            item = self.rnd.randn(1).item(0)
            data.append(item)
            cfunc_heappush(heap, item)
            self.check_invariant(heap)
        results = []
        while heap:
            item = cfunc_heappop(heap)
            self.check_invariant(heap)
            results.append(item)
        data_sorted = data[:]
        data_sorted.sort()
        self.assertPreciseEqual(list(data_sorted), results)
        self.check_invariant(results)

    def test_heapify(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc = heapify
        cfunc = jit(nopython=True)(pyfunc)

        for size in list(range(1, 30)) + [20000]:
            heap = self.listimpl(self.rnd.random_sample(size))
            cfunc(heap)
            self.check_invariant(heap)

    def test_heapify_exceptions(self):
        pyfunc = heapify
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc((1, 5, 4))

        msg = 'heap argument must be a list'
        self.assertIn(msg, str(e.exception))

        with self.assertTypingError() as e:
            cfunc(self.listimpl([1 + 1j, 2 - 3j]))

        msg = ("'<' not supported between instances "
               "of 'complex' and 'complex'")
        self.assertIn(msg, str(e.exception))

    def test_heappop_basic_sanity(self):
        pyfunc = heappop
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield [1, 3, 5, 7, 9, 2, 4, 6, 8, 0]
            yield [(3, 33), (1, 111), (2, 2222)]
            yield np.full(5, fill_value=np.nan).tolist()
            yield np.linspace(-10, -5, 100).tolist()

        for a in a_variations():
            heapify(a)
            b = self.listimpl(a)

            for i in range(len(a)):
                val_py = pyfunc(a)
                val_c = cfunc(b)
                self.assertPreciseEqual(a, list(b))
                self.assertPreciseEqual(val_py, val_c)

    def test_heappop_exceptions(self):
        pyfunc = heappop
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc((1, 5, 4))

        msg = 'heap argument must be a list'
        self.assertIn(msg, str(e.exception))

    def iterables(self):
        yield self.listimpl([1, 3, 5, 7, 9, 2, 4, 6, 8, 0])
        a = np.linspace(-10, 2, 23)
        yield self.listimpl(a)
        yield self.listimpl(a[::-1])
        self.rnd.shuffle(a)
        yield self.listimpl(a)

    def test_heappush_basic(self):
        pyfunc_push = heappush
        cfunc_push = jit(nopython=True)(pyfunc_push)

        pyfunc_pop = heappop
        cfunc_pop = jit(nopython=True)(pyfunc_pop)

        for iterable in self.iterables():
            expected = sorted(iterable)
            heap = self.listimpl([iterable.pop(0)])  # must initialise heap

            for value in iterable:
                cfunc_push(heap, value)

            got = [cfunc_pop(heap) for _ in range(len(heap))]
            self.assertPreciseEqual(expected, got)

    def test_heappush_exceptions(self):
        pyfunc = heappush
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc((1, 5, 4), 6)

        msg = 'heap argument must be a list'
        self.assertIn(msg, str(e.exception))

        with self.assertTypingError() as e:
            cfunc(self.listimpl([1, 5, 4]), 6.0)

        msg = 'heap type must be the same as item type'
        self.assertIn(msg, str(e.exception))

    def test_nsmallest_basic(self):
        pyfunc = nsmallest
        cfunc = jit(nopython=True)(pyfunc)

        for iterable in self.iterables():
            for n in range(-5, len(iterable) + 3):
                expected = pyfunc(1, iterable)
                got = cfunc(1, iterable)
                self.assertPreciseEqual(expected, got)

        # n is boolean
        out = cfunc(False, self.listimpl([3, 2, 1]))
        self.assertPreciseEqual(out, [])

        out = cfunc(True, self.listimpl([3, 2, 1]))
        self.assertPreciseEqual(out, [1])

        # iterable is not a list
        out = cfunc(2, (6, 5, 4, 3, 2, 1))
        self.assertPreciseEqual(out, [1, 2])

        out = cfunc(3, np.arange(6))
        self.assertPreciseEqual(out, [0, 1, 2])

    def test_nlargest_basic(self):
        pyfunc = nlargest
        cfunc = jit(nopython=True)(pyfunc)

        for iterable in self.iterables():
            for n in range(-5, len(iterable) + 3):
                expected = pyfunc(1, iterable)
                got = cfunc(1, iterable)
                self.assertPreciseEqual(expected, got)

        # n is boolean
        out = cfunc(False, self.listimpl([3, 2, 1]))
        self.assertPreciseEqual(out, [])

        out = cfunc(True, self.listimpl([3, 2, 1]))
        self.assertPreciseEqual(out, [3])

        # iterable is not a list
        out = cfunc(2, (6, 5, 4, 3, 2, 1))
        self.assertPreciseEqual(out, [6, 5])

        out = cfunc(3, np.arange(6))
        self.assertPreciseEqual(out, [5, 4, 3])

    def _assert_typing_error(self, cfunc):

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc(2.2, self.listimpl([3, 2, 1]))

        msg = "First argument 'n' must be an integer"
        self.assertIn(msg, str(e.exception))

        with self.assertTypingError() as e:
            cfunc(2, 100)

        msg = "Second argument 'iterable' must be iterable"
        self.assertIn(msg, str(e.exception))

    def test_nsmallest_exceptions(self):
        pyfunc = nsmallest
        cfunc = jit(nopython=True)(pyfunc)
        self._assert_typing_error(cfunc)

    def test_nlargest_exceptions(self):
        pyfunc = nlargest
        cfunc = jit(nopython=True)(pyfunc)
        self._assert_typing_error(cfunc)

    def test_heapreplace_basic(self):
        pyfunc = heapreplace
        cfunc = jit(nopython=True)(pyfunc)

        a = [1, 3, 5, 7, 9, 2, 4, 6, 8, 0]

        heapify(a)
        b = self.listimpl(a)

        for item in [-4, 4, 14]:
            pyfunc(a, item)
            cfunc(b, item)
            self.assertPreciseEqual(a, list(b))

        a = np.linspace(-3, 13, 20)
        a[4] = np.nan
        a[-1] = np.inf
        a = a.tolist()

        heapify(a)
        b = self.listimpl(a)

        for item in [-4.0, 3.142, -np.inf, np.inf]:
            pyfunc(a, item)
            cfunc(b, item)
            self.assertPreciseEqual(a, list(b))

    def test_heapreplace_exceptions(self):
        pyfunc = heapreplace
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc((1, 5, 4), -1)

        msg = 'heap argument must be a list'
        self.assertIn(msg, str(e.exception))

        with self.assertTypingError() as e:
            cfunc(self.listimpl([1, 5, 4]), -1.0)

        msg = 'heap type must be the same as item type'
        self.assertIn(msg, str(e.exception))

    def heapiter(self, heap):
        try:
            while 1:
                yield heappop(heap)
        except IndexError:
            pass

    def test_nbest(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        cfunc_heapify = jit(nopython=True)(heapify)
        cfunc_heapreplace = jit(nopython=True)(heapreplace)

        data = self.rnd.choice(range(2000), 1000).tolist()
        heap = self.listimpl(data[:10])
        cfunc_heapify(heap)

        for item in data[10:]:
            if item > heap[0]:
                cfunc_heapreplace(heap, item)

        self.assertPreciseEqual(list(self.heapiter(list(heap))),
                                sorted(data)[-10:])

    def test_heapsort(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        cfunc_heapify = jit(nopython=True)(heapify)
        cfunc_heappush = jit(nopython=True)(heappush)
        cfunc_heappop = jit(nopython=True)(heappop)

        for trial in range(100):
            # Ensure consistency of typing, use float64 as it's double
            # everywhere
            values = np.arange(5, dtype=np.float64)
            data = self.listimpl(self.rnd.choice(values, 10))
            if trial & 1:
                heap = data[:]
                cfunc_heapify(heap)
            else:
                heap = self.listimpl([data[0]])
                for item in data[1:]:
                    cfunc_heappush(heap, item)
            heap_sorted = [cfunc_heappop(heap) for _ in range(10)]
            self.assertPreciseEqual(heap_sorted, sorted(data))

    def test_nsmallest(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc = nsmallest
        cfunc = jit(nopython=True)(pyfunc)

        data = self.listimpl(self.rnd.choice(range(2000), 1000))

        for n in (0, 1, 2, 10, 100, 400, 999, 1000, 1100):
            self.assertPreciseEqual(list(cfunc(n, data)), sorted(data)[:n])

    def test_nlargest(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc = nlargest
        cfunc = jit(nopython=True)(pyfunc)

        data = self.listimpl(self.rnd.choice(range(2000), 1000))

        for n in (0, 1, 2, 10, 100, 400, 999, 1000, 1100):
            self.assertPreciseEqual(list(cfunc(n, data)),
                                    sorted(data, reverse=True)[:n])

    def test_nbest_with_pushpop(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc_heappushpop = heappushpop
        cfunc_heappushpop = jit(nopython=True)(pyfunc_heappushpop)

        pyfunc_heapify = heapify
        cfunc_heapify = jit(nopython=True)(pyfunc_heapify)

        # Ensure consistency of typing, use float64 as it's double everywhere
        values = np.arange(2000, dtype=np.float64)
        data = self.listimpl(self.rnd.choice(values, 1000))
        heap = data[:10]
        cfunc_heapify(heap)

        for item in data[10:]:
            cfunc_heappushpop(heap, item)

        self.assertPreciseEqual(list(self.heapiter(list(heap))),
                                sorted(data)[-10:])

    def test_heappushpop(self):
        # inspired by
        # https://github.com/python/cpython/blob/e42b7051/Lib/test/test_heapq.py
        pyfunc = heappushpop
        cfunc = jit(nopython=True)(pyfunc)

        h = self.listimpl([1.0])
        x = cfunc(h, 10.0)
        self.assertPreciseEqual((list(h), x), ([10.0], 1.0))
        self.assertPreciseEqual(type(h[0]), float)
        self.assertPreciseEqual(type(x), float)

        h = self.listimpl([10])
        x = cfunc(h, 9)
        self.assertPreciseEqual((list(h), x), ([10], 9))

        h = self.listimpl([10])
        x = cfunc(h, 11)
        self.assertPreciseEqual((list(h), x), ([11], 10))

    def test_heappushpop_exceptions(self):
        pyfunc = heappushpop
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc((1, 5, 4), -1)

        msg = 'heap argument must be a list'
        self.assertIn(msg, str(e.exception))

        with self.assertTypingError() as e:
            cfunc(self.listimpl([1, 5, 4]), False)

        msg = 'heap type must be the same as item type'
        self.assertIn(msg, str(e.exception))


class TestHeapqReflectedList(_TestHeapq, TestCase):
    """Test heapq with reflected lists"""

    listimpl = list


class TestHeapqTypedList(_TestHeapq, TestCase):
    """Test heapq with typed lists"""

    listimpl = typed.List
