import collections

import numpy as np

from numba.core import types, config


QuicksortImplementation = collections.namedtuple(
    'QuicksortImplementation',
    (# The compile function itself
     'compile',
     # All subroutines exercised by test_sort
     'partition', 'partition3', 'insertion_sort',
     # The top-level function
     'run_quicksort',
     ))


Partition = collections.namedtuple('Partition', ('start', 'stop'))

# Under this size, switch to a simple insertion sort
SMALL_QUICKSORT = 15

MAX_STACK = 100


def make_quicksort_impl(wrap, lt=None, is_argsort=False, is_list=False, is_np_array=False):

    if config.USE_LEGACY_TYPE_SYSTEM:
        intp = types.intp
    else:
        intp = types.py_int
    zero = intp(0)

    # Two subroutines to make the core algorithm generic wrt. argsort
    # or normal sorting.  Note the genericity may make basic sort()
    # slightly slower (~5%)
    if is_argsort:
        if is_list:
            @wrap
            def make_res(A):
                return [x for x in range(len(A))]
        else:
            @wrap
            def make_res(A):
                return np.arange(A.size)

        @wrap
        def GET(A, idx_or_val):
            return A[idx_or_val]

    else:
        @wrap
        def make_res(A):
            return A

        @wrap
        def GET(A, idx_or_val):
            return idx_or_val

    def default_lt(a, b):
        """
        Trivial comparison function between two keys.
        """
        return a < b

    LT = wrap(lt if lt is not None else default_lt)

    @wrap
    def insertion_sort(A, R, low, high):
        """
        Insertion sort A[low:high + 1]. Note the inclusive bounds.
        """
        assert low >= 0
        if high <= low:
            return

        for i in range(low + 1, high + 1):
            k = R[i]
            v = GET(A, k)
            # Insert v into A[low:i]
            j = i
            while j > low and LT(v, GET(A, R[j - 1])):
                # Make place for moving A[i] downwards
                R[j] = R[j - 1]
                j -= 1
            R[j] = k

    @wrap
    def partition(A, R, low, high):
        """
        Partition A[low:high + 1] around a chosen pivot.  The pivot's index
        is returned.
        """
        assert low >= 0
        assert high > low

        mid = (low + high) >> 1
        # NOTE: the pattern of swaps below for the pivot choice and the
        # partitioning gives good results (i.e. regular O(n log n))
        # on sorted, reverse-sorted, and uniform arrays.  Subtle changes
        # risk breaking this property.

        # median of three {low, middle, high}
        if LT(GET(A, R[mid]), GET(A, R[low])):
            R[low], R[mid] = R[mid], R[low]
        if LT(GET(A, R[high]), GET(A, R[mid])):
            R[high], R[mid] = R[mid], R[high]
        if LT(GET(A, R[mid]), GET(A, R[low])):
            R[low], R[mid] = R[mid], R[low]
        pivot = GET(A, R[mid])

        # Temporarily stash the pivot at the end
        R[high], R[mid] = R[mid], R[high]
        i = low
        j = high - 1
        while True:
            while i < high and LT(GET(A, R[i]), pivot):
                i += 1
            while j >= low and LT(pivot, GET(A, R[j])):
                j -= 1
            if i >= j:
                break
            R[i], R[j] = R[j], R[i]
            i += 1
            j -= 1
        # Put the pivot back in its final place (all items before `i`
        # are smaller than the pivot, all items at/after `i` are larger)
        R[i], R[high] = R[high], R[i]
        return i

    @wrap
    def partition3(A, low, high):
        """
        Three-way partition [low, high) around a chosen pivot.
        A tuple (lt, gt) is returned such that:
            - all elements in [low, lt) are < pivot
            - all elements in [lt, gt] are == pivot
            - all elements in (gt, high] are > pivot
        """
        mid = (low + high) >> 1
        # median of three {low, middle, high}
        if LT(A[mid], A[low]):
            A[low], A[mid] = A[mid], A[low]
        if LT(A[high], A[mid]):
            A[high], A[mid] = A[mid], A[high]
        if LT(A[mid], A[low]):
            A[low], A[mid] = A[mid], A[low]
        pivot = A[mid]

        A[low], A[mid] = A[mid], A[low]
        lt = low
        gt = high
        i = low + 1
        while i <= gt:
            if LT(A[i], pivot):
                A[lt], A[i] = A[i], A[lt]
                lt += 1
                i += 1
            elif LT(pivot, A[i]):
                A[gt], A[i] = A[i], A[gt]
                gt -= 1
            else:
                i += 1
        return lt, gt

    @wrap
    def run_quicksort1(A):
        R = make_res(A)

        if len(A) < 2:
            return R

        stack = [Partition(zero, zero)] * MAX_STACK
        stack[0] = Partition(zero, len(A) - 1)
        n = 1

        while n > 0:
            n -= 1
            low, high = stack[n]
            # Partition until it becomes more efficient to do an insertion sort
            while high - low >= SMALL_QUICKSORT:
                assert n < MAX_STACK
                i = partition(A, R, low, high)
                # Push largest partition on the stack
                if high - i > i - low:
                    # Right is larger
                    if high > i:
                        stack[n] = Partition(i + 1, high)
                        n += 1
                    high = i - 1
                else:
                    if i > low:
                        stack[n] = Partition(low, i - 1)
                        n += 1
                    low = i + 1

            insertion_sort(A, R, low, high)

        return R

    if is_np_array:
        @wrap
        def run_quicksort(A):
            if A.ndim == 1:
                return run_quicksort1(A)
            else:
                for idx in np.ndindex(A.shape[:-1]):
                    run_quicksort1(A[idx])
                return A
    else:
        @wrap
        def run_quicksort(A):
            return run_quicksort1(A)

    # Unused quicksort implementation based on 3-way partitioning; the
    # partitioning scheme turns out exhibiting bad behaviour on sorted arrays.
    @wrap
    def _run_quicksort(A):
        stack = [Partition(zero, zero)] * 100
        stack[0] = Partition(zero, len(A) - 1)
        n = 1

        while n > 0:
            n -= 1
            low, high = stack[n]
            # Partition until it becomes more efficient to do an insertion sort
            while high - low >= SMALL_QUICKSORT:
                assert n < MAX_STACK
                l, r = partition3(A, low, high)
                # One trivial (empty) partition => iterate on the other
                if r == high:
                    high = l - 1
                elif l == low:
                    low = r + 1
                # Push largest partition on the stack
                elif high - r > l - low:
                    # Right is larger
                    stack[n] = Partition(r + 1, high)
                    n += 1
                    high = l - 1
                else:
                    stack[n] = Partition(low, l - 1)
                    n += 1
                    low = r + 1

            insertion_sort(A, low, high)


    return QuicksortImplementation(wrap,
                                   partition, partition3, insertion_sort,
                                   run_quicksort)


def make_py_quicksort(*args, **kwargs):
    return make_quicksort_impl((lambda f: f), *args, **kwargs)

def make_jit_quicksort(*args, **kwargs):
    from numba.core.extending import register_jitable
    return make_quicksort_impl((lambda f: register_jitable(f)),
                               *args, **kwargs)
