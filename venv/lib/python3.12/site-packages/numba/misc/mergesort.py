"""
The same algorithm as translated from numpy.
See numpy/core/src/npysort/mergesort.c.src.
The high-level numba code is adding a little overhead comparing to
the pure-C implementation in numpy.
"""
import numpy as np
from collections import namedtuple

# Array size smaller than this will be sorted by insertion sort
SMALL_MERGESORT = 20


MergesortImplementation = namedtuple('MergesortImplementation', [
    'run_mergesort',
])


def make_mergesort_impl(wrap, lt=None, is_argsort=False):
    kwargs_lite = dict(no_cpython_wrapper=True, _nrt=False)

    # The less than
    if lt is None:
        @wrap(**kwargs_lite)
        def lt(a, b):
            return a < b
    else:
        lt = wrap(**kwargs_lite)(lt)

    if is_argsort:
        @wrap(**kwargs_lite)
        def lessthan(a, b, vals):
            return lt(vals[a], vals[b])
    else:
        @wrap(**kwargs_lite)
        def lessthan(a, b, vals):
            return lt(a, b)

    @wrap(**kwargs_lite)
    def argmergesort_inner(arr, vals, ws):
        """The actual mergesort function

        Parameters
        ----------
        arr : array [read+write]
            The values being sorted inplace.  For argsort, this is the
            indices.
        vals : array [readonly]
            ``None`` for normal sort.  In argsort, this is the actual array values.
        ws : array [write]
            The workspace.  Must be of size ``arr.size // 2``
        """
        if arr.size > SMALL_MERGESORT:
            # Merge sort
            mid = arr.size // 2

            argmergesort_inner(arr[:mid], vals, ws)
            argmergesort_inner(arr[mid:], vals, ws)

            # Copy left half into workspace so we don't overwrite it
            for i in range(mid):
                ws[i] = arr[i]

            # Merge
            left = ws[:mid]
            right = arr[mid:]
            out = arr

            i = j = k = 0
            while i < left.size and j < right.size:
                if not lessthan(right[j], left[i], vals):
                    out[k] = left[i]
                    i += 1
                else:
                    out[k] = right[j]
                    j += 1
                k += 1

            # Leftovers
            while i < left.size:
                out[k] = left[i]
                i += 1
                k += 1

            while j < right.size:
                out[k] = right[j]
                j += 1
                k += 1
        else:
            # Insertion sort
            i = 1
            while i < arr.size:
                j = i
                while j > 0 and lessthan(arr[j], arr[j - 1], vals):
                    arr[j - 1], arr[j] = arr[j], arr[j - 1]
                    j -= 1
                i += 1

    # The top-level entry points

    @wrap(no_cpython_wrapper=True)
    def mergesort(arr):
        "Inplace"
        ws = np.empty(arr.size // 2, dtype=arr.dtype)
        argmergesort_inner(arr, None, ws)
        return arr


    @wrap(no_cpython_wrapper=True)
    def argmergesort(arr):
        "Out-of-place"
        idxs = np.arange(arr.size)
        ws = np.empty(arr.size // 2, dtype=idxs.dtype)
        argmergesort_inner(idxs, arr, ws)
        return idxs

    return MergesortImplementation(
        run_mergesort=(argmergesort if is_argsort else mergesort)
        )


def make_jit_mergesort(*args, **kwargs):
    from numba import njit
    # NOTE: wrap with njit to allow recursion
    #       because @register_jitable => @overload doesn't support recursion
    return make_mergesort_impl(njit, *args, **kwargs)
