from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import TestCase

ndims = 2
size = 10
shape = tuple([size] * ndims)


def _indices_for_nelems(nelems):
    """Returns slices of length nelems, from start onwards, in direction sign."""

    if nelems == 0:
        return [size // 2]  # int index

    res = []
    for step in (1, 2):
        for sign in (-1, 1):
            start = size // 2 - nelems * step * sign // 2
            stop = start + nelems * step * sign
            res.append(slice(start, stop, step * sign))

    return res


def _indices_for_axis():
    """Returns (src, dst) pairs of indices."""

    res = []
    for nelems in (0, 2, 3):
        ind = _indices_for_nelems(nelems)

        # no itertools.product available in Py2.4
        res.extend([(a, b) for a in ind for b in ind])  # all assignments of size "nelems"

    return res


def _indices(ndims):
    """Returns ((axis0_src, axis0_dst), (axis1_src, axis1_dst), ... ) index pairs."""

    ind = _indices_for_axis()

    # no itertools.product available in Py2.4

    res = [[]]
    for i in range(ndims):
        newres = []
        for elem in ind:
            for others in res:
                newres.append([elem] + others)
        res = newres

    return res


def _check_assignment(srcidx, dstidx):
    """Check assignment arr[dstidx] = arr[srcidx] works."""

    arr = np.arange(np.product(shape)).reshape(shape)

    cpy = arr.copy()

    cpy[dstidx] = arr[srcidx]
    arr[dstidx] = arr[srcidx]

    assert np.all(arr == cpy), 'assigning arr[%s] = arr[%s]' % (dstidx, srcidx)


def test_overlapping_assignments():
    """Test automatically generated assignments which overlap in memory."""

    inds = _indices(ndims)

    for ind in inds:
        srcidx = tuple([a[0] for a in ind])
        dstidx = tuple([a[1] for a in ind])

        yield _check_assignment, srcidx, dstidx
