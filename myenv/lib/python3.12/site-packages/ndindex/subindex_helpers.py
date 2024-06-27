"""
This file has the main algorithm for Slice.as_subindex(Slice)

Since Integer can use the same algorithm via Slice(i, i+1), and IntegerArray
needs to do this but in a way that only uses array friendly operations, we
need to have this factored out into a separately callable function.
"""

import sys

from ._crt import crt, ilcm

def _crt(m1, m2, v1, v2):
    """
    Chinese Remainder Theorem

    Returns x such that x = v1 (mod m1) and x = v2 (mod m2), or None if no
    such solution exists.

    """
    # Avoid calling crt in the cases where the inputs would be arrays.
    if m1 == 1:
        return v2 % m2
    if m2 == 1:
        return v1 % m1

    assert m1 > 0
    assert m2 > 0
    res = crt([m1, m2], [v1, v2])
    if res is None:
        return res
    return res

def _ilcm(a, b):
    # Avoid calling ilcm in the cases where the inputs would be arrays.
    if a == 1:
        return b
    if b == 1:
        return a

    assert a > 0
    assert b > 0

    return ilcm(a, b)

def where(cond, x, y):
    if 'numpy' in sys.modules:
        from numpy import where
        return where(cond, x, y)
    return x if cond else y # pragma: no cover

def ceiling(a, b):
    """
    Returns ceil(a/b)
    """
    return -(-a//b)

def _max(a, b):
    if isinstance(a, int) and isinstance(b, int):
        return max(a, b)

    from numpy import broadcast_arrays, amax

    return amax(broadcast_arrays(a, b), axis=0)

def _min(a, b):
    if isinstance(a, int) and isinstance(b, int):
        return min(a, b)

    from numpy import broadcast_arrays, amin

    return amin(broadcast_arrays(a, b), axis=0)

def _smallest(x, a, m):
    """
    Gives the smallest integer >= x that equals a (mod m)

    Assumes x >= 0, m >= 1, and 0 <= a < m.
    """
    n = ceiling(x - a, m)
    return a + n*m

def subindex_slice(s_start, s_stop, s_step, i_start, i_stop, i_step):
    """
    Computes s.as_subindex(i) for slices s and i in a way that is (mostly)
    compatible with NumPy arrays.

    Returns (start, stop, step).

    """
    # Chinese Remainder Theorem. We are looking for a solution to
    #
    # x = s.start (mod s.step)
    # x = index.start (mod index.step)
    #
    # If crt() returns None, then there are no solutions (the slices do
    # not overlap).
    common = _crt(s_step, i_step, s_start, i_start)

    if common is None:
        return (0, 0, 1)
    lcm = _ilcm(s_step, i_step)
    start = _max(s_start, i_start)

    # Get the smallest lcm multiple of common that is >= start
    start = _smallest(start, common, lcm)
    # Finally, we need to shift start so that it is relative to index
    start = (start - i_start)//i_step

    stop = ceiling((_min(s_stop, i_stop) - i_start), i_step)
    stop = where(stop < 0, 0, stop)

    step = lcm//i_step # = s_step//igcd(s_step, i_step)

    return (start, stop, step)
