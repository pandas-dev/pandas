"""
Implementation of the Chinese Remainder Theorem.

The code in this file is taken from SymPy (see the license below). We could
(and used to) just import sympy.ntheory.modular.crt, but SymPy is a very heavy
dependency for a rather basic algorithm.

Note: although many of the helper functions from SymPy are copied here, these
functions should not be used by code outside of ndindex. If you want to use
one of these algorithms, you should import them from SymPy (note that the
functions here have had their functionality stripped down from the SymPy
versions).

License
-------

Copyright (c) 2006-2021 SymPy Development Team

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
"""

from math import gcd, prod

def gcdex(a, b):
    """Returns x, y, g such that g = x*a + y*b = gcd(a, b).

    Examples
    ========

    >>> from ndindex._crt import gcdex
    >>> gcdex(2, 3)
    (-1, 1, 1)
    >>> gcdex(10, 12)
    (-1, 1, 2)

    >>> x, y, g = gcdex(100, 2004)
    >>> x, y, g
    (-20, 1, 4)
    >>> x*100 + y*2004
    4

    """
    if (not a) and (not b):
        return (0, 1, 0)

    if not a:
        return (0, b//abs(b), abs(b))
    if not b:
        return (a//abs(a), 0, abs(a))

    if a < 0:
        a, x_sign = -a, -1
    else:
        x_sign = 1

    if b < 0:
        b, y_sign = -b, -1
    else:
        y_sign = 1

    x, y, r, s = 1, 0, 0, 1

    while b:
        (c, q) = (a % b, a // b)
        (a, b, r, s, x, y) = (b, c, x - q*r, y - q*s, r, s)

    return (x*x_sign, y*y_sign, a)

def _crt(U, M):
    """
    Chinese Remainder Theorem.

    Given a set of integer residues ``u_0,...,u_n`` and a set of
    co-prime integer moduli ``m_0,...,m_n``, returns an integer
    ``u``, such that ``u = u_i mod m_i`` for ``i = ``0,...,n``.

    Examples
    ========

    Consider a set of residues ``U = [49, 76, 65]``
    and a set of moduli ``M = [99, 97, 95]``. Then we have::

       >>> from ndindex._crt import _crt

       >>> _crt([49, 76, 65], [99, 97, 95])
       639985

    This is the correct result because::

       >>> [639985 % m for m in [99, 97, 95]]
       [49, 76, 65]

    Note: this is a low-level routine with no error checking.
    """
    p = prod(M)
    v = 0

    for u, m in zip(U, M):
        e = p // m
        s, _, _ = gcdex(e, m)
        v += e*(u*s % m)

    return v % p


def solve_congruence(*remainder_modulus_pairs):
    """Compute the integer ``n`` that has the residual ``ai`` when it is
    divided by ``mi`` where the ``ai`` and ``mi`` are given as pairs to
    this function: ((a1, m1), (a2, m2), ...). If there is no solution,
    return None. Otherwise return ``n`` and its modulus.

    The ``mi`` values need not be co-prime.

    Examples
    ========

    >>> from ndindex._crt import solve_congruence

    What number is 2 mod 3, 3 mod 5 and 2 mod 7?

    >>> solve_congruence((2, 3), (3, 5), (2, 7))
    23
    >>> [23 % m for m in [3, 5, 7]]
    [2, 3, 2]

    If you prefer to work with all remainder in one list and
    all moduli in another, send the arguments like this:

    >>> solve_congruence(*zip((2, 3, 2), (3, 5, 7)))
    23

    The moduli need not be co-prime; in this case there may or
    may not be a solution:

    >>> solve_congruence((2, 3), (4, 6)) is None
    True

    >>> solve_congruence((2, 3), (5, 6))
    5

    """
    def combine(c1, c2):
        """Return the tuple (a, m) which satisfies the requirement
        that n = a + i*m satisfy n = a1 + j*m1 and n = a2 = k*m2.

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Method_of_successive_substitution
        """
        a1, m1 = c1
        a2, m2 = c2
        a, b, c = m1, a2 - a1, m2
        g = gcd(a, b, c)
        a, b, c = [i//g for i in [a, b, c]]
        if a != 1:
            inv_a, _, g = gcdex(a, c)
            if g != 1:
                return None
            b *= inv_a
        a, m = a1 + m1*b, m1*c
        return a, m

    rm = remainder_modulus_pairs

    rv = (0, 1)
    for rmi in rm:
        rv = combine(rv, rmi)
        if rv is None:
            break
        n, m = rv
        n = n % m
    else:
        return n

def crt(m, v, check=True):
    r"""Chinese Remainder Theorem.

    The moduli in m are assumed to be pairwise coprime.  The output
    is then an integer f, such that f = v_i mod m_i for each pair out
    of v and m.

    If the moduli are not co-prime the correct result will be returned
    if/when the test of the result is found to be incorrect. This result
    will be None if there is no solution.

    The keyword ``check`` can be set to False if it is known that the moduli
    are coprime.

    Examples
    ========

    As an example consider a set of residues ``U = [49, 76, 65]``
    and a set of moduli ``M = [99, 97, 95]``. Then we have::

       >>> from ndindex._crt import crt

       >>> crt([99, 97, 95], [49, 76, 65])
       639985

    This is the correct result because::

       >>> [639985 % m for m in [99, 97, 95]]
       [49, 76, 65]

    If the moduli are not co-prime, you may receive an incorrect result
    if you use ``check=False``:

       >>> crt([12, 6, 17], [3, 4, 2], check=False)
       954
       >>> [954 % m for m in [12, 6, 17]]
       [6, 0, 2]
       >>> crt([12, 6, 17], [3, 4, 2]) is None
       True
       >>> crt([3, 6], [2, 5])
       5

    Note: the order of gf_crt's arguments is reversed relative to crt,
    and that solve_congruence takes residue, modulus pairs.

    Programmer's note: rather than checking that all pairs of moduli share
    no GCD (an O(n**2) test) and rather than factoring all moduli and seeing
    that there is no factor in common, a check that the result gives the
    indicated residuals is performed -- an O(n) operation.
    """
    result = _crt(v, m)

    if check:
        if not all(v % m == result % m for v, m in zip(v, m)):
            result = solve_congruence(*list(zip(v, m)))

    return result

def ilcm(a, b):
    """Computes integer least common multiple.

    Examples
    ========

    >>> from ndindex._crt import ilcm
    >>> ilcm(5, 10)
    10
    >>> ilcm(7, 3)
    21

    """
    if 0 in [a, b]:
        return 0
    return a // gcd(a, b) * b # since gcd(a,b) | a
