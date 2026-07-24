from math import gcd

from .._crt import crt, ilcm, gcdex

from hypothesis import given, example
from hypothesis.strategies import integers, lists, shared

from sympy.ntheory.modular import crt as crt_sympy
size = shared(integers(min_value=1, max_value=10))
@given(
    size.flatmap(lambda s: lists(integers(min_value=1), min_size=s, max_size=s)),
    size.flatmap(lambda s: lists(integers(), min_size=s, max_size=s)),
)
def test_crt(m, v):
    res = crt(m, v)

    if res is not None:
        for m_i, v_i in zip(m, v):
            assert v_i % m_i == res % m_i

        assert res == crt_sympy(m, v)[0]
    else:
        assert crt_sympy(m, v) is None

@example(1, 2)
@given(integers(min_value=0), integers(min_value=0))
def test_ilcm(x, y):
    L = ilcm(x, y)

    if 0 in [x, y]:
        assert L == 0
        return

    assert L >= x
    assert L >= y

    # L is a common multiple
    assert L % x == 0
    assert L % y == 0

    if L - min(x, y) <= 1000:
        # L is the least common multiple
        for i in range(min(x, y), L):
            assert i % x != 0 or i % y != 0

@example(0, 3)
@example(3, 0)
@example(0, 0)
@given(integers(), integers())
def test_gcdex(a, b):
    x, y, g = gcdex(a, b)

    assert g == gcd(a, b)
    assert x*a + y*b == g
