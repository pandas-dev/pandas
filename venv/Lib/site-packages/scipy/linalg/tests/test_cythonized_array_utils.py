import numpy as np
from scipy.linalg import bandwidth, issymmetric, ishermitian
from scipy.conftest import skip_xp_invalid_arg
import pytest
from pytest import raises
from numpy.testing import assert_equal

BANDWIDTH_DTYPES = (
    np.bool,
    np.int8, np.int16, np.int32, np.int64,
    np.uint8, np.uint16, np.uint32, np.uint64,
    np.float32, np.float64,
    np.complex64, np.complex128,
)


@skip_xp_invalid_arg
@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_bandwidth_dtypes():
    n = 5
    for dt in BANDWIDTH_DTYPES:
        _ = bandwidth(np.zeros([n, n], dtype=dt))

    unsupported_dtypes = [np.float16, np.longdouble, np.clongdouble,
                          np.bytes_, np.str_, np.void, np.object_,
                          np.datetime64, np.timedelta64]
    for dt in unsupported_dtypes:
        raises(TypeError, bandwidth, np.zeros([n, n], dtype=dt))


@pytest.mark.parametrize('T', BANDWIDTH_DTYPES)
def test_bandwidth_square_inputs(T):
    n = 20
    k = 4
    R = np.zeros([n, n], dtype=T, order='F')
    # form a banded matrix inplace
    R[[x for x in range(n)], [x for x in range(n)]] = 1
    R[[x for x in range(n-k)], [x for x in range(k, n)]] = 1
    R[[x for x in range(1, n)], [x for x in range(n-1)]] = 1
    R[[x for x in range(k, n)], [x for x in range(n-k)]] = 1
    assert bandwidth(R) == (k, k)
    A = np.array([
        [1, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
    ])
    assert bandwidth(A) == (2, 2)

    A = np.array(
        [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # diagonal
            [[0, 0, 1], [0, 0, 0], [0, 0, 0]],  # upper triangular
            [[0, 0, 0], [0, 0, 0], [1, 0, 0]],  # lower triangular
            [[0, 0, 1], [0, 0, 0], [1, 0, 0]],  # full
            [[0, 0, 1], [0, 0, 0], [0, 1, 0.]],  # upper hessenberg
        ]
    )
    lo, hi = bandwidth(A)
    assert_equal(lo, [0, 0, 2, 2, 1])
    assert_equal(hi, [0, 2, 0, 2, 2])


@skip_xp_invalid_arg
@pytest.mark.parametrize('T', BANDWIDTH_DTYPES)
def test_bandwidth_rect_inputs(T):
    n, m = 10, 20
    k = 5
    R = np.zeros([n, m], dtype=T, order='F')
    # form a banded matrix inplace
    R[[x for x in range(n)], [x for x in range(n)]] = 1
    R[[x for x in range(n-k)], [x for x in range(k, n)]] = 1
    R[[x for x in range(1, n)], [x for x in range(n-1)]] = 1
    R[[x for x in range(k, n)], [x for x in range(n-k)]] = 1
    assert bandwidth(R) == (k, k)

    R2 = np.tril(np.ones((2, 10, 2), dtype=T))
    lo, hi = bandwidth(R2)
    assert_equal(lo, [9, 9])
    assert_equal(hi, [0, 0])

    R3 = np.triu(np.ones((2, 10, 2), dtype=T))
    lo, hi = bandwidth(R3)
    assert_equal(lo, [0, 0])
    assert_equal(hi, [1, 1])


@skip_xp_invalid_arg
@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_issymetric_ishermitian_dtypes():
    n = 5
    for t in np.typecodes['All']:
        A = np.zeros([n, n], dtype=t)
        if t in 'eUVOMm':
            raises(TypeError, issymmetric, A)
            raises(TypeError, ishermitian, A)
        elif t == 'G':  # No-op test. On win these pass on others fail.
            pass
        else:
            assert issymmetric(A)
            assert ishermitian(A)


def test_issymmetric_ishermitian_invalid_input():
    A = np.array([1, 2, 3])
    raises(ValueError, issymmetric, A)
    raises(ValueError, ishermitian, A)
    A = np.array([[[1, 2, 3], [4, 5, 6]]])
    raises(ValueError, issymmetric, A)
    raises(ValueError, ishermitian, A)
    A = np.array([[1, 2, 3], [4, 5, 6]])
    raises(ValueError, issymmetric, A)
    raises(ValueError, ishermitian, A)


def test_issymetric_complex_decimals():
    A = np.arange(1, 10).astype(complex).reshape(3, 3)
    A += np.arange(-4, 5).astype(complex).reshape(3, 3)*1j
    # make entries decimal
    A /= np.pi
    A = A + A.T
    assert issymmetric(A)


def test_ishermitian_complex_decimals():
    A = np.arange(1, 10).astype(complex).reshape(3, 3)
    A += np.arange(-4, 5).astype(complex).reshape(3, 3)*1j
    # make entries decimal
    A /= np.pi
    A = A + A.T.conj()
    assert ishermitian(A)


def test_issymmetric_approximate_results():
    n = 20
    rng = np.random.RandomState(123456789)
    x = rng.uniform(high=5., size=[n, n])
    y = x @ x.T  # symmetric
    p = rng.standard_normal([n, n])
    z = p @ y @ p.T
    assert issymmetric(z, atol=1e-10)
    assert issymmetric(z, atol=1e-10, rtol=0.)
    assert issymmetric(z, atol=0., rtol=1e-12)
    assert issymmetric(z, atol=1e-13, rtol=1e-12)


def test_ishermitian_approximate_results():
    n = 20
    rng = np.random.RandomState(987654321)
    x = rng.uniform(high=5., size=[n, n])
    y = x @ x.T  # symmetric
    p = rng.standard_normal([n, n]) + rng.standard_normal([n, n])*1j
    z = p @ y @ p.conj().T
    assert ishermitian(z, atol=1e-10)
    assert ishermitian(z, atol=1e-10, rtol=0.)
    assert ishermitian(z, atol=0., rtol=1e-12)
    assert ishermitian(z, atol=1e-13, rtol=1e-12)
