from numpy.testing import assert_allclose
from scipy.linalg import cython_lapack as cython_lapack
from scipy.linalg import lapack


class TestLamch:

    def test_slamch(self):
        lapack_slamch = lapack.get_lapack_funcs(
            'lamch', dtype='float32', ilp64='preferred'
        )
        for c in [b'e', b's', b'b', b'p', b'n', b'r', b'm', b'u', b'l', b'o']:
            assert_allclose(cython_lapack._test_slamch(c),
                            lapack_slamch(c))

    def test_dlamch(self):
        lapack_dlamch = lapack.get_lapack_funcs(
            'lamch', dtype='float64', ilp64='preferred'
        )

        for c in [b'e', b's', b'b', b'p', b'n', b'r', b'm', b'u', b'l', b'o']:
            assert_allclose(cython_lapack._test_dlamch(c),
                            lapack_dlamch(c))

    def test_complex_ladiv(self):
        cx = .5 + 1.j
        cy = .875 + 2.j
        assert_allclose(cython_lapack._test_zladiv(cy, cx), 1.95+0.1j)
        assert_allclose(cython_lapack._test_cladiv(cy, cx), 1.95+0.1j)


class TestBlasInt:
    """Test for blas_int/blas_bint types used in cython_blas/cython_lapack.

    Note similar test is in test_blas.py: both cython_blas.pyx and cython_lapack.pyx
    define the types, so we smoke test that `blas_int` and `blas_bint` types
    can be cimported from both `cython_blas` and `cython_lapack` extensions.
    """

    def test_blas_int_size(self):
        """Verify blas_int and blas_bint sizes matches the build configuration."""
        from scipy.__config__ import CONFIG
        size = cython_lapack._blas_int_size()
        cython_blas_ilp64 = CONFIG['Build Dependencies']['blas']['cython blas ilp64']
        if cython_blas_ilp64:
            assert size == 8, f"ILP64 build but blas_int is {size} bytes"
        else:
            assert size == 4, f"LP64 build but blas_int is {size} bytes"

        # also verify blas_bint size
        # bint_sizes is a 2-tuple, (sizeof(blas_bint), sizeof(bint))
        bint_sizes = cython_lapack._blas_bint_bint_sizes()
        if cython_blas_ilp64:
            assert bint_sizes == (8, 4)
        else:
            assert bint_sizes == (4, 4)
