import numpy as np

from numba import cuda, complex64, int32, float64
from numba.cuda.testing import unittest, CUDATestCase
from numba.core.config import ENABLE_CUDASIM

CONST_EMPTY = np.array([])
CONST1D = np.arange(10, dtype=np.float64) / 2.
CONST2D = np.asfortranarray(
    np.arange(100, dtype=np.int32).reshape(10, 10))
CONST3D = ((np.arange(5 * 5 * 5, dtype=np.complex64).reshape(5, 5, 5) + 1j) /
           2j)
CONST3BYTES = np.arange(3, dtype=np.uint8)

CONST_RECORD_EMPTY = np.array(
    [],
    dtype=[('x', float), ('y', int)])
CONST_RECORD = np.array(
    [(1.0, 2), (3.0, 4)],
    dtype=[('x', float), ('y', int)])
CONST_RECORD_ALIGN = np.array(
    [(1, 2, 3, 0xDEADBEEF, 8), (4, 5, 6, 0xBEEFDEAD, 10)],
    dtype=np.dtype(
        dtype=[
            ('a', np.uint8),
            ('b', np.uint8),
            ('x', np.uint8),
            ('y', np.uint32),
            ('z', np.uint8),
        ],
        align=True))


def cuconstEmpty(A):
    C = cuda.const.array_like(CONST_EMPTY)
    i = cuda.grid(1)
    A[i] = len(C)


def cuconst(A):
    C = cuda.const.array_like(CONST1D)
    i = cuda.grid(1)

    # +1 or it'll be loaded & stored as a u32
    A[i] = C[i] + 1.0


def cuconst2d(A):
    C = cuda.const.array_like(CONST2D)
    i, j = cuda.grid(2)
    A[i, j] = C[i, j]


def cuconst3d(A):
    C = cuda.const.array_like(CONST3D)
    i = cuda.threadIdx.x
    j = cuda.threadIdx.y
    k = cuda.threadIdx.z
    A[i, j, k] = C[i, j, k]


def cuconstRecEmpty(A):
    C = cuda.const.array_like(CONST_RECORD_EMPTY)
    i = cuda.grid(1)
    A[i] = len(C)


def cuconstRec(A, B):
    C = cuda.const.array_like(CONST_RECORD)
    i = cuda.grid(1)
    A[i] = C[i]['x']
    B[i] = C[i]['y']


def cuconstRecAlign(A, B, C, D, E):
    Z = cuda.const.array_like(CONST_RECORD_ALIGN)
    i = cuda.grid(1)
    A[i] = Z[i]['a']
    B[i] = Z[i]['b']
    C[i] = Z[i]['x']
    D[i] = Z[i]['y']
    E[i] = Z[i]['z']


def cuconstAlign(z):
    a = cuda.const.array_like(CONST3BYTES)
    b = cuda.const.array_like(CONST1D)
    i = cuda.grid(1)
    z[i] = a[i] + b[i]


class TestCudaConstantMemory(CUDATestCase):
    def test_const_array(self):
        sig = (float64[:],)
        jcuconst = cuda.jit(sig)(cuconst)
        A = np.zeros_like(CONST1D)
        jcuconst[2, 5](A)
        self.assertTrue(np.all(A == CONST1D + 1))

        if not ENABLE_CUDASIM:
            self.assertIn(
                'ld.const.f64',
                jcuconst.inspect_asm(sig),
                "as we're adding to it, load as a double")

    def test_const_empty(self):
        jcuconstEmpty = cuda.jit('void(int64[:])')(cuconstEmpty)
        A = np.full(1, fill_value=-1, dtype=np.int64)
        jcuconstEmpty[1, 1](A)
        self.assertTrue(np.all(A == 0))

    def test_const_align(self):
        jcuconstAlign = cuda.jit('void(float64[:])')(cuconstAlign)
        A = np.full(3, fill_value=np.nan, dtype=float)
        jcuconstAlign[1, 3](A)
        self.assertTrue(np.all(A == (CONST3BYTES + CONST1D[:3])))

    def test_const_array_2d(self):
        sig = (int32[:,:],)
        jcuconst2d = cuda.jit(sig)(cuconst2d)
        A = np.zeros_like(CONST2D, order='C')
        jcuconst2d[(2, 2), (5, 5)](A)
        self.assertTrue(np.all(A == CONST2D))

        if not ENABLE_CUDASIM:
            self.assertIn(
                'ld.const.u32',
                jcuconst2d.inspect_asm(sig),
                "load the ints as ints")

    def test_const_array_3d(self):
        sig = (complex64[:,:,:],)
        jcuconst3d = cuda.jit(sig)(cuconst3d)
        A = np.zeros_like(CONST3D, order='F')
        jcuconst3d[1, (5, 5, 5)](A)
        self.assertTrue(np.all(A == CONST3D))

        if not ENABLE_CUDASIM:
            asm = jcuconst3d.inspect_asm(sig)
            complex_load = 'ld.const.v2.f32'
            description = 'Load the complex as a vector of 2x f32'
            self.assertIn(complex_load, asm, description)

    def test_const_record_empty(self):
        jcuconstRecEmpty = cuda.jit('void(int64[:])')(cuconstRecEmpty)
        A = np.full(1, fill_value=-1, dtype=np.int64)
        jcuconstRecEmpty[1, 1](A)
        self.assertTrue(np.all(A == 0))

    def test_const_record(self):
        A = np.zeros(2, dtype=float)
        B = np.zeros(2, dtype=int)
        jcuconst = cuda.jit(cuconstRec).specialize(A, B)

        jcuconst[2, 1](A, B)
        np.testing.assert_allclose(A, CONST_RECORD['x'])
        np.testing.assert_allclose(B, CONST_RECORD['y'])

    def test_const_record_align(self):
        A = np.zeros(2, dtype=np.float64)
        B = np.zeros(2, dtype=np.float64)
        C = np.zeros(2, dtype=np.float64)
        D = np.zeros(2, dtype=np.float64)
        E = np.zeros(2, dtype=np.float64)
        jcuconst = cuda.jit(cuconstRecAlign).specialize(A, B, C, D, E)

        jcuconst[2, 1](A, B, C, D, E)
        np.testing.assert_allclose(A, CONST_RECORD_ALIGN['a'])
        np.testing.assert_allclose(B, CONST_RECORD_ALIGN['b'])
        np.testing.assert_allclose(C, CONST_RECORD_ALIGN['x'])
        np.testing.assert_allclose(D, CONST_RECORD_ALIGN['y'])
        np.testing.assert_allclose(E, CONST_RECORD_ALIGN['z'])


if __name__ == '__main__':
    unittest.main()
