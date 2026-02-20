import numpy as np
from numba import cuda, int32, int64, float32, float64
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim
from numba.core import config


def useful_syncwarp(ary):
    i = cuda.grid(1)
    if i == 0:
        ary[0] = 42
    cuda.syncwarp(0xffffffff)
    ary[i] = ary[0]


def use_shfl_sync_idx(ary, idx):
    i = cuda.grid(1)
    val = cuda.shfl_sync(0xffffffff, i, idx)
    ary[i] = val


def use_shfl_sync_up(ary, delta):
    i = cuda.grid(1)
    val = cuda.shfl_up_sync(0xffffffff, i, delta)
    ary[i] = val


def use_shfl_sync_down(ary, delta):
    i = cuda.grid(1)
    val = cuda.shfl_down_sync(0xffffffff, i, delta)
    ary[i] = val


def use_shfl_sync_xor(ary, xor):
    i = cuda.grid(1)
    val = cuda.shfl_xor_sync(0xffffffff, i, xor)
    ary[i] = val


def use_shfl_sync_with_val(ary, into):
    i = cuda.grid(1)
    val = cuda.shfl_sync(0xffffffff, into, 0)
    ary[i] = val


def use_vote_sync_all(ary_in, ary_out):
    i = cuda.grid(1)
    pred = cuda.all_sync(0xffffffff, ary_in[i])
    ary_out[i] = pred


def use_vote_sync_any(ary_in, ary_out):
    i = cuda.grid(1)
    pred = cuda.any_sync(0xffffffff, ary_in[i])
    ary_out[i] = pred


def use_vote_sync_eq(ary_in, ary_out):
    i = cuda.grid(1)
    pred = cuda.eq_sync(0xffffffff, ary_in[i])
    ary_out[i] = pred


def use_vote_sync_ballot(ary):
    i = cuda.threadIdx.x
    ballot = cuda.ballot_sync(0xffffffff, True)
    ary[i] = ballot


def use_match_any_sync(ary_in, ary_out):
    i = cuda.grid(1)
    ballot = cuda.match_any_sync(0xffffffff, ary_in[i])
    ary_out[i] = ballot


def use_match_all_sync(ary_in, ary_out):
    i = cuda.grid(1)
    ballot, pred = cuda.match_all_sync(0xffffffff, ary_in[i])
    ary_out[i] = ballot if pred else 0


def use_independent_scheduling(arr):
    i = cuda.threadIdx.x
    if i % 4 == 0:
        ballot = cuda.ballot_sync(0x11111111, True)
    elif i % 4 == 1:
        ballot = cuda.ballot_sync(0x22222222, True)
    elif i % 4 == 2:
        ballot = cuda.ballot_sync(0x44444444, True)
    elif i % 4 == 3:
        ballot = cuda.ballot_sync(0x88888888, True)
    arr[i] = ballot


def _safe_cc_check(cc):
    if config.ENABLE_CUDASIM:
        return True
    else:
        return cuda.get_current_device().compute_capability >= cc


@skip_on_cudasim("Warp Operations are not yet implemented on cudasim")
class TestCudaWarpOperations(CUDATestCase):
    def test_useful_syncwarp(self):
        compiled = cuda.jit("void(int32[:])")(useful_syncwarp)
        nelem = 32
        ary = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary)
        self.assertTrue(np.all(ary == 42))

    def test_shfl_sync_idx(self):
        compiled = cuda.jit("void(int32[:], int32)")(use_shfl_sync_idx)
        nelem = 32
        idx = 4
        ary = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary, idx)
        self.assertTrue(np.all(ary == idx))

    def test_shfl_sync_up(self):
        compiled = cuda.jit("void(int32[:], int32)")(use_shfl_sync_up)
        nelem = 32
        delta = 4
        ary = np.empty(nelem, dtype=np.int32)
        exp = np.arange(nelem, dtype=np.int32)
        exp[delta:] -= delta
        compiled[1, nelem](ary, delta)
        self.assertTrue(np.all(ary == exp))

    def test_shfl_sync_down(self):
        compiled = cuda.jit("void(int32[:], int32)")(use_shfl_sync_down)
        nelem = 32
        delta = 4
        ary = np.empty(nelem, dtype=np.int32)
        exp = np.arange(nelem, dtype=np.int32)
        exp[:-delta] += delta
        compiled[1, nelem](ary, delta)
        self.assertTrue(np.all(ary == exp))

    def test_shfl_sync_xor(self):
        compiled = cuda.jit("void(int32[:], int32)")(use_shfl_sync_xor)
        nelem = 32
        xor = 16
        ary = np.empty(nelem, dtype=np.int32)
        exp = np.arange(nelem, dtype=np.int32) ^ xor
        compiled[1, nelem](ary, xor)
        self.assertTrue(np.all(ary == exp))

    def test_shfl_sync_types(self):
        types = int32, int64, float32, float64
        values = (np.int32(-1), np.int64(1 << 42),
                  np.float32(np.pi), np.float64(np.pi))
        for typ, val in zip(types, values):
            compiled = cuda.jit((typ[:], typ))(use_shfl_sync_with_val)
            nelem = 32
            ary = np.empty(nelem, dtype=val.dtype)
            compiled[1, nelem](ary, val)
            self.assertTrue(np.all(ary == val))

    def test_vote_sync_all(self):
        compiled = cuda.jit("void(int32[:], int32[:])")(use_vote_sync_all)
        nelem = 32
        ary_in = np.ones(nelem, dtype=np.int32)
        ary_out = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))
        ary_in[-1] = 0
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))

    def test_vote_sync_any(self):
        compiled = cuda.jit("void(int32[:], int32[:])")(use_vote_sync_any)
        nelem = 32
        ary_in = np.zeros(nelem, dtype=np.int32)
        ary_out = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))
        ary_in[2] = 1
        ary_in[5] = 1
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))

    def test_vote_sync_eq(self):
        compiled = cuda.jit("void(int32[:], int32[:])")(use_vote_sync_eq)
        nelem = 32
        ary_in = np.zeros(nelem, dtype=np.int32)
        ary_out = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))
        ary_in[1] = 1
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))
        ary_in[:] = 1
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))

    def test_vote_sync_ballot(self):
        compiled = cuda.jit("void(uint32[:])")(use_vote_sync_ballot)
        nelem = 32
        ary = np.empty(nelem, dtype=np.uint32)
        compiled[1, nelem](ary)
        self.assertTrue(np.all(ary == np.uint32(0xffffffff)))

    @unittest.skipUnless(_safe_cc_check((7, 0)),
                         "Matching requires at least Volta Architecture")
    def test_match_any_sync(self):
        compiled = cuda.jit("void(int32[:], int32[:])")(use_match_any_sync)
        nelem = 10
        ary_in = np.arange(nelem, dtype=np.int32) % 2
        ary_out = np.empty(nelem, dtype=np.int32)
        exp = np.tile((0b0101010101, 0b1010101010), 5)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == exp))

    @unittest.skipUnless(_safe_cc_check((7, 0)),
                         "Matching requires at least Volta Architecture")
    def test_match_all_sync(self):
        compiled = cuda.jit("void(int32[:], int32[:])")(use_match_all_sync)
        nelem = 10
        ary_in = np.zeros(nelem, dtype=np.int32)
        ary_out = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0b1111111111))
        ary_in[1] = 4
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))

    @unittest.skipUnless(_safe_cc_check((7, 0)),
                         "Independent scheduling requires at least Volta "
                         "Architecture")
    def test_independent_scheduling(self):
        compiled = cuda.jit("void(uint32[:])")(use_independent_scheduling)
        arr = np.empty(32, dtype=np.uint32)
        exp = np.tile((0x11111111, 0x22222222, 0x44444444, 0x88888888), 8)
        compiled[1, 32](arr)
        self.assertTrue(np.all(arr == exp))

    def test_activemask(self):
        @cuda.jit
        def use_activemask(x):
            i = cuda.grid(1)
            if (i % 2) == 0:
                # Even numbered threads fill in even numbered array entries
                # with binary "...01010101"
                x[i] = cuda.activemask()
            else:
                # Odd numbered threads fill in odd numbered array entries
                # with binary "...10101010"
                x[i] = cuda.activemask()

        out = np.zeros(32, dtype=np.uint32)
        use_activemask[1, 32](out)

        # 0x5 = 0101: The pattern from even-numbered threads
        # 0xA = 1010: The pattern from odd-numbered threads
        expected = np.tile((0x55555555, 0xAAAAAAAA), 16)
        np.testing.assert_equal(expected, out)

    def test_lanemask_lt(self):
        @cuda.jit
        def use_lanemask_lt(x):
            i = cuda.grid(1)
            x[i] = cuda.lanemask_lt()

        out = np.zeros(32, dtype=np.uint32)
        use_lanemask_lt[1, 32](out)

        # A string of 1s that grows from the LSB for each entry:
        # 0, 1, 3, 7, F, 1F, 3F, 7F, FF, 1FF, etc.
        # or in binary:
        # ...0001, ....0011, ...0111, etc.
        expected = np.asarray([(2 ** i) - 1 for i in range(32)],
                              dtype=np.uint32)
        np.testing.assert_equal(expected, out)


if __name__ == '__main__':
    unittest.main()
