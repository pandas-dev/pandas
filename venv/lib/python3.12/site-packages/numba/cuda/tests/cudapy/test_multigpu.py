from numba import cuda
import numpy as np
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import threading
import unittest


class TestMultiGPUContext(CUDATestCase):
    @unittest.skipIf(len(cuda.gpus) < 2, "need more than 1 gpus")
    def test_multigpu_context(self):
        @cuda.jit("void(float64[:], float64[:])")
        def copy_plus_1(inp, out):
            i = cuda.grid(1)
            if i < out.size:
                out[i] = inp[i] + 1

        def check(inp, out):
            np.testing.assert_equal(inp + 1, out)

        N = 32
        A = np.arange(N, dtype=np.float64)
        B = np.arange(N, dtype=np.float64)

        with cuda.gpus[0]:
            copy_plus_1[1, N](A, B)

        check(A, B)

        copy_plus_1[1, N](A, B)
        check(A, B)

        with cuda.gpus[0]:
            A0 = np.arange(N, dtype=np.float64)
            B0 = np.arange(N, dtype=np.float64)
            copy_plus_1[1, N](A0, B0)

            with cuda.gpus[1]:
                A1 = np.arange(N, dtype=np.float64)
                B1 = np.arange(N, dtype=np.float64)
                copy_plus_1[1, N](A1, B1)

        check(A0, B0)
        check(A1, B1)

        A = np.arange(N, dtype=np.float64)
        B = np.arange(N, dtype=np.float64)
        copy_plus_1[1, N](A, B)
        check(A, B)

    @skip_on_cudasim('Simulator does not support multiple threads')
    def test_multithreaded(self):
        def work(gpu, dA, results, ridx):
            try:
                with gpu:
                    arr = dA.copy_to_host()

            except Exception as e:
                results[ridx] = e

            else:
                results[ridx] = np.all(arr == np.arange(10))

        dA = cuda.to_device(np.arange(10))

        nthreads = 10
        results = [None] * nthreads
        threads = [threading.Thread(target=work, args=(cuda.gpus.current,
                                                       dA, results, i))
                   for i in range(nthreads)]
        for th in threads:
            th.start()

        for th in threads:
            th.join()

        for r in results:
            if isinstance(r, BaseException):
                raise r
            else:
                self.assertTrue(r)

    @unittest.skipIf(len(cuda.gpus) < 2, "need more than 1 gpus")
    def test_with_context(self):

        @cuda.jit
        def vector_add_scalar(arr, val):
            i = cuda.grid(1)
            if i < arr.size:
                arr[i] += val

        hostarr = np.arange(10, dtype=np.float32)
        with cuda.gpus[0]:
            arr1 = cuda.to_device(hostarr)

        with cuda.gpus[1]:
            arr2 = cuda.to_device(hostarr)

        with cuda.gpus[0]:
            vector_add_scalar[1, 10](arr1, 1)

        with cuda.gpus[1]:
            vector_add_scalar[1, 10](arr2, 2)

        with cuda.gpus[0]:
            np.testing.assert_equal(arr1.copy_to_host(), (hostarr + 1))

        with cuda.gpus[1]:
            np.testing.assert_equal(arr2.copy_to_host(), (hostarr + 2))

    @unittest.skipIf(len(cuda.gpus) < 2, "need more than 1 gpus")
    def test_with_context_peer_copy(self):
        # Peer access is not always possible - for example, with one GPU in TCC
        # mode and one in WDDM - if that is the case, this test would fail so
        # we need to skip it.
        with cuda.gpus[0]:
            ctx = cuda.current_context()
            if not ctx.can_access_peer(1):
                self.skipTest('Peer access between GPUs disabled')

        # 1. Create a range in an array
        hostarr = np.arange(10, dtype=np.float32)

        # 2. Copy range array from host -> GPU 0
        with cuda.gpus[0]:
            arr1 = cuda.to_device(hostarr)

        # 3. Initialize a zero-filled array on GPU 1
        with cuda.gpus[1]:
            arr2 = cuda.to_device(np.zeros_like(hostarr))

        with cuda.gpus[0]:
            # 4. Copy range from GPU 0 -> GPU 1
            arr2.copy_to_device(arr1)

            # 5. Copy range from GPU 1 -> host and check contents
            np.testing.assert_equal(arr2.copy_to_host(), hostarr)


if __name__ == '__main__':
    unittest.main()
